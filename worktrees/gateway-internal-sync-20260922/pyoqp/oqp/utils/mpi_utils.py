import contextlib
import os
import sys
from collections.abc import Iterator

try:
    from mpi4py import MPI
    mpi4py_available = True
except (ImportError, OSError) as e:
    mpi4py_available = False

def plan_task_groups(nworld, ntask, ranks_per_group=0):
    """How many rank groups ``nworld`` ranks should form for ``ntask`` tasks.

    ``ranks_per_group <= 0`` means "as much fan-out as there is work for":
    one group per task, capped by the world size.  A positive value is the
    requested group size, and the group count is capped so that no group is
    left without a task (a group with nothing to do would be pure idle,
    while a slightly larger group still contributes to its own SCF).
    """
    nworld = max(int(nworld), 1)
    ntask = max(int(ntask), 0)
    if ranks_per_group and int(ranks_per_group) > 0:
        ngroup = max(1, nworld // int(ranks_per_group))
    else:
        ngroup = nworld
    return max(1, min(ngroup, max(ntask, 1), nworld))


def group_of_rank(rank, nworld, ngroup):
    """Block (not round-robin) assignment of world ranks to groups.

    Ranks of one group are then adjacent in COMM_WORLD, which is where a
    launcher puts ranks that share a socket, and group sizes differ by at
    most one when ``ngroup`` does not divide ``nworld``.
    """
    return (int(rank) * int(ngroup)) // int(nworld)


class MPIManager:

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            if mpi4py_available:
                cls._instance = super(MPIManager, cls).__new__(cls)
                cls._instance.comm = MPI.COMM_WORLD
                cls._instance.rank = cls._instance.comm.Get_rank()
                cls._instance.size = cls._instance.comm.Get_size()
                cls._instance.use_mpi = int(cls._instance.size > 1)
                cls._instance.mpi4py_available = True
                # world_* never follow a sub-communicator swap (see
                # task_groups).  Everything with a side effect outside this
                # process -- log appends, molden/json writes, scratch cleanup
                # -- must be keyed off world_rank, because inside a split
                # EVERY group root has rank == 0 and they would all write the
                # same file at the same time.
                cls._instance.world_comm = MPI.COMM_WORLD
                cls._instance.world_rank = cls._instance.rank
                cls._instance.world_size = cls._instance.size
                if cls._instance.rank != 0:
                    sys.stdout = open(os.devnull, 'w')
            else:
                print(f"Failed to import mpi4py")
                cls._instance = super(MPIManager, cls).__new__(cls)
                cls._instance.rank = 0
                cls._instance.size = 1
                cls._instance.use_mpi = 0
                cls._instance.mpi4py_available = False
                cls._instance.world_comm = None
                cls._instance.world_rank = 0
                cls._instance.world_size = 1
        return cls._instance

    def mpi4py_available(self):
        return self.mpi4py_available

    def barrier(self):

        if self.use_mpi:
            self.comm.Barrier()

    def finalize_mpi(self):
        if self.use_mpi:
            try:
                self.comm.Barrier()
                MPI.Finalize()
            except Exception as err:
                print(f"Failed to finalize MPI: {str(err)}")

    def set_mpi_comm(self, data):
        """Publish the CURRENT communicator to the native control struct.

        ``data["comm"]`` is ``infos%mpiinfo%comm`` and every native module
        re-reads it through ``pe%init`` on entry, so calling this again after
        ``self.comm`` changes is what makes a sub-communicator take effect
        inside liboqp.  The handle is written even when MPI is off, so the
        struct can never be left holding a communicator that has since been
        freed.
        """
        if self.use_mpi:
            data["usempi"] = 1
            data["debug_mode"] = 0
            comm_handle = self.comm.py2f()
            data["comm"] = comm_handle
        else:
            data["usempi"] = 0
            if self.mpi4py_available:
                data["comm"] = self.comm.py2f()

    @contextlib.contextmanager
    def task_groups(self, ntask, ranks_per_group=0, enabled=True, on_switch=None):
        """Split COMM_WORLD into groups, one per concurrent task.

        Existing OQP MPI parallelism lives INSIDE one energy evaluation: the
        native layer receives a communicator and decomposes the SCF, the
        integrals and the grid collectively over it.  A caller that owns a
        list of independent whole-pipeline evaluations therefore cannot hand
        different work to different ranks of COMM_WORLD -- doing that leaves
        each rank waiting in a collective its peers never enter.

        This splits COMM_WORLD into ``ngroup`` sub-communicators and installs
        the sub-communicator as ``self.comm`` for the duration of the block,
        so that both the Python-level collectives (bcast, mpi_get_attr, the
        basis-set broadcast) and the native layer become collective over the
        GROUP.  Each group then runs whole tasks by itself and MPI inside one
        evaluation keeps working, one level down.

        ``on_switch`` is invoked (with no arguments) right after the swap and
        right after the restore.  Callers holding a native control struct
        pass ``lambda: manager.set_mpi_comm(mol.data)`` so liboqp sees the
        same communicator the Python layer does.

        The restore is unconditional, so a task that raises still leaves
        COMM_WORLD installed -- otherwise the next unrelated collective
        would hang instead of failing.  A failure flag is reduced over
        COMM_WORLD on the way out, so a group that failed cannot leave the
        other groups blocked in the reassembly reduction.

        Parameters
        ----------
        ntask : int
            Number of independent tasks to be spread over the groups.
        ranks_per_group : int
            Ranks per group; 0 (default) picks ``ngroup = min(world_size,
            ntask)``, i.e. as much fan-out as there is work for.  Groups need
            not be equal in size and ``ntask`` need not divide by ``ngroup``.
        enabled : bool
            False forces the single-group (unsplit) path.  Used by callers
            that must honour a per-molecule ``usempi=False``.
        """
        ntask = max(int(ntask), 0)
        world = getattr(self, 'world_comm', None)
        if (not self.mpi4py_available or world is None
                or not enabled or self.world_size <= 1 or ntask <= 1):
            yield MPITaskGroups(None, ntask, 1, 0, True)
            return

        nworld = self.world_size
        ngroup = plan_task_groups(nworld, ntask, ranks_per_group)

        if ngroup == 1:
            yield MPITaskGroups(world, ntask, 1, 0, self.world_rank == 0)
            return

        group = group_of_rank(self.world_rank, nworld, ngroup)
        subcomm = world.Split(color=group, key=self.world_rank)

        saved = (self.comm, self.rank, self.size, self.use_mpi)
        self.comm = subcomm
        self.rank = subcomm.Get_rank()
        self.size = subcomm.Get_size()
        self.use_mpi = int(self.size > 1)
        if on_switch is not None:
            on_switch()

        groups = MPITaskGroups(world, ntask, ngroup, group, self.rank == 0)
        failed = 0
        try:
            yield groups
        except BaseException:
            failed = 1
            raise
        finally:
            self.comm, self.rank, self.size, self.use_mpi = saved
            if on_switch is not None:
                on_switch()
            try:
                subcomm.Free()
            except Exception:  # pragma: no cover - communicator already gone
                pass
            # Rendezvous on COMM_WORLD before anyone leaves, so a group that
            # raised takes the others down with it instead of stranding them
            # in groups.reduce_sum.
            anyfail = world.allreduce(failed, op=MPI.MAX)
            if anyfail and not failed:
                raise RuntimeError(
                    'MPI task groups: another rank group raised while '
                    'evaluating its share of the tasks; this rank group is '
                    'aborting so the run fails instead of deadlocking. The '
                    'original traceback is on the rank group that failed.')

    def bcast(self, data, root=0, barrier=True):
        if self.use_mpi:
            data = self.comm.bcast(data, root=root)
            if barrier:
                self.comm.Barrier()
        return data

    def allreduce_sum(self, array):
        """Sum a numpy array across all ranks, in place on a copy.

        Used to collect work distributed by disjoint slot: every rank fills
        only the entries it owns and leaves the rest zero, so a plain sum
        reassembles the whole array without any explicit index bookkeeping.
        """
        if not self.use_mpi:
            return array
        import numpy as np

        contiguous = np.ascontiguousarray(array)
        total = np.zeros_like(contiguous)
        self.comm.Allreduce(contiguous, total, op=MPI.SUM)
        return total

    def allgather_list(self, items):
        """Concatenate per-rank lists into one list held by every rank."""
        if not self.use_mpi:
            return list(items)
        gathered = self.comm.allgather(list(items))
        return [item for chunk in gathered for item in chunk]

    def split_indices(self, count):
        """Round-robin the integers ``0..count-1`` over the ranks."""
        if not self.use_mpi:
            return range(count)
        return range(self.rank, count, self.size)


class MPITaskGroups:
    """The layout handed back by ``MPIManager.task_groups``.

    Holds COMM_WORLD, not the sub-communicator, so the reassembly below is
    still correct after the context manager has restored the world.
    """

    def __init__(self, world_comm, ntask, ngroup, group, is_group_root):
        self.world_comm = world_comm
        self.ntask = ntask
        self.ngroup = ngroup
        self.group = group
        self.is_group_root = bool(is_group_root)
        self.split = world_comm is not None and ngroup > 1

    @property
    def indices(self):
        """The task indices this rank's group owns.

        Round robin over groups, so a task count that does not divide by
        ``ngroup`` spreads the remainder one task per group rather than
        piling it on the last one.
        """
        return range(self.group, self.ntask, self.ngroup)

    def reduce_sum(self, array):
        """Reassemble a per-task array over COMM_WORLD.

        THE DOUBLE-COUNTING TRAP: every rank of a group runs the group's
        tasks collectively and ends up holding the SAME values, so a plain
        sum over COMM_WORLD would multiply each task by its group size.  The
        non-root ranks of each group are zeroed first, which makes the sum
        pick up exactly one copy per task -- and also makes the result
        independent of whether the group members really did agree.
        """
        if not self.split:
            return array
        import numpy as np

        # np.array(copy=True): ascontiguousarray can hand back the caller's
        # own buffer, and the non-root branch zeroes it in place.
        contiguous = np.array(array, dtype=float, order='C', copy=True)
        if not self.is_group_root:
            contiguous[...] = 0.0
        total = np.zeros_like(contiguous)
        kwargs = {'op': MPI.SUM} if mpi4py_available else {}
        self.world_comm.Allreduce(contiguous, total, **kwargs)
        return total

    def gather_list(self, items):
        """Concatenate the group roots' lists into one list held by everyone."""
        if not self.split:
            return list(items)
        mine = list(items) if self.is_group_root else []
        gathered = self.world_comm.allgather(mine)
        return [item for chunk in gathered for item in chunk]


class MPIPool:

    def __init__(self, processes):
        self.processes = processes

    def imap_unordered(self, func, inp):
        return MPIMap(func, inp, self.processes)

    def close(self):
        pass

class MPIMap(Iterator):

    def __init__(self, func, var, processes):
        self.a = -1
        self.func = func
        self.var = var
        self.processes = processes
        self.mpi_manager = MPIManager()

    def __next__(self):
        self.a += 1

        if self.a == len(self.var):
            raise StopIteration

        slot = self.a % self.processes

        if slot == self.mpi_manager.rank:
            data = self.func(self.var[self.a])

        else:
            data = None

        if self.mpi_manager.rank == 0:
            if slot != 0:
                data = self.mpi_manager.comm.irecv(source=slot).wait()

            return data
        else:
            if slot == self.mpi_manager.rank:
                self.mpi_manager.comm.isend(data, dest=0).wait()

def mpi_get_attr(func):
    # mpi decorator to get values then broadcast before return values
    def wrapper(self, *args):
        if self.mpi_manager.rank == 0 or not self.usempi:
            attr = func(self, *args)
        else:
            attr = None

        if self.usempi:
            attr = self.mpi_manager.bcast(attr)
        return attr
    return wrapper

def mpi_update_attr(func):
    # mpi decorator to update values then broadcast without return values
    def wrapper(self, *args):
        if self.mpi_manager.rank == 0 or not self.usempi:
            func(self, *args)

        if self.usempi:
            self.data = self.mpi_manager.bcast(self.data)
    return wrapper

def mpi_dump(func):
    # mpi decorator to write logs
    def wrapper(mol, *args, **kwargs):
        # world_rank, NOT rank: inside an MPIManager.task_groups split every
        # group root has rank == 0, and these functions all append to one
        # shared log / molden / json path.  world_rank keeps exactly one
        # writer.  Outside a split the two are the same number.
        if MPIManager().world_rank == 0 or not mol.usempi:
            func(mol, *args, **kwargs)

    return wrapper
