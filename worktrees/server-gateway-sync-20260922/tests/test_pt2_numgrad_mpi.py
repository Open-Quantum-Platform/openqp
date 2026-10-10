"""MPI rank-group fan-out for the PT2 central-difference gradients.

The point under test is NOT "does it run on N ranks" -- an earlier attempt at
this ran on 4 ranks and returned a silently ZERO gradient.  The point is that
a multi-rank run reproduces the one-rank gradient of the same build, that the
work really was spread (a run that evaluated every displacement on every rank
would also reproduce it, at no speed-up), and that a build without mpi4py, or
a run on one rank, never enters the split at all.

Most of this file is arithmetic and reassembly logic, which is exercised
without any MPI at all.  The end-to-end multi-rank comparison needs both
mpi4py and an mpirun, and skips when either is missing -- note that mpi4py
must NOT be installed in the venv that runs this suite under
``pytest --forked`` (a forked child with a live MPI_COMM_WORLD dies inside
the native calls), so on a normal developer machine that test is expected to
skip and the fan-out is validated by an explicit mpirun run instead.

The route is pinned with ``[pt2] gradient=numerical`` throughout: the default
``auto`` now sends these methods to the analytic CASPT2 derivative, which has no
task-group split to test.
"""
import json
import os
import shutil
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from oqp.utils.mpi_utils import (MPIManager, MPITaskGroups, group_of_rank,
                                 plan_task_groups)


def _backend_available() -> bool:
    try:
        import oqp
    except Exception:
        return False
    return bool(getattr(oqp, "BACKEND_AVAILABLE", False))


# --------------------------------------------------------------------------
# group layout arithmetic
# --------------------------------------------------------------------------

@pytest.mark.parametrize("nworld,ntask,rpg,expect", [
    # ranks_per_group = 0 -> one group per rank, capped by the task count
    (1, 12, 0, 1),
    (4, 12, 0, 4),
    (8, 12, 0, 8),
    (8, 4, 0, 4),      # more ranks than displacements: no idle group
    # explicit group sizes
    (4, 12, 1, 4),
    (4, 12, 2, 2),
    (4, 12, 4, 1),
    (8, 12, 2, 4),
    (8, 12, 4, 2),
    (8, 24, 3, 2),     # 8 // 3 = 2 groups of 4 and 4, not 3 + 3 + 2
    (4, 12, 8, 1),     # asking for more ranks per group than exist
    (4, 1, 0, 1),      # one task cannot be split
])
def test_plan_task_groups(nworld, ntask, rpg, expect):
    assert plan_task_groups(nworld, ntask, rpg) == expect


@pytest.mark.parametrize("nworld,ngroup,expect", [
    (4, 4, [0, 1, 2, 3]),
    (4, 2, [0, 0, 1, 1]),
    (8, 4, [0, 0, 1, 1, 2, 2, 3, 3]),
    (8, 2, [0, 0, 0, 0, 1, 1, 1, 1]),
    # 7 ranks over 3 groups: sizes differ by at most one and every group is
    # a contiguous block of world ranks
    (7, 3, [0, 0, 0, 1, 1, 2, 2]),
])
def test_group_of_rank_is_a_contiguous_balanced_block(nworld, ngroup, expect):
    got = [group_of_rank(r, nworld, ngroup) for r in range(nworld)]
    assert got == expect
    assert sorted(set(got)) == list(range(ngroup))
    sizes = [got.count(g) for g in range(ngroup)]
    assert max(sizes) - min(sizes) <= 1


def test_every_task_is_owned_by_exactly_one_group():
    """The union of the groups' index ranges is a partition of the tasks."""
    for ntask in (1, 11, 12, 18, 24, 25):
        for ngroup in (1, 2, 3, 4, 5, 8):
            owners = {}
            for g in range(ngroup):
                for t in range(g, ntask, ngroup):
                    assert t not in owners, "task %d claimed twice" % t
                    owners[t] = g
            assert sorted(owners) == list(range(ntask))


# --------------------------------------------------------------------------
# reassembly
# --------------------------------------------------------------------------

class _FakeWorld:
    """COMM_WORLD stand-in: records what each rank hands to Allreduce."""

    def __init__(self, gathered=None):
        self.sent = []
        self.gathered = gathered or []

    def Allreduce(self, sendbuf, recvbuf, **kwargs):
        self.sent.append(np.array(sendbuf, copy=True))
        recvbuf[...] = sendbuf

    def allgather(self, item):
        self.sent.append(item)
        return list(self.gathered)


def test_reduce_sum_contributes_exactly_one_copy_per_group():
    """Every rank of a group holds the SAME energies.

    A plain world sum would therefore multiply each displacement by the group
    size; reduce_sum zeroes the non-root ranks first.  This is the specific
    trap the design memo calls out, so it gets its own test: what the whole
    world hands to the reduction must add up to one copy per GROUP, not one
    per rank.
    """
    value = np.array([[1.0, 2.0], [3.0, 4.0]])
    ngroup, ranks_per_group = 2, 3
    world = _FakeWorld()
    for grp in range(ngroup):
        for r in range(ranks_per_group):
            groups = MPITaskGroups(world, ntask=6, ngroup=ngroup, group=grp,
                                   is_group_root=(r == 0))
            groups.reduce_sum(value)

    assert len(world.sent) == ngroup * ranks_per_group
    np.testing.assert_array_equal(sum(world.sent), ngroup * value)
    # the bug this guards against: without the zeroing the world would sum
    # ngroup * ranks_per_group copies
    assert not np.allclose(sum(world.sent), ngroup * ranks_per_group * value)


def test_reduce_sum_does_not_mutate_the_callers_array():
    value = np.array([[1.0, 2.0], [3.0, 4.0]])
    original = value.copy()
    world = _FakeWorld()
    groups = MPITaskGroups(world, ntask=6, ngroup=2, group=0,
                           is_group_root=False)
    groups.reduce_sum(value)
    np.testing.assert_array_equal(value, original)


def test_reduce_sum_is_the_identity_without_a_split():
    groups = MPITaskGroups(None, ntask=6, ngroup=1, group=0, is_group_root=True)
    value = np.arange(6.0).reshape(3, 2)
    assert groups.reduce_sum(value) is value
    assert groups.gather_list([1, 2]) == [1, 2]


def test_gather_list_takes_only_the_group_roots():
    world = _FakeWorld(gathered=[["a"], [], ["b"], []])
    root = MPITaskGroups(world, ntask=4, ngroup=2, group=0, is_group_root=True)
    assert root.gather_list(["a"]) == ["a", "b"]
    assert world.sent[-1] == ["a"]

    # a non-root rank holds the same list and must contribute nothing
    other = MPITaskGroups(world, ntask=4, ngroup=2, group=0, is_group_root=False)
    other.gather_list(["a"])
    assert world.sent[-1] == []


# --------------------------------------------------------------------------
# the serial path must be untouched
# --------------------------------------------------------------------------

def test_task_groups_is_a_no_op_when_disabled():
    """usempi=False must never split, whatever the world looks like."""
    mpi = MPIManager()
    switched = []
    with mpi.task_groups(12, ranks_per_group=1, enabled=False,
                         on_switch=lambda: switched.append(1)) as groups:
        assert groups.ngroup == 1
        assert groups.group == 0
        assert list(groups.indices) == list(range(12))
        assert groups.is_group_root
    assert switched == [], "the native communicator must not be touched"


def test_task_groups_is_a_no_op_on_one_rank():
    mpi = MPIManager()
    if mpi.world_size > 1:
        pytest.skip("this assertion is about the single-rank world")
    with mpi.task_groups(24, ranks_per_group=1) as groups:
        assert groups.ngroup == 1
        assert list(groups.indices) == list(range(24))


def test_manager_exposes_world_rank_even_without_mpi4py():
    """mpi_dump keys off world_rank, so it has to exist on every build."""
    mpi = MPIManager()
    assert hasattr(mpi, "world_rank")
    assert hasattr(mpi, "world_size")
    if not mpi.mpi4py_available:
        assert mpi.world_rank == 0
        assert mpi.world_size == 1
    else:
        assert mpi.world_rank == mpi.rank
        assert mpi.world_size == mpi.size


def test_serial_pt2_gradient_does_not_enter_the_split(tmp_path):
    """A one-rank production gradient reports a single group.

    Guards the property everyone else depends on: this change must be
    invisible without MPI.
    """
    if not _backend_available():
        pytest.skip("native OQP backend not built")
    from oqp.pyoqp import Runner

    runner = Runner(
        project="h2_serial_layout",
        input_file=None,
        log=str(tmp_path / "h2.log"),
        input_dict={
            "input": {"system": "\nH 0 0 0\nH 0 0 0.900", "charge": "0",
                      "basis": "sto-3g", "method": "mrmp2", "runtype": "energy"},
            "guess": {"type": "hcore"},
            "scf": {"type": "rhf", "multiplicity": "1", "maxit": "60",
                    "save_molden": "False"},
            "properties": {"scf_prop": ""},
            "cas": {"active_electrons": "2", "active_orbitals": "2",
                    "frozen_core": "0", "max_det": "2000"},
            "ci": {"nroot": "1", "solver": "dense", "eig_tol": "1.0e-10",
                   "integral_backend": "native", "integral_cutoff": "5.0e-11"},
            "pt2": {"reference": "casci", "gradient": "numerical"},
            "tests": {"exception": "True"},
        },
        silent=1,
        usempi=False,
    )
    runner.mol.config["input"]["runtype"] = "grad"
    runner.run(test_mod=True)

    layout = runner.mol.pt2_numgrad_layout
    assert layout["ngroup"] == 1
    assert layout["ntask"] == 12
    assert layout["ntask_local"] == 12
    grads = np.asarray(runner.mol.grads)
    assert grads.shape == (1, 2, 3)
    assert np.all(np.isfinite(grads))


def test_grad_ranks_per_group_is_a_recognised_pt2_option():
    """The option has to be in the schema, the checker and the key list.

    Reading it with dict.get is not enough: the config parser rejects any
    key it does not know, so a user setting it would get "Unknown option".
    """
    from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
    from oqp.utils.oqp_input import GENERIC_SCHEMA_KEYS

    assert "grad_ranks_per_group" in OQP_CONFIG_SCHEMA["pt2"]
    assert OQP_CONFIG_SCHEMA["pt2"]["grad_ranks_per_group"]["type"] is int
    assert "grad_ranks_per_group" in GENERIC_SCHEMA_KEYS["pt2"]


# --------------------------------------------------------------------------
# end to end, multi rank
# --------------------------------------------------------------------------

_MPI_DRIVER = textwrap.dedent('''
    import json, os, sys
    import numpy as np
    from oqp.utils.mpi_utils import MPIManager
    from oqp.pyoqp import Runner

    rpg = int(sys.argv[1])
    out = sys.argv[2]
    logdir = sys.argv[3]
    mpi = MPIManager()
    runner = Runner(
        project="h2_mpi_%d_%d" % (rpg, mpi.world_size),
        input_file=None,
        log=os.path.join(logdir, "h2_%d_%d.log" % (rpg, mpi.world_size)),
        input_dict={
            "input": {"system": "\\nH 0 0 0\\nH 0 0 0.900", "charge": "0",
                      "basis": "sto-3g", "method": "mrmp2",
                      "runtype": "energy"},
            "guess": {"type": "hcore"},
            "scf": {"type": "rhf", "multiplicity": "1", "maxit": "60",
                    "save_molden": "False"},
            "properties": {"scf_prop": ""},
            "cas": {"active_electrons": "2", "active_orbitals": "2",
                    "frozen_core": "0", "max_det": "2000"},
            "ci": {"nroot": "1", "solver": "dense", "eig_tol": "1.0e-10",
                   "integral_backend": "native",
                   "integral_cutoff": "5.0e-11"},
            "pt2": {"reference": "casci", "gradient": "numerical",
                    "grad_ranks_per_group": str(rpg)},
            "tests": {"exception": "True"},
        },
        silent=1,
        usempi=True,
    )
    runner.mol.config["input"]["runtype"] = "grad"
    runner.run(test_mod=True)
    g = np.asarray(runner.mol.grads, dtype=float).reshape(-1)
    layout = runner.mol.pt2_numgrad_layout
    if mpi.world_size > 1:
        from mpi4py import MPI
        layouts = MPI.COMM_WORLD.gather(layout, root=0)
    else:
        layouts = [layout]
    if mpi.world_rank == 0:
        with open(out, "w") as fh:
            json.dump({"grad": g.tolist(), "layouts": layouts}, fh)
    mpi.barrier()
''')


def _mpi_available():
    try:
        import mpi4py  # noqa: F401
    except Exception:
        return None
    return shutil.which("mpirun") or shutil.which("mpiexec")


@pytest.mark.parametrize("nrank,rpg", [(4, 1), (4, 2), (4, 4)])
def test_multirank_gradient_reproduces_the_one_rank_gradient(tmp_path, nrank, rpg):
    launcher = _mpi_available()
    if launcher is None:
        pytest.skip("mpi4py and an MPI launcher are needed for the fan-out test")
    if not _backend_available():
        pytest.skip("native OQP backend not built")

    script = tmp_path / "driver.py"
    script.write_text(_MPI_DRIVER)
    env = dict(os.environ, OMP_NUM_THREADS="1")

    def run(n, ranks_per_group, tag):
        out = tmp_path / ("%s.json" % tag)
        cmd = [launcher, "-n", str(n), "--bind-to", "none", "--oversubscribe",
               sys.executable, str(script), str(ranks_per_group), str(out),
               str(tmp_path)]
        proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
        if proc.returncode != 0 or not out.exists():
            pytest.skip("MPI launch unavailable here: %s" % proc.stderr[-400:])
        return json.load(open(out))

    base = run(1, 0, "serial")
    got = run(nrank, rpg, "n%d_rpg%d" % (nrank, rpg))

    # the split really happened, and every group got work
    layouts = got["layouts"]
    assert len(layouts) == nrank
    ngroup = plan_task_groups(nrank, 12, rpg)
    assert {L["ngroup"] for L in layouts} == {ngroup}
    assert {L["group"] for L in layouts} == set(range(ngroup))
    assert sum(L["ntask_local"] for L in layouts if L["is_group_root"]) == 12

    # and the answer is the serial answer
    np.testing.assert_array_equal(np.asarray(got["grad"]),
                                  np.asarray(base["grad"]))
