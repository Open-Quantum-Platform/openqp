import os
import sys
from collections.abc import Iterator

try:
    from mpi4py import MPI
    mpi4py_available = True
except (ImportError, OSError) as e:
    mpi4py_available = False

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
                if cls._instance.rank != 0:
                    sys.stdout = open(os.devnull, 'w')
            else:
                print(f"Failed to import mpi4py")
                cls._instance = super(MPIManager, cls).__new__(cls)
                cls._instance.rank = 0
                cls._instance.size = 1
                cls._instance.use_mpi = 0
                cls._instance.mpi4py_available = False
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
        if self.use_mpi:
            data["usempi"] = 1
            data["debug_mode"] = 0
            comm_handle = self.comm.py2f()
            data["comm"] = comm_handle
        else:
            data["usempi"] = 0

    def bcast(self, data, root=0, barrier=True):
        if self.use_mpi:
            data = self.comm.bcast(data, root=root)
            if barrier:
                self.comm.Barrier()
        return data


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

def mpi_write(func):
    # mpi decorator to write files
    def wrapper(self, *args):
        if self.mpi_manager.rank == 0 or not self.usempi:
            func(self, *args)
    return wrapper

def mpi_dump(func):
    # mpi decorator to write logs
    def wrapper(mol, *args, **kwargs):
        if MPIManager().rank == 0 or not mol.usempi:
            func(mol, *args, **kwargs)
    return wrapper
