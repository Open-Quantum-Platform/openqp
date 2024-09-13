import platform
import sys
import os
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
            except Exception as e:
                print(f"Failed to finalize MPI: {str(e)}")

    def set_mpi_comm(self, data):
        if self.use_mpi:
            data["usempi"] = 1
            data["debug_mode"] = 0
            comm_handle = self.comm.py2f()
            data["comm"] = comm_handle
        else:
            data["usempi"] = 0

    def bcast(self, data, root=0):
        if self.use_mpi :
            data = self.comm.bcast(data, root=root)
        return data
