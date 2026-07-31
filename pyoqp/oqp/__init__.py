"""OQP instance"""

import os
import re
import sys


def _configure_threads_before_native_import():
    """Configure OpenMP before importing anything that can load liboqp.

    A console-script entry point imports :mod:`oqp` before it imports
    ``oqp.pyoqp``.  Consequently, handling ``--omp`` in ``pyoqp.main`` is too
    late for OpenMP runtimes that cache ``OMP_NUM_THREADS`` when liboqp is
    loaded.  Keep this small pre-import parser here, ahead of MPI, DFT-D4,
    NumPy/BLAS, and CFFI.
    """
    requested = None
    for index, argument in enumerate(sys.argv[1:], start=1):
        if argument == '--omp' and index + 1 < len(sys.argv):
            requested = sys.argv[index + 1]
            break
        if argument.startswith('--omp='):
            requested = argument.split('=', 1)[1]
            break

    if requested is None:
        input_file = next((argument for argument in sys.argv[1:]
                           if not argument.startswith('-')
                           and os.path.isfile(argument)), None)
        if input_file:
            try:
                with open(input_file, encoding='utf-8', errors='ignore') as handle:
                    match = re.search(
                        r'(?mi)^[ \t]*omp_threads[ \t]*=[ \t]*(\d+)',
                        handle.read(),
                    )
                if match:
                    requested = match.group(1)
            except OSError:
                pass

    try:
        requested_threads = int(requested) if requested is not None else None
    except (TypeError, ValueError):
        requested_threads = None
    if requested_threads is not None and requested_threads > 0:
        os.environ['OMP_NUM_THREADS'] = str(requested_threads)

    # OpenQP owns the outer parallel regions; threaded BLAS inside those
    # regions would oversubscribe the machine.
    for variable in (
        'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS',
    ):
        os.environ.setdefault(variable, '1')
    os.environ.setdefault('OMP_STACKSIZE', '256M')
    os.environ.setdefault('GOMP_STACKSIZE', '256M')

    try:
        if int(os.environ.get('OMP_NUM_THREADS', '0')) < 1:
            raise ValueError
    except ValueError:
        os.environ['OMP_NUM_THREADS'] = '1'


_configure_threads_before_native_import()

from oqp.utils.mpi_utils import MPIManager
from oqp.runtime import resolve_oqp_root
MPIManager()
# we must import dftd4 ffi lib before oqp to load library correctly
try:
    import dftd4.interface
except ModuleNotFoundError:
    print('\nPyOQP: dftd4 is not available')

def _oqp_wrapper(func):
    """Decorator for OQP library functions"""

    def wrapper(molecule, *args):
        return func(molecule.data._data, *args)

    return wrapper


if os.environ.get('OQP_RTLD'):
    RTLD = str(os.environ.get('OQP_RTLD')).lower() in ('true', '1', 't', 'y', 'yes', 'on')
else:
    RTLD = True

oqp_root, suffix = resolve_oqp_root()

if RTLD:
    from cffi import FFI

    ffi = FFI()

    with open(f"{oqp_root}/include/oqp.h", "r", encoding="ascii") as oqp_header:
        defs = oqp_header.read().replace("#include", "//#include")

    ffi.cdef(defs)
    lib = ffi.dlopen(f"{oqp_root}/lib/liboqp.{suffix}")

else:
    from _oqp import ffi, lib

for attr_name in dir(lib):
    attr_value = getattr(lib, attr_name)
    if callable(attr_value):
        if attr_name not in ('oqp_init', 'oqp_clean', 'oqp_set_atoms'):
            globals()[attr_name] = _oqp_wrapper(attr_value)
        else:
            globals()[attr_name] = attr_value
