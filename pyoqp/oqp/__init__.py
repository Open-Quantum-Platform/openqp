"""OQP instance"""

import os
from oqp.utils.mpi_utils import MPIManager
from oqp.runtime import resolve_oqp_root
MPIManager()
# we must import dftd4 ffi lib before oqp to load library correctly
try:
    import dftd4.interface
except ModuleNotFoundError:
    print('\nPyOQP: dftd4 is not available')


try:
    int(os.environ['OMP_NUM_THREADS'])
except (KeyError, ValueError):
    os.environ['OMP_NUM_THREADS'] = '1'


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
