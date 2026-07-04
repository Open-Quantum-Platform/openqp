"""OQP instance"""

import os
from oqp.utils.mpi_utils import MPIManager
from oqp.runtime import resolve_oqp_root
MPIManager()
# DFT-D4 is linked natively into liboqp (source/dftd4_interface.F90); the
# external `dftd4` Python package is no longer required or imported.


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

def _missing_symbol_stub(name):
    """Placeholder for a native symbol declared in include/oqp.h but absent from
    this liboqp build -- i.e. an optional feature module compiled out via
    ENABLE_NAMD / ENABLE_QMMM = OFF. Importing ``oqp`` and running energy-only
    jobs still work; only *using* the missing feature raises, and with an
    actionable message instead of a bare cffi symbol-not-found error."""

    def _stub(*args, **kwargs):
        raise RuntimeError(
            f"OpenQP native symbol '{name}' is not present in this liboqp build. "
            f"It belongs to an optional feature module that was compiled out "
            f"(e.g. -DENABLE_NAMD=OFF or -DENABLE_QMMM=OFF). Rebuild OpenQP with "
            f"that feature enabled to use it."
        )

    return _stub


for attr_name in dir(lib):
    try:
        attr_value = getattr(lib, attr_name)
    except AttributeError:
        # Declared in oqp.h but not exported by this liboqp (feature compiled
        # out via ENABLE_NAMD / ENABLE_QMMM = OFF). Bind a stub so `import oqp`
        # and energy-only jobs still succeed; the stub raises a clear error only
        # if the missing feature is actually used.
        globals()[attr_name] = _missing_symbol_stub(attr_name)
        continue
    if callable(attr_value):
        if attr_name not in ('oqp_init', 'oqp_clean', 'oqp_set_atoms'):
            globals()[attr_name] = _oqp_wrapper(attr_value)
        else:
            globals()[attr_name] = attr_value
