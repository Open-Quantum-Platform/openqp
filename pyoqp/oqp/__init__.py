"""OQP instance"""

import os
import platform
from oqp.utils.mpi_utils import MPIManager
MPIManager()
# we must import dftd4 ffi lib before oqp to load library correctly
try:
    import dftd4.interface
except ModuleNotFoundError:
    print('\nPyOQP: dftd4 is not available')

def _library_suffix():
    if platform.uname()[0] == "Windows":
        return "dll"
    if platform.uname()[0] == "Darwin":
        return "dylib"
    return "so"


def _is_oqp_root(path, suffix):
    return (
        path
        and os.path.exists(os.path.join(path, "include", "oqp.h"))
        and os.path.exists(os.path.join(path, "lib", f"liboqp.{suffix}"))
    )


def _resolve_oqp_root():
    """Choose a coherent root containing both oqp.h and liboqp.

    Installed wheels put the Python package, header, and native library under the
    package directory. Prefer that self-contained root over OPENQP_ROOT so a
    leftover development environment variable cannot mix an installed Python
    package with a different source-tree library. Source-tree imports do not have
    package-local include/lib directories, so they still use OPENQP_ROOT.
    """
    suffix = _library_suffix()
    package_root = os.path.abspath(os.path.dirname(__file__))
    env_root = os.environ.get("OPENQP_ROOT")

    if _is_oqp_root(package_root, suffix):
        os.environ["OPENQP_ROOT"] = package_root
        return package_root, suffix
    if _is_oqp_root(env_root, suffix):
        return env_root, suffix

    if env_root:
        raise RuntimeError(
            "OPENQP_ROOT does not contain matching include/oqp.h and "
            f"lib/liboqp.{suffix}: {env_root}"
        )
    raise RuntimeError(
        "Cannot locate OpenQP runtime files. Install OpenQP as a package or set "
        "OPENQP_ROOT to a tree containing include/oqp.h and lib/liboqp."
    )


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

if RTLD:
    from cffi import FFI

    ffi = FFI()
    oqp_root, suffix = _resolve_oqp_root()

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
