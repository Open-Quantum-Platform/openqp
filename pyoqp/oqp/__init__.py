"""OQP instance"""

import os
from oqp.utils.mpi_utils import MPIManager
from oqp.runtime import library_path, resolve_oqp_root
MPIManager()
# DFT-D4 is linked dynamically into liboqp (source/dftd4_interface.F90); the
# package-local shared libraries replace the external Python `dftd4` package.


#: True when OMP_NUM_THREADS carried a usable value before OpenQP was imported,
#: i.e. the thread count is the caller's request rather than the default below.
#: Kernels that take an explicit thread argument (the native FCI/CASSCF/RDM
#: entry points, see oqp.library.fci._fci_lib_threads) need to tell the two
#: apart: a defaulted '1' means "nobody asked", not "run me serially".
#: pyoqp.Runner sets this to True when it applies `[input] omp_threads`.
OMP_THREADS_FROM_ENV = True
try:
    int(os.environ['OMP_NUM_THREADS'])
except (KeyError, ValueError):
    os.environ['OMP_NUM_THREADS'] = '1'
    OMP_THREADS_FROM_ENV = False


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
    lib = ffi.dlopen(str(library_path(oqp_root, suffix)))

else:
    from _oqp import ffi, lib

for attr_name in dir(lib):
    try:
        attr_value = getattr(lib, attr_name)
    except AttributeError:
        # include/oqp.h declares every C entry point, but a build may leave an
        # optional backend out (e.g. ENABLE_DFTD4=OFF, the default on Windows,
        # drops oqp_dftd4_disp/_v2).  In cffi's ABI mode the declaration exists
        # while the symbol does not, so skip it rather than failing the import;
        # callers already probe with hasattr(oqp.lib, ...).
        continue
    if callable(attr_value):
        if attr_name not in (
            'oqp_init', 'oqp_clean', 'oqp_set_atoms',
            'oqp_namd_counter_random', 'oqp_namd_counter_normal_fill',
            'oqp_namd_baeck_an_tdc', 'oqp_namd_nacme_gate',
            'oqp_odp_umbrella_evaluate',
            'oqp_namd_droplet_boundary', 'oqp_namd_com_restraint',
            'oqp_namd_langevin_thermostat',
            'oqp_maximum_overlap_assignment', 'oqp_diagonal_phase_tracking',
            'oqp_simplex_qp_solve', 'oqp_simplex_qp_solve_avoid',
        ):
            globals()[attr_name] = _oqp_wrapper(attr_value)
        else:
            globals()[attr_name] = attr_value

#: True once the compiled liboqp backend loaded successfully (the import above
#: fails loudly otherwise). WF/test backend guards probe this flag.
BACKEND_AVAILABLE = True
