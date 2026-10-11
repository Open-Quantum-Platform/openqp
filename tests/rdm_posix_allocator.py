"""GNU/glibc Linux and macOS compiler-allocator regression lane.

Explicitly consumed by .github/workflows/CI.yml on GNU Linux. This filename
keeps a platform-specific interposer out of Windows' portable test collection;
invoke it explicitly on supported managed macOS builds too. Unsupported hosts
fail if asked to run it. No synthetic controller is armed in these tests.
"""

import ctypes
import os
import subprocess
import sys

import pytest

from test_rdm_fallback import native

@pytest.fixture
def allocation_interposer(native, tmp_path):
    """Fail each of the 18 owned allocations, including both link helpers.

    The interposer is confined to a child and fails malloc calls originating
    from liboqp only. It never changes the installed library or compiler flags.
    """
    assert sys.platform in ("linux", "darwin"), "POSIX allocator lane requires GNU/glibc Linux or macOS"
    source = tmp_path / "fail_rdm_alloc.c"
    library = tmp_path / ("fail_rdm_alloc.dylib" if sys.platform == "darwin" else "fail_rdm_alloc.so")
    source.write_text(r'''
#define _GNU_SOURCE
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>
#include <errno.h>
static int selected, count, active, strings_only;
static size_t total_bytes, largest_bytes;
void rdm_fail_allocation(int n) { selected=n; count=0; active=1; strings_only=0; total_bytes=largest_bytes=0; }
void rdm_fail_strings_allocation(int n) { rdm_fail_allocation(n); strings_only=1; }
void rdm_fail_transform_allocation(int n) { rdm_fail_allocation(n); strings_only=2; }
size_t rdm_total_bytes(void) { return total_bytes; }
size_t rdm_largest_bytes(void) { return largest_bytes; }
int rdm_allocation_count(void) { active=0; return count; }
static int fail(void *caller, size_t n) {
    Dl_info info;
    if (active && dladdr(caller, &info) && info.dli_fname &&
        strstr(info.dli_fname, "liboqp") && (!strings_only ||
        (info.dli_sname && strings_only==1 && strstr(info.dli_sname, "_MOD_rdm12_strings")) ||
        (info.dli_sname && strings_only==2 && !strcmp(info.dli_sname, "mo_transform_eri")))) {
        total_bytes += n; if (n > largest_bytes) largest_bytes=n;
        if (++count == selected) { active=0; errno=ENOMEM; return 1; }
    }
    return 0;
}
#ifdef __APPLE__
static void *replacement_malloc(size_t n) {
    if (fail(__builtin_return_address(0), n)) return NULL;
    return malloc(n);
}
__attribute__((used)) static struct { const void *replacement, *original; }
interpose __attribute__((section("__DATA,__interpose"))) = {
    (const void *)replacement_malloc, (const void *)malloc
};
#else
extern void *__libc_malloc(size_t);
void *malloc(size_t n) {
    if (fail(__builtin_return_address(0), n)) return NULL;
    return __libc_malloc(n);
}
#endif
''')
    if sys.platform == "darwin":
        compiler = "/opt/homebrew/bin/gcc-15" if os.uname().machine == "arm64" else "/usr/local/bin/gcc-15"
        command = [compiler, "-dynamiclib", "-O0", str(source), "-o", str(library)]
        preload_key = "DYLD_INSERT_LIBRARIES"
    else:
        command = ["cc", "-shared", "-fPIC", "-O0", str(source), "-ldl", "-o", str(library)]
        preload_key = "LD_PRELOAD"
    subprocess.run(command, check=True, capture_output=True, text=True)
    env = os.environ.copy()
    env[preload_key] = str(library)
    return library, env


@pytest.mark.parametrize("threads", [1, 4])
def test_each_real_allocation_failure_restores_caller(native, allocation_interposer, threads):
    library, env = allocation_interposer
    child = r'''
import ctypes, sys
import numpy as np
from oqp.library import rdm
from oqp import runtime
lib, ffi = rdm._rdm_backend()
ctl = ctypes.CDLL(str(runtime.library_path(*runtime.resolve_oqp_root())))
ctl.oqp_blas_thread_count.restype = ctypes.c_int64
ctl.oqp_blas_thread_set.argtypes = [ctypes.c_int64]
inject = ctypes.CDLL(sys.argv[1])
inject.rdm_fail_allocation.argtypes = [ctypes.c_int]
saved = ctl.oqp_blas_thread_count()
if saved > 0:
    ctl.oqp_blas_thread_set(3)
initial = ctl.oqp_blas_thread_count()
dets = np.array(rdm.determinant_basis(4, (2, 1)), dtype=np.int64)
c = np.ones(len(dets))/np.sqrt(len(dets))
d1 = np.empty((4,4)); d2 = np.empty((4,)*4)
args = (4, len(dets), ffi.cast("int64_t *", dets.ctypes.data),
        ffi.cast("double *", c.ctypes.data), ffi.cast("double *", d1.ctypes.data),
        ffi.cast("double *", d2.ctypes.data), int(sys.argv[2]))
assert lib.rdm12_strings_c(*args) == 0
ref1, ref2 = d1.copy(), d2.copy()
try:
    for ordinal in range(1, 19):
        d1.fill(23); d2.fill(29)
        inject.rdm_fail_allocation(ordinal)
        info = lib.rdm12_strings_c(*args)
        count = inject.rdm_allocation_count()
        assert count == ordinal, (ordinal, count)
        assert info == -1, (ordinal, info)
        assert np.all(d1 == 23) and np.all(d2 == 29), ordinal
        assert ctl.oqp_blas_thread_count() == initial, ordinal
        # Reuse caller-owned output buffers; native locals allocate fresh on recovery.
        for repeat in range(2):
            assert lib.rdm12_strings_c(*args) == 0, ordinal
            assert ctl.oqp_blas_thread_count() == initial, ordinal
            assert d1.tobytes() == ref1.tobytes()
            assert d2.tobytes() == ref2.tobytes()
    print("18 allocation failures propagated; caller BLAS restored; two recoveries per allocation passed")
finally:
    ctl.oqp_blas_thread_set(saved)
'''
    result = subprocess.run([sys.executable, "-c", child, str(library), str(threads)], env=env,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "two recoveries per allocation passed" in result.stdout
    print(result.stdout)



@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("failure_target", ["rdm", "transform"])
def test_persistent_casscf_real_allocator_recovery(native, allocation_interposer, tmp_path, threads,
                                                 failure_target):
    """One Runner/Molecule: success, native failure, then two top-level recoveries."""
    library, env = allocation_interposer
    child = r'''
import ctypes, pathlib, sys
import numpy as np
from oqp.library import casscf as cs
from oqp.pyoqp import Runner
requested_threads = int(sys.argv[3])
failure_target = sys.argv[4]
shim = ctypes.CDLL(sys.argv[1])
shim.rdm_fail_strings_allocation.argtypes = [ctypes.c_int]
shim.rdm_fail_transform_allocation.argtypes = [ctypes.c_int]
real_backend = cs._casscf_energy_backend
lib, ffi = real_backend()
for name in ("casscf_energy",):
    assert hasattr(lib, name), "Missing native entry point: " + name
initial_mo = None
controls = None
seen = []
phase = "success"
ordinal = 0
optimized_mo = None
transform_allocations = None
class FailRDM:
    def casscf_energy(self, *args):
        global initial_mo, controls, optimized_mo, transform_allocations
        # All four calls share a Molecule. Pin the native input orbitals after
        # the top-level SCF, so separate SCF reductions cannot mask recovery.
        if initial_mo is None:
            initial_mo = runner.mol.data["OQP::VEC_MO_A"].copy()
        runner.mol.data["OQP::VEC_MO_A"][...] = initial_mo
        current = tuple(bytes(ffi.buffer(args[i], n))
                        for i, n in ((1, len(cs._CAS_IOPT)*4),
                                     (2, len(cs._CAS_DOPT)*8), (3, 8), (4, 4)))
        iopt = np.frombuffer(current[0], dtype=np.int32)
        assert iopt[cs._CAS_IOPT_INDEX["nthreads"]] == requested_threads
        if controls is None:
            controls = current
        assert current == controls, "Recovery changed native controls"
        if phase == "failure":
            if failure_target == "transform":
                shim.rdm_fail_transform_allocation(ordinal)
            else:
                shim.rdm_fail_strings_allocation(ordinal)
            try:
                status = int(lib.casscf_energy(*args))
            finally:
                count = shim.rdm_allocation_count()
            assert count == ordinal, (status, count, ordinal)
            assert status == -2, status
            # Reacquire the native record; never keep a borrowed view across calls.
            post_commit = failure_target == "transform" and ordinal == transform_allocations
            expected = optimized_mo if post_commit else initial_mo
            assert runner.mol.data["OQP::VEC_MO_A"].tobytes() == expected.tobytes()
        else:
            if failure_target == "transform":
                shim.rdm_fail_transform_allocation(0)
            try:
                status = int(lib.casscf_energy(*args))
            finally:
                if failure_target == "transform":
                    count = shim.rdm_allocation_count()
            assert status == 0, status
            if failure_target == "transform":
                # A converged control may build multiple analytic Hessians.
                # Each owns one transform, followed by canonicalization and
                # post-commit final CI; select those sites from native stats.
                assert args[8][2] == 1, "Control must converge at the unchanged tolerance"
                assert args[8][4] >= 1, "Control must enter the analytic Hessian"
                assert count == int(args[8][4]) + 2, (count, args[8][4])
                if transform_allocations is None:
                    transform_allocations = count
                assert count == transform_allocations, (count, transform_allocations)
            native_energies.append(bytes(ffi.buffer(args[5], 8)))
            if optimized_mo is None:
                optimized_mo = runner.mol.data["OQP::VEC_MO_A"].copy()
        seen.append(status)
        return status
config = {
    "input": {"system": "\nLi 0 0 0\nH 0 0 1.6", "basis": "sto-3g",
              "method": "casscf", "runtype": "energy", "omp_threads": str(requested_threads)},
    "guess": {"type": "hcore"},
    "scf": {"type": "rhf", "multiplicity": "1", "save_molden": "False"},
    "cas": {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1"},
    "ci": {"nroot": "1", "solver": "dense", "target_spin": "any",
           "integral_backend": "native", "root_tracking": "energy"},
    "casscf": {"max_macro_iterations": "50", "canonicalize": "true"},
    "tests": {"exception": "True"},
}
if failure_target == "transform":
    config["casscf"].update({"converger": "ah", "hessian": "analytic"})
runner = Runner(project="persistent-recovery", input_file=None, silent=1, usempi=False,
                input_dict=config, log=str(pathlib.Path(sys.argv[2])/"recovery.log"))
molecule = runner.mol
native_energies = []
real_optimize = cs._optimize
def forbidden_fallback(*args, **kwargs):
    raise AssertionError("Allocation failure entered the Python optimizer")
cs._optimize = forbidden_fallback
cs._casscf_energy_backend = lambda: (FailRDM(), ffi)
try:
    runner.run(test_mod=True)
    assert runner.mol is molecule
    if failure_target == "transform":
        # First Hessian, canonicalization, and final CI, independent
        # of the number of Hessian builds needed for convergence.
        failure_ordinals = (1, transform_allocations - 1, transform_allocations)
        assert len(set(failure_ordinals)) == 3
    else:
        failure_ordinals = (1,)
    phases = [(name, failure_ordinal) for failure_ordinal in failure_ordinals
              for name in ("failure", "recovery1", "recovery2")]
    for phase, ordinal in phases:
        if phase == "failure":
            try:
                runner.run(test_mod=True)
            except MemoryError as exc:
                assert "Native CASSCF scratch" in str(exc)
            else:
                raise AssertionError("Native allocation failure was swallowed")
        else:
            runner.run(test_mod=True)
        assert runner.mol is molecule
finally:
    cs._casscf_energy_backend = real_backend
    cs._optimize = real_optimize
assert seen == [0] + [-2, 0, 0]*(3 if failure_target == "transform" else 1), seen
assert all(value == native_energies[0] for value in native_energies), native_energies
print("native CASSCF status -2 propagated; persistent top-level recoveries compared",
      "target=", failure_target, "threads=", requested_threads,
      "ordinals=", failure_ordinals, "statuses=", seen)
'''
    env["OMP_NUM_THREADS"] = str(threads)
    env["OQP_FCI_THREADS"] = str(threads)
    result = subprocess.run([sys.executable, "-c", child, str(library), str(tmp_path), str(threads),
                             failure_target],
                            env=env, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "persistent top-level recoveries compared" in result.stdout
    print(result.stdout)
