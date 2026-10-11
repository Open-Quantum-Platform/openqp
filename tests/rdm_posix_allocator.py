"""GNU/glibc Linux and macOS compiler-allocator regression lane.

Explicitly consumed by .github/workflows/CI.yml on GNU Linux. This filename
keeps a platform-specific interposer out of Windows' portable test collection;
invoke it explicitly on supported managed macOS builds too. Unsupported hosts
fail if asked to run it. No synthetic controller is armed in these tests.
"""

import ctypes
import os
from pathlib import Path
import subprocess
import sys

import pytest

from test_rdm_fallback import native


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("norb,nelec", [(4, (2, 1)), (8, (4, 4))])
def test_rdm_never_changes_foreign_blas_policy(native, tmp_path, threads, norb, nelec):
    """Count setter calls while real RDM arithmetic uses the installed BLAS.

    A separate process resolves the control fixture before loading liboqp.
    The fixture changes only the control API; no numerical BLAS is replaced.
    Both sides of the T-construction OpenMP threshold must leave a foreign
    caller's temporary setting untouched, including overlapping RDM calls.
    """
    assert sys.platform in ("linux", "darwin"), "BLAS control lane requires POSIX"
    source = Path(__file__).with_name("ecp_thread_control_fixture.c")
    library = tmp_path / ("controls.dylib" if sys.platform == "darwin" else "controls.so")
    if sys.platform == "darwin":
        compiler = "/opt/homebrew/bin/gcc-15" if os.uname().machine == "arm64" else "/usr/local/bin/gcc-15"
        command = [compiler, "-dynamiclib", str(source), "-o", str(library)]
    else:
        command = ["cc", "-shared", "-fPIC", str(source), "-o", str(library)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    child = r'''
import ctypes, sys
fixture = ctypes.CDLL(sys.argv[1], mode=ctypes.RTLD_GLOBAL)
fixture.openblas_set_num_threads.argtypes = [ctypes.c_int]
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier
from oqp.library import rdm
from oqp import runtime
lib, ffi = rdm._rdm_backend()
assert hasattr(lib, "rdm12_strings_c"), "Missing string RDM export"
for name in ("oqp_rdm_test_alloc_begin", "oqp_rdm_test_alloc_end",
             "oqp_rdm_test_contraction_record"):
    assert hasattr(lib, name), "Missing contraction observation export: " + name
ctl = ctypes.CDLL(str(runtime.library_path(*runtime.resolve_oqp_root())))
ctl.oqp_blas_thread_count.restype = ctypes.c_int64
ctl.oqp_have_openmp.restype = ctypes.c_int
have_openmp = ctl.oqp_have_openmp()
threads, norb, na, nb = map(int, sys.argv[2:])
dets = np.ascontiguousarray(rdm.determinant_basis(norb, (na, nb)), dtype=np.int64)
vectors = [np.random.default_rng(742+i).standard_normal(len(dets)) for i in range(2)]
vectors = [np.ascontiguousarray(c / np.linalg.norm(c)) for c in vectors]
team_log = []
def evaluate(c):
    d1 = np.empty((norb, norb))
    d2 = np.empty((norb,) * 4)
    policy = ctl.oqp_blas_thread_count()
    assert lib.oqp_rdm_test_alloc_begin(0) == 0
    stats = ffi.new("struct oqp_rdm_test_alloc_stats *")
    try:
        status = lib.rdm12_strings_c(norb, len(dets), ffi.cast("int64_t *", dets.ctypes.data),
            ffi.cast("double *", c.ctypes.data), ffi.cast("double *", d1.ctypes.data),
            ffi.cast("double *", d2.ctypes.data), threads)
    finally:
        assert lib.oqp_rdm_test_alloc_end(stats) == 0
    assert status == 0, status
    expected_team = threads if have_openmp and policy == 1 and norb == 8 else 1
    assert stats.max_contraction_team == expected_team, (policy, stats.max_contraction_team, expected_team)
    team_log.append((policy, int(stats.max_contraction_team)))
    return d1, d2
assert ctl.oqp_blas_thread_count() == 4, "Control fixture was not resolved"
assert fixture.ecp_fixture_calls() == 0
reference = [evaluate(c) for c in vectors]
assert fixture.ecp_fixture_calls() == 0, "RDM mutated process-wide BLAS policy"
# Another native kernel has saved 4 and temporarily selected 2. RDM must
# participate in neither that save nor its later restore, in any build.
fixture.openblas_set_num_threads(2)
barrier = Barrier(2)
def overlapping(c):
    barrier.wait(timeout=10)
    return evaluate(c)
with ThreadPoolExecutor(max_workers=2) as pool:
    actual = list(pool.map(overlapping, vectors))
assert fixture.ecp_fixture_calls() == 1, "RDM interfered with the foreign transaction"
assert ctl.oqp_blas_thread_count() == 2
fixture.openblas_set_num_threads(4)
assert fixture.ecp_fixture_calls() == 2
assert ctl.oqp_blas_thread_count() == 4
results = [actual]
# Exercise the OpenMP contraction dispatch for a confirmed one-thread BLAS
# and the conservative unknown-policy branch, on the same input intervals.
for policy, calls in ((1, 3), (-1, 4)):
    fixture.openblas_set_num_threads(policy)
    assert ctl.oqp_blas_thread_count() == policy
    with ThreadPoolExecutor(max_workers=2) as pool:
        results.append(list(pool.map(overlapping, vectors)))
    assert fixture.ecp_fixture_calls() == calls, (policy, calls)
    assert ctl.oqp_blas_thread_count() == policy
fixture.openblas_set_num_threads(4)
assert fixture.ecp_fixture_calls() == 5
for c, previous, *variants in zip(vectors, reference, *results):
    for result in variants:
        assert result[0].tobytes() == previous[0].tobytes()
        assert result[1].tobytes() == previous[1].tobytes()
    result = variants[-1]
    spin = rdm._lib_rdm2(c, dets, 2*norb)
    assert spin is not None
    signed = np.zeros((norb,) * 4)
    for a in (0, norb):
        for b in (0, norb):
            signed += spin[a:a+norb,b:b+norb,a:a+norb,b:b+norb].transpose(0,2,1,3)
    assert np.max(np.abs(signed - signed.swapaxes(0,1))) > 1e-4
    np.testing.assert_allclose(result[1], signed, rtol=0, atol=2e-12)
print("RDM setter calls=0; foreign BLAS policy and overlapping signed pins passed",
      "threads=", threads, "norb=", norb, "control_widths=", [4, 2, 1, -1],
      "contraction_teams=", team_log)
'''
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["OMP_DYNAMIC"] = "FALSE"
    env["OQP_RDM_TEST_ALLOCATIONS"] = "1"
    result = subprocess.run([sys.executable, "-c", child, str(library), str(threads),
                             str(norb), str(nelec[0]), str(nelec[1])],
                            env=env, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    print(result.stdout)


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
#include <stdint.h>
static int selected, count, active, strings_only;
static int final_transform, transform_count;
static int evaluate_phase, evaluate_markers;
static int32_t (*original_fault_request)(int32_t);
static int32_t replacement_fault_request(int32_t site);
int32_t oqp_casscf_test_fault_request(int32_t site);
int rdm_set_fault_forward(void *original) {
#ifdef __APPLE__
    /* References from the interposer image retain the original dyld binding,
     * as does malloc below. The owning dylib is a linked dependency. */
    original = (void *)oqp_casscf_test_fault_request;
#endif
    if (!original || original == (void *)replacement_fault_request) return 1;
#ifndef __APPLE__
    if (original == (void *)oqp_casscf_test_fault_request) return 1;
#endif
    original_fault_request = (int32_t (*)(int32_t))original;
    return 0;
}
static int32_t replacement_fault_request(int32_t site) {
    if (!original_fault_request) abort();
    if (active && strings_only==4 && site==4 && !evaluate_phase) {
        evaluate_phase=1; evaluate_markers++;
    }
    return original_fault_request(site);
}
static size_t final_ci_bytes[3];
static size_t total_bytes, largest_bytes;
void rdm_fail_allocation(int n) { selected=n; count=0; active=1; strings_only=0; total_bytes=largest_bytes=0; }
void rdm_fail_strings_allocation(int n) { rdm_fail_allocation(n); strings_only=1; }
void rdm_fail_transform_allocation(int n) { rdm_fail_allocation(n); strings_only=2; }
void rdm_fail_final_ci_allocation(int transform, int n) {
    rdm_fail_allocation(n); strings_only=3;
    final_transform=transform; transform_count=0;
    memset(final_ci_bytes, 0, sizeof final_ci_bytes);
}
void rdm_fail_evaluate_ci_allocation(int n) {
    rdm_fail_final_ci_allocation(0, n); strings_only=4;
    evaluate_phase=evaluate_markers=0;
}
int rdm_evaluate_ci_markers(void) { return evaluate_markers; }
int rdm_final_ci_transform_count(void) { return transform_count; }
size_t rdm_final_ci_allocation_bytes(int i) { return i>=0 && i<3 ? final_ci_bytes[i] : 0; }
size_t rdm_total_bytes(void) { return total_bytes; }
size_t rdm_largest_bytes(void) { return largest_bytes; }
int rdm_allocation_count(void) { active=0; return count; }
static int fail(void *caller, size_t n) {
    Dl_info info;
    if (active && dladdr(caller, &info) && info.dli_fname && strstr(info.dli_fname, "liboqp")) {
        if (strings_only==3) {
            /* The control counts one allocation per transform. After its
             * selected transform, fci_solve allocates xor_code, orb_code, dets.
             * Reject a changed allocation layout instead of failing an
             * unchecked allocation: Python requires the actual failure. */
            if (transform_count < final_transform) {
                if (info.dli_sname && !strcmp(info.dli_sname, "mo_transform_eri")) transform_count++;
                return 0;
            }
            if (count < 3) final_ci_bytes[count]=n;
        } else if (strings_only==4) {
            /* The caller-owned observation opens the first active CI phase.
             * Private allocator symbols and inlining do not select it. */
            if (!evaluate_phase) return 0;
            if (count < 3) final_ci_bytes[count]=n;
        } else if (strings_only==1 && (!info.dli_sname || !strstr(info.dli_sname, "_MOD_rdm12_strings"))) {
            return 0;
        } else if (strings_only==2 && (!info.dli_sname || strcmp(info.dli_sname, "mo_transform_eri"))) {
            return 0;
        }
        total_bytes += n; if (n > largest_bytes) largest_bytes=n;
        if (++count == selected) {
            active=0;
            /* ILP64 CAS(2,2): two one-integer arrays, then four int64 dets.
             * A layout mismatch fails the Python guard without injecting. */
            if ((strings_only==3 || strings_only==4) && (selected!=3 || final_ci_bytes[0]!=8 ||
                final_ci_bytes[1]!=8 || final_ci_bytes[2]!=32)) return 0;
            errno=ENOMEM; return 1;
        }
    }
    return 0;
}
#ifdef __APPLE__
static void *replacement_malloc(size_t n) {
    if (fail(__builtin_return_address(0), n)) return NULL;
    return malloc(n);
}
__attribute__((used)) static struct { const void *replacement, *original; }
interpose[] __attribute__((section("__DATA,__interpose"))) = {
    {(const void *)replacement_malloc, (const void *)malloc},
    {(const void *)replacement_fault_request, (const void *)oqp_casscf_test_fault_request}
};
#else
extern void *__libc_malloc(size_t);
int32_t oqp_casscf_test_fault_request(int32_t site) {
    return replacement_fault_request(site);
}
void *malloc(size_t n) {
    if (fail(__builtin_return_address(0), n)) return NULL;
    return __libc_malloc(n);
}
#endif
''')
    if sys.platform == "darwin":
        from oqp import runtime
        owner = Path(runtime.library_path(*runtime.resolve_oqp_root()))
        compiler = "/opt/homebrew/bin/gcc-15" if os.uname().machine == "arm64" else "/usr/local/bin/gcc-15"
        # Load the hook's actual owner before dyld resolves the interpose tuple.
        # An unresolved dynamic-lookup symbol cannot satisfy this address
        # relocation when liboqp would otherwise be loaded later by Python.
        command = [compiler, "-dynamiclib", "-O0", str(source), str(owner),
                   "-Wl,-rpath," + str(owner.parent), "-o", str(library)]
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
@pytest.mark.parametrize("failure_target", ["rdm", "transform", "ci-evaluate", "ci-canonical", "final-ci"])
def test_persistent_casscf_real_allocator_recovery(native, allocation_interposer, tmp_path, threads,
                                                 failure_target):
    """One Runner/Molecule: success, native failure, then two top-level recoveries."""
    library, env = allocation_interposer
    child = r'''
import ctypes, pathlib, sys
import numpy as np
from oqp.library import casscf as cs
from oqp import runtime
from oqp.pyoqp import Runner
requested_threads = int(sys.argv[3])
failure_target = sys.argv[4]
shim = ctypes.CDLL(sys.argv[1])
shim.rdm_fail_strings_allocation.argtypes = [ctypes.c_int]
shim.rdm_fail_transform_allocation.argtypes = [ctypes.c_int]
shim.rdm_fail_final_ci_allocation.argtypes = [ctypes.c_int, ctypes.c_int]
shim.rdm_fail_evaluate_ci_allocation.argtypes = [ctypes.c_int]
shim.rdm_final_ci_allocation_bytes.argtypes = [ctypes.c_int]
shim.rdm_final_ci_allocation_bytes.restype = ctypes.c_size_t
real_backend = cs._casscf_energy_backend
lib, ffi = real_backend()
ctl = ctypes.CDLL(str(runtime.library_path(*runtime.resolve_oqp_root())))
shim.rdm_set_fault_forward.argtypes = [ctypes.c_void_p]
# Resolve the original hook from its owning library before malloc is armed.
assert shim.rdm_set_fault_forward(ctypes.cast(ctl.oqp_casscf_test_fault_request,
                                            ctypes.c_void_p)) == 0
for name in ("casscf_energy",):
    assert hasattr(lib, name), "Missing native entry point: " + name
initial_mo = None
controls = None
seen = []
phase = "success"
ordinal = 0
optimized_mo = None
transform_allocations = None
final_ci_layouts = []
ci_targets = {"ci-evaluate", "ci-canonical", "final-ci"}
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
            args[5][0], args[6][0] = 13.25, -3.5
            result_before = tuple(bytes(ffi.buffer(args[i], 8)) for i in (5, 6))
            if failure_target == "transform":
                shim.rdm_fail_transform_allocation(ordinal)
            elif failure_target == "ci-evaluate":
                shim.rdm_fail_evaluate_ci_allocation(ordinal)
            elif failure_target in ci_targets:
                selected_transform = transform_allocations - (failure_target == "ci-canonical")
                shim.rdm_fail_final_ci_allocation(selected_transform, ordinal)
            else:
                shim.rdm_fail_strings_allocation(ordinal)
            try:
                status = int(lib.casscf_energy(*args))
            finally:
                count = shim.rdm_allocation_count()
            assert count == ordinal, (status, count, ordinal)
            if failure_target in ci_targets:
                if failure_target == "ci-evaluate":
                    assert shim.rdm_evaluate_ci_markers() == 1, "Missing first-CI observation"
                expected_transform = (0 if failure_target == "ci-evaluate" else
                                      transform_allocations - (failure_target == "ci-canonical"))
                assert shim.rdm_final_ci_transform_count() == expected_transform
                sizes = tuple(shim.rdm_final_ci_allocation_bytes(i) for i in range(3))
                assert sizes == (8, 8, 32), ("CI allocation layout changed", failure_target, sizes, status)
                final_ci_layouts.append((failure_target, expected_transform, sizes))
            assert status == -2, status
            # Reacquire the native record; never keep a borrowed view across calls.
            assert runner.mol.data["OQP::VEC_MO_A"].tobytes() == initial_mo.tobytes(), (
                "Failed native CASSCF committed orbitals", failure_target, ordinal)
            assert tuple(bytes(ffi.buffer(args[i], 8)) for i in (5, 6)) == result_before, (
                "Failed native CASSCF published energies or spin", failure_target, ordinal)
        else:
            if failure_target != "rdm":
                shim.rdm_fail_transform_allocation(0)
            try:
                status = int(lib.casscf_energy(*args))
            finally:
                if failure_target != "rdm":
                    count = shim.rdm_allocation_count()
            assert status == 0, status
            if failure_target != "rdm":
                # A converged control may build multiple analytic Hessians.
                # Each owns one transform, followed by canonicalization and
                # final CI; select those sites from native stats.
                assert args[8][2] == 1, "Control must converge at the unchanged tolerance"
                assert args[8][4] >= 1, "Control must enter the analytic Hessian"
                assert count == int(args[8][4]) + 2, (count, args[8][4])
                if transform_allocations is None:
                    transform_allocations = count
                assert count == transform_allocations, (count, transform_allocations)
            native_results.append(tuple(bytes(ffi.buffer(args[i], 8)) for i in (5, 6)))
            if optimized_mo is None:
                optimized_mo = runner.mol.data["OQP::VEC_MO_A"].copy()
                if failure_target != "rdm":
                    assert optimized_mo.tobytes() != initial_mo.tobytes(), (
                        "Control must change orbitals to expose a premature commit")
            assert runner.mol.data["OQP::VEC_MO_A"].tobytes() == optimized_mo.tobytes(), (
                "Recovery changed optimized orbitals", failure_target, ordinal)
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
if failure_target != "rdm":
    config["casscf"].update({"converger": "ah", "hessian": "analytic"})
runner = Runner(project="persistent-recovery", input_file=None, silent=1, usempi=False,
                input_dict=config, log=str(pathlib.Path(sys.argv[2])/"recovery.log"))
molecule = runner.mol
native_results = []
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
    elif failure_target in ci_targets:
        # Real malloc failure on dets(ndet), the first checked CI allocation.
        failure_ordinals = (3,)
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
assert all(value == native_results[0] for value in native_results), native_results
print("native CASSCF status -2 propagated; persistent top-level recoveries compared",
      "target=", failure_target, "threads=", requested_threads,
      "ordinals=", failure_ordinals, "statuses=", seen,
      "final_ci_layouts=", final_ci_layouts)
'''
    env["OMP_NUM_THREADS"] = str(threads)
    env["OQP_FCI_THREADS"] = str(threads)
    result = subprocess.run([sys.executable, "-c", child, str(library), str(tmp_path), str(threads),
                             failure_target],
                            env=env, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "persistent top-level recoveries compared" in result.stdout
    print(result.stdout)
