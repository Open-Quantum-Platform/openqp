"""Native string RDM acceptance, numerical pin, and caller-state lifetime."""

import ctypes
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager

import numpy as np
import pytest

from oqp.library import rdm
from oqp import runtime


@pytest.fixture
def native():
    backend = rdm._rdm_backend()
    assert backend is not None, "RDM tests require the installed native library"
    lib, ffi = backend
    assert hasattr(lib, "rdm12_strings_c"), "Missing string RDM entry point"
    assert hasattr(lib, "rdm2_spatial"), "Missing Gram RDM entry point"
    path = runtime.library_path(*runtime.resolve_oqp_root())
    ctl = ctypes.CDLL(str(path))
    ctl.oqp_blas_thread_count.restype = ctypes.c_int64
    ctl.oqp_blas_thread_set.argtypes = [ctypes.c_int64]
    return lib, ffi, ctl


@contextmanager
def _caller_blas(native):
    ctl = native[2]
    saved = ctl.oqp_blas_thread_count()
    try:
        if saved > 0:
            ctl.oqp_blas_thread_set(3)
        requested = ctl.oqp_blas_thread_count()
        if os.environ.get("OQP_RDM_REQUIRE_BLAS3") == "1":
            assert requested == 3, "This validation lane requires BLAS threads >1"
        yield
        assert ctl.oqp_blas_thread_count() == requested
    finally:
        ctl.oqp_blas_thread_set(saved)


def _ci(norb, nelec, seed=742):
    dets = np.ascontiguousarray(rdm.determinant_basis(norb, nelec), dtype=np.int64)
    c = np.random.default_rng(seed).standard_normal(len(dets))
    c /= np.linalg.norm(c)
    return dets, np.ascontiguousarray(c)


def _strings(native, norb, dets, c, threads, d1=None, d2=None):
    lib, ffi, _ = native
    if d1 is None:
        d1 = np.full((norb, norb), 23.0)
    if d2 is None:
        d2 = np.full((norb,) * 4, 29.0)
    info = lib.rdm12_strings_c(
        norb, len(dets), ffi.cast("int64_t *", dets.ctypes.data),
        ffi.cast("double *", c.ctypes.data),
        ffi.cast("double *", d1.ctypes.data), ffi.cast("double *", d2.ctypes.data), threads)
    return int(info), d1, d2


def _spin_gram(c, dets, norb):
    spin = rdm._lib_rdm2(c, dets, 2 * norb)
    assert spin is not None, "Independent spin-orbital Gram pin must execute"
    ref = np.zeros((norb,) * 4)
    for a in (0, norb):
        for b in (0, norb):
            ref += spin[a:a+norb, b:b+norb, a:a+norb, b:b+norb].transpose(0, 2, 1, 3)
    return ref


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("norb,nelec", [(4, (2, 1)), (8, (4, 4))])
def test_canonical_strings_match_signed_gram(native, threads, norb, nelec):
    # 24*16 < 65536, 4900*64 > 65536: both sides of the OpenMP condition.
    dets, c = _ci(norb, nelec)
    with _caller_blas(native):
        info, d1, d2 = _strings(native, norb, dets, c, threads)
    assert info == 0
    ref1 = rdm._lib_rdm1_spatial(c, dets, norb)
    ref2 = _spin_gram(c, dets, norb)
    assert np.max(np.abs(ref2 - ref2.swapaxes(0, 1))) > 1e-4
    np.testing.assert_allclose(d1, ref1, rtol=0, atol=2e-12)
    np.testing.assert_allclose(d2, ref2, rtol=0, atol=2e-12)
    np.testing.assert_allclose(np.einsum("pqrr->pq", d2), (sum(nelec)-1)*d1,
                               rtol=0, atol=2e-12)
    # The default Gram route retains its stronger, pre-existing bit contract.
    assert np.array_equal(rdm.make_rdm2_spatial(c, dets, norb), ref2)


@pytest.mark.parametrize("kind", ["truncated", "permuted"])
@pytest.mark.parametrize("threads", [1, 4])
def test_noncanonical_declines_without_mutating_caller(native, kind, threads):
    dets, c = _ci(4, (2, 2))
    if kind == "truncated":
        dets, c = dets[:7].copy(), c[:7].copy()
        c /= np.linalg.norm(c)
    else:
        order = np.random.default_rng(63).permutation(len(dets))
        dets, c = dets[order].copy(), c[order].copy()
    with _caller_blas(native):
        info, d1, d2 = _strings(native, 4, dets, c, threads)
        assert rdm.make_rdm12_spatial_strings(c, dets, 4) is None
    assert info > 0
    assert np.all(d1 == 23) and np.all(d2 == 29)
    assert np.array_equal(rdm.make_rdm2_spatial(c, dets, 4), _spin_gram(c, dets, 4))


@pytest.mark.parametrize("norb,dets", [(0, [0]), (-1, [0]), (32, [0]),
                                      (64, [0]), (2, [1 << 10]), (2, [])])
def test_native_declining_boundaries_preserve_buffers_and_threads(native, norb, dets):
    dets = np.ascontiguousarray(dets, dtype=np.int64)
    c = np.ones(len(dets))
    d1, d2 = np.array([23.0]), np.array([29.0])
    with _caller_blas(native):
        info, _, _ = _strings(native, norb, dets, c, 4, d1, d2)
    assert info > 0
    assert d1[0] == 23 and d2[0] == 29


@pytest.mark.parametrize("nelec", [(0, 0), (1, 0), (31, 31)])
def test_supported_upper_boundary(native, nelec):
    dets, c = _ci(31, nelec)
    with _caller_blas(native):
        info, d1, d2 = _strings(native, 31, dets, c, 4)
    assert info == 0
    assert np.trace(d1) == pytest.approx(sum(nelec), abs=2e-12)
    if sum(nelec) <= 1:
        np.testing.assert_allclose(d2, 0, rtol=0, atol=2e-12)
    else:
        eye = np.eye(31)
        ref = 4*np.einsum("pq,rs->pqrs", eye, eye) - 2*np.einsum("ps,rq->pqrs", eye, eye)
        assert np.array_equal(d1, 2*eye)
        assert np.array_equal(d2, ref)


def test_python_propagates_native_allocation_failure(native, monkeypatch):
    class AllocationFailure:
        def rdm12_strings_c(self, *args):
            return -1

    monkeypatch.setattr(rdm, "_rdm_backend", lambda: (AllocationFailure(), native[1]))
    dets, c = _ci(4, (2, 1))
    with pytest.raises(MemoryError, match="RDM scratch"):
        rdm.make_rdm12_spatial_strings(c, dets, 4)


@pytest.mark.parametrize("threads", [1, 2])
def test_overlapping_rdm_calls_restore_blas(native, threads):
    inputs = [_ci(8, (4, 4), seed=742 + i) for i in range(4)]
    with _caller_blas(native):
        reference = [_strings(native, 8, dets, c, threads) for dets, c in inputs]
        with ThreadPoolExecutor(max_workers=2) as pool:
            runs = list(pool.map(lambda pair: _strings(native, 8, *pair, threads), inputs))
    assert all(result[0] == 0 for result in runs)
    assert all(result[0] == 0 for result in reference)
    for result, ref in zip(runs, reference):
        assert result[1].tobytes() == ref[1].tobytes()
        assert result[2].tobytes() == ref[2].tobytes()


def test_repeated_dimensions_and_short_final_block(native):
    lib = native[0]
    saved = lib.fci_get_work_bytes_cap()
    try:
        for norb, nelec in [(4, (2, 1)), (10, (5, 5)), (2, (1, 1)), (8, (4, 4)), (4, (2, 1))]:
            dets, c = _ci(norb, nelec)
            lib.fci_set_work_bytes_cap(0)
            with _caller_blas(native):
                info, ref1, ref2 = _strings(native, norb, dets, c, 1)
            assert info == 0
            # The capped CAS(10,10) call traverses several alpha-string
            # blocks, including a short final block; smaller cases fit whole.
            lib.fci_set_work_bytes_cap(16*1024**2)
            with _caller_blas(native):
                info, d1, d2 = _strings(native, norb, dets, c, 4)
            assert info == 0
            np.testing.assert_allclose(d1, ref1, rtol=0, atol=2e-12)
            np.testing.assert_allclose(d2, ref2, rtol=0, atol=2e-12)
    finally:
        lib.fci_set_work_bytes_cap(saved)
    assert lib.fci_get_work_bytes_cap() == saved


@pytest.fixture
def allocation_interposer(native, tmp_path):
    """Fail each of the 18 owned allocations, including both link helpers.

    The interposer is confined to a child and fails malloc calls originating
    from liboqp only. It never changes the installed library or compiler flags.
    """
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
size_t rdm_total_bytes(void) { return total_bytes; }
size_t rdm_largest_bytes(void) { return largest_bytes; }
int rdm_allocation_count(void) { active=0; return count; }
static int fail(void *caller, size_t n) {
    Dl_info info;
    if (active && dladdr(caller, &info) && info.dli_fname &&
        strstr(info.dli_fname, "liboqp") && (!strings_only ||
        (info.dli_sname && strstr(info.dli_sname, "_MOD_rdm12_strings")))) {
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


def test_each_native_allocation_failure_restores_caller(native, allocation_interposer):
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
        ffi.cast("double *", d2.ctypes.data), 4)
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
        # Reacquire all buffers on the next successful call in the same process.
        assert lib.rdm12_strings_c(*args) == 0, ordinal
    print("18 allocation failures propagated; caller BLAS restored; 18 recoveries passed")
finally:
    ctl.oqp_blas_thread_set(saved)
'''
    result = subprocess.run([sys.executable, "-c", child, str(library)], env=env,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "18 recoveries passed" in result.stdout


def test_requested_width_respects_scratch_budget(native, allocation_interposer):
    """Observe actual native allocation sizes, including a serial width of 64."""
    library, env = allocation_interposer
    child = r'''
import ctypes, os, sys
import numpy as np
from oqp.library import rdm
lib, ffi = rdm._rdm_backend()
shim = ctypes.CDLL(sys.argv[1])
shim.rdm_total_bytes.restype = ctypes.c_size_t
shim.rdm_largest_bytes.restype = ctypes.c_size_t
dets = np.array(rdm.determinant_basis(31, (1, 1)), dtype=np.int64)
c = np.ones(len(dets))/np.sqrt(len(dets))
d1 = np.empty((31,31)); d2 = np.empty((31,)*4)
args = (31, len(dets), ffi.cast("int64_t *", dets.ctypes.data),
        ffi.cast("double *", c.ctypes.data), ffi.cast("double *", d1.ctypes.data),
        ffi.cast("double *", d2.ctypes.data))
serial = os.environ.get("OQP_RDM_SERIAL_BUILD") == "1"
cap = 512*1024**2 if serial else 16*1024**2
lib.fci_set_work_bytes_cap(cap)
observed = []
for width in (1,64):
    shim.rdm_fail_allocation(0)
    assert lib.rdm12_strings_c(*args, width) == 0
    shim.rdm_allocation_count()
    observed.append((shim.rdm_total_bytes(), shim.rdm_largest_bytes()))
    assert observed[-1][1] <= 31**4*8, observed
    if width == 1:
        ref1, ref2 = d1.copy(), d2.copy()
    else:
        np.testing.assert_allclose(d1, ref1, rtol=0, atol=2e-12)
        np.testing.assert_allclose(d2, ref2, rtol=0, atol=2e-12)
if serial:
    assert observed[0] == observed[1], observed
print("native allocation sizes", observed, "serial", serial)
'''
    result = subprocess.run([sys.executable, "-c", child, str(library)], env=env,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr


def test_minimum_block_budget_failure_preserves_caller(native):
    lib = native[0]
    saved = lib.fci_get_work_bytes_cap()
    try:
        # One beta interval, Gram, and one partial exceed 16 MiB here.
        dets, c = _ci(31, (1, 2))
        lib.fci_set_work_bytes_cap(16*1024**2)
        with _caller_blas(native):
            info, d1, d2 = _strings(native, 31, dets, c, 4)
        assert info == -1
        assert np.all(d1 == 23) and np.all(d2 == 29)
        dets, c = _ci(4, (2, 1))
        assert _strings(native, 4, dets, c, 4)[0] == 0
    finally:
        lib.fci_set_work_bytes_cap(saved)


def test_native_casscf_propagates_rdm_allocation_failure(native, allocation_interposer, tmp_path):
    """Fail at the actual string RDM allocation inside the native optimizer."""
    library, env = allocation_interposer
    child = r'''
import ctypes, pathlib, sys
import numpy as np
from oqp.library import casscf as cs
from oqp.pyoqp import Runner
shim = ctypes.CDLL(sys.argv[1])
shim.rdm_fail_strings_allocation.argtypes = [ctypes.c_int]
real_backend = cs._casscf_energy_backend
lib, ffi = real_backend()
assert hasattr(lib, "casscf_energy"), "Missing native CASSCF entry point"
seen = []
class FailRDM:
    def casscf_energy(self, *args):
        # Pin recovery to the same SCF orbitals. Independent threaded SCF
        # runs need not return bit-identical orbitals or final energies.
        initial_mo = current_runner.mol.data["OQP::VEC_MO_A"].copy()
        assert lib.casscf_energy(*args) == 0
        reference_energy = bytes(ffi.buffer(args[5], 8))
        current_runner.mol.data["OQP::VEC_MO_A"][...] = initial_mo
        shim.rdm_fail_strings_allocation(1)
        status = int(lib.casscf_energy(*args))
        count = shim.rdm_allocation_count()
        assert count == 1, (status, count)
        assert status == -2, (status, count)
        assert current_runner.mol.data["OQP::VEC_MO_A"].tobytes() == initial_mo.tobytes()
        recovered = []
        for repeat in range(2):
            # Reacquire the borrowed record after each native call.
            current_runner.mol.data["OQP::VEC_MO_A"][...] = initial_mo
            assert lib.casscf_energy(*args) == 0
            recovered.append(bytes(ffi.buffer(args[5], 8)))
        assert recovered[0] == recovered[1] == reference_energy, recovered
        current_runner.mol.data["OQP::VEC_MO_A"][...] = initial_mo
        seen.append(status)
        return status
config = {
    "input": {"system": "\nLi 0 0 0\nH 0 0 1.6", "basis": "sto-3g",
              "method": "casscf", "runtype": "energy"},
    "guess": {"type": "hcore"},
    "scf": {"type": "rhf", "multiplicity": "1", "save_molden": "False"},
    "cas": {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1"},
    "ci": {"nroot": "1", "solver": "dense", "target_spin": "any",
           "integral_backend": "native", "root_tracking": "energy"},
    "casscf": {"max_macro_iterations": "50", "canonicalize": "true"},
    "tests": {"exception": "True"},
}
def run(name):
    global current_runner
    runner = Runner(project=name, input_file=None, silent=1, usempi=False,
                    input_dict=config, log=str(pathlib.Path(sys.argv[2])/f"{name}.log"))
    current_runner = runner
    runner.run(test_mod=True)
    return float(runner.mol.energies[0])
for attempt in range(2):
    cs._casscf_energy_backend = lambda: (FailRDM(), ffi)
    try:
        run(f"failed-{attempt}")
    except MemoryError as exc:
        assert "Native CASSCF scratch" in str(exc)
    else:
        raise AssertionError("Native allocation failure was swallowed")
    finally:
        cs._casscf_energy_backend = real_backend
assert seen == [-2, -2], seen
print("native CASSCF allocation status -2 propagated; two exact recovery comparisons")
'''
    result = subprocess.run([sys.executable, "-c", child, str(library), str(tmp_path)],
                            env=env, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "two exact recovery comparisons" in result.stdout
