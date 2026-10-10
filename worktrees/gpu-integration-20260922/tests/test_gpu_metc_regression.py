"""CPU-vs-GPU numerical regression test for the experimental METC kernel.

This locks the *kernel and indexing correctness* of ``oqp_gpu_metc_contract``
(``source/gpu_metc_cuda.cu``) against a pure-NumPy reimplementation of the exact
CPU contraction performed by ``int2_mrsf_data_t_update`` /
``int2_umrsf_data_t_update`` (``source/tdhf_mrsf_lib.F90``).

Why a unit test and not only an end-to-end run:
  * It is deterministic and pinpoints indexing/stride bugs directly. The kernel's
    ``idx4(f, m, row, col) = f + nf*(m + nmatrix*(row + nbf*col))`` is exactly the
    column-major layout of the Fortran array ``f3(nfocks, nmatrix, nbf, nbf)``
    (allocated in tdhf_mrsf_lib.F90), so a flat buffer round-trips losslessly and
    the NumPy reference below can mirror the Fortran statements one-for-one.
  * It exercises the 1-based -> 0-based index shift, every ``m``-range mask
    (``m<4``/``m<7``/``m<8``/``m in {8,9}``/``m==10``/``m==6``), both kernels
    (MRSF/UMRSF), both Davidson passes, and accumulation into a non-zero f3.

Tolerance, not bit-exactness: the GPU uses ``atomicAdd``, which sums each f3 cell
in a nondeterministic order. Floating-point addition is non-associative, so the
result differs from the CPU at the ULP level (the CPU path is itself not
bit-reproducible across OpenMP thread counts -- parallel_stop sums over threads).
We therefore compare with a tight relative tolerance.

The test auto-skips unless a standalone METC library is importable *and* a
working GPU is present (``oqp_gpu_metc_contract`` returns non-zero -> no device).
"""

import ctypes
import os
import unittest
from ctypes import POINTER, c_bool, c_double, c_int
from pathlib import Path

try:
    import numpy as np
except ImportError:  # pragma: no cover - numpy always available in the test env
    np = None

ROOT = Path(__file__).resolve().parents[1]

# Matrix-column count of the f3/d3 tensors (ubound(d3, 2)): 7 for MRSF, 11 for UMRSF.
NMATRIX = {"mrsf": 7, "umrsf": 11}


def _find_library():
    """Return a loaded OpenQP shared library exporting the METC symbol, or None."""
    candidates = []
    if os.environ.get("OPENQP_GPU_METC_LIB"):
        candidates.append(Path(os.environ["OPENQP_GPU_METC_LIB"]))
    for name in ("libopenqp_gpu_metc.so", "libopenqp_gpu_metc.dylib"):
        candidates.extend([ROOT / "build" / name, ROOT / "lib" / name])
    for path in candidates:
        if not path or not path.is_file():
            continue
        try:
            lib = ctypes.CDLL(str(path))
        except OSError:
            continue
        if hasattr(lib, "oqp_gpu_metc_contract"):
            return lib
    return None


def _bind(lib):
    fn = lib.oqp_gpu_metc_contract
    fn.restype = c_int
    fn.argtypes = [
        POINTER(c_int),     # ids   (4*ncur, 1-based, int32)
        POINTER(c_double),  # ints  (ncur)
        c_int,              # ncur
        POINTER(c_double),  # f3    (nf*nmatrix*nbf*nbf, in/out)
        POINTER(c_double),  # d3    (nf*nmatrix*nbf*nbf)
        c_int, c_int, c_int,  # nf, nmatrix, nbf
        c_int,              # cur_pass
        c_double, c_double,  # scale_exchange, scale_coulomb
        c_bool,             # is_umrsf
    ]
    return fn


# --- The contraction stencils, transcribed verbatim from tdhf_mrsf_lib.F90 -----
# Each entry is (dst_row, dst_col, src_row, src_col) given local indices i,j,k,l.

def _coulomb8(i, j, k, l):
    return [(i, j, k, l), (k, l, i, j), (i, j, l, k), (l, k, i, j),
            (j, i, k, l), (k, l, j, i), (j, i, l, k), (l, k, j, i)]


def _exchange8(i, j, k, l):
    return [(i, k, j, l), (k, i, l, j), (i, l, j, k), (l, i, k, j),
            (j, k, i, l), (k, j, l, i), (j, l, i, k), (l, j, k, i)]


def _block910(i, j, k, l):  # the UMRSF f3(:,9:10,...) stencil
    return [(i, l, k, j), (l, i, j, k), (k, j, i, l), (j, k, l, i),
            (i, k, l, j), (k, i, j, l), (l, j, i, k), (j, l, k, i)]


def reference_contract(mode, cur_pass, ids, ints, f3_init, d3, nf, nmatrix, nbf,
                       scale_exchange, scale_coulomb):
    """Pure-NumPy mirror of the CPU f3 update. Returns the updated flat f3."""
    f3 = f3_init.astype(np.float64).copy()
    ncur = len(ints)

    def idx(f, m, r, c):
        return f + nf * (m + nmatrix * (r + nbf * c))

    def apply(m, stencil, coeff):
        for (a, b, cc, dd) in stencil:
            for f in range(nf):
                f3[idx(f, m, a, b)] += coeff * d3[idx(f, m, cc, dd)]

    for n in range(ncur):
        i, j, k, l = (int(ids[4 * n + t]) - 1 for t in range(4))
        val = float(ints[n])
        xval = val * scale_exchange
        cval = val * scale_coulomb

        if mode == "mrsf":
            if cur_pass == 1:
                for m in range(0, 4):
                    apply(m, _coulomb8(i, j, k, l), cval)
                for m in range(0, 7):
                    apply(m, _exchange8(i, j, k, l), -xval)
            elif cur_pass == 2:
                apply(6, _exchange8(i, j, k, l), -xval)
        else:  # umrsf
            if cur_pass == 1:
                for m in range(0, 8):
                    apply(m, _coulomb8(i, j, k, l), cval)
                    apply(m, _exchange8(i, j, k, l), -xval)
                for m in (8, 9):
                    apply(m, _block910(i, j, k, l), -xval)
                apply(10, _exchange8(i, j, k, l), -xval)
            elif cur_pass == 2:
                apply(10, _exchange8(i, j, k, l), -xval)
    return f3


def _deterministic_inputs(mode, nf, nbf):
    """Build a fixed, reproducible (ids, ints, f3_init, d3) set.

    The integral list deliberately includes aliased quartets (i==k, j==l,
    diagonal i==j==k==l) to stress the in-thread overlap handling, plus index 1
    and index nbf to catch off-by-one stride errors.
    """
    nmatrix = NMATRIX[mode]
    # (i, j, k, l) 1-based; mix of distinct, partially-aliased and diagonal.
    quartets = [
        (1, 2, 3, 4),
        (2, 2, 1, 1),
        (1, 1, 1, 1),
        (nbf, nbf - 1, 2, 1),
        (3, 1, 3, 1),   # i==k, j==l
        (4, 4, 4, 1),
    ]
    quartets = [q for q in quartets if all(1 <= x <= nbf for x in q)]
    ncur = len(quartets)

    ids = np.empty(4 * ncur, dtype=np.int32)
    for n, q in enumerate(quartets):
        ids[4 * n:4 * n + 4] = q
    # Deterministic, well-separated magnitudes; no RNG (keeps the test stable).
    ints = np.array([0.37 - 0.11 * n + 0.013 * n * n for n in range(ncur)],
                    dtype=np.float64)

    size = nf * nmatrix * nbf * nbf
    t = np.arange(size, dtype=np.float64)
    d3 = np.cos(0.017 * t) * 0.5 + 0.001 * (t % 7)        # arbitrary but fixed
    f3_init = np.sin(0.013 * t) * 0.25                    # non-zero -> tests accumulation
    return ids, ints, f3_init.copy(), d3


@unittest.skipIf(np is None, "numpy is required")
class TestGpuMetcContractKernel(unittest.TestCase):
    """oqp_gpu_metc_contract must reproduce the CPU contraction within tolerance."""

    @classmethod
    def setUpClass(cls):
        cls.lib = _find_library()
        if cls.lib is None:
            if os.environ.get("OPENQP_GPU_METC_REQUIRE") == "1":
                raise AssertionError("required METC CUDA library unavailable")
            raise unittest.SkipTest(
                "no standalone METC library found "
                "(build with -DENABLE_CUDA=ON; set OPENQP_GPU_METC_LIB)")
        cls.lib.oqp_gpu_metc_device_count.restype = c_int
        if cls.lib.oqp_gpu_metc_device_count() < 1:
            if os.environ.get("OPENQP_GPU_METC_REQUIRE") == "1":
                raise AssertionError("required CUDA device unavailable")
            raise unittest.SkipTest("no CUDA device")
        cls.fn = _bind(cls.lib)

    def _run_case(self, mode, cur_pass, nf=2, nbf=4):
        nmatrix = NMATRIX[mode]
        ids, ints, f3_init, d3 = _deterministic_inputs(mode, nf, nbf)

        ref = reference_contract(mode, cur_pass, ids, ints, f3_init, d3,
                                 nf, nmatrix, nbf, 0.71, 1.13)

        gpu_f3 = f3_init.astype(np.float64).copy()
        ierr = self.fn(
            ids.ctypes.data_as(POINTER(c_int)),
            ints.ctypes.data_as(POINTER(c_double)),
            c_int(len(ints)),
            gpu_f3.ctypes.data_as(POINTER(c_double)),
            d3.ctypes.data_as(POINTER(c_double)),
            c_int(nf), c_int(nmatrix), c_int(nbf),
            c_int(cur_pass), c_double(0.71), c_double(1.13),
            c_bool(mode == "umrsf"),
        )
        self.assertEqual(ierr, 0, f"CUDA contraction failed for {mode} pass {cur_pass}")

        # Tolerance, not bit-exact: see module docstring (atomicAdd reordering).
        np.testing.assert_allclose(
            gpu_f3, ref, rtol=1e-9, atol=1e-11,
            err_msg=f"GPU {mode} pass {cur_pass} diverged from CPU reference")

    def test_mrsf_pass1(self):
        self._run_case("mrsf", 1)

    def test_mrsf_pass2(self):
        self._run_case("mrsf", 2)

    def test_umrsf_pass1(self):
        self._run_case("umrsf", 1)

    def test_umrsf_pass2(self):
        self._run_case("umrsf", 2)

    def test_mrsf_pass1_multi_fock(self):
        # nf > 1 stride check on the leading dimension.
        self._run_case("mrsf", 1, nf=3, nbf=5)



# Exercise every preserved kernel variant and the resident session ABI. These
# tests are kernel equivalence checks, not a validation of an electronic method.
class TestGpuMetcVariants(TestGpuMetcContractKernel):
    def test_all_variants(self):
        from unittest.mock import patch
        for variant in ('reference', 'combined_coulomb', 'combined_coulomb_warp', 'two_phase_accum'):
            with patch.dict(os.environ, {'OQP_GPU_METC_VARIANT': variant}):
                for mode in ('mrsf', 'umrsf'):
                    for cur_pass in (1, 2):
                        with self.subTest(variant=variant, mode=mode, cur_pass=cur_pass):
                            self._run_case(mode, cur_pass, nf=3, nbf=5)

    def test_resident_two_flushes_two_threads(self):
        lib = self.lib
        ptr_d = POINTER(c_double)
        ptr_i = POINTER(c_int)
        bindings = {
            'oqp_gpu_metc_session_begin': [c_int] * 5,
            'oqp_gpu_metc_session_zero_f3': [c_int],
            'oqp_gpu_metc_session_upload_d3': [c_int, ptr_d],
            'oqp_gpu_metc_session_download_f3': [c_int, ptr_d],
            'oqp_gpu_metc_session_contract': [c_int, c_int, ptr_i, ptr_d, c_int, c_int, c_double, c_double, c_int],
            'oqp_gpu_metc_session_end': [c_int],
        }
        for name, args in bindings.items():
            getattr(lib, name).argtypes = args
            getattr(lib, name).restype = c_int
        for mode in ('mrsf', 'umrsf'):
            nf, nbf, nm = 2, 4, NMATRIX[mode]
            ids, ints, f0, d3 = _deterministic_inputs(mode, nf, nbf)
            session = lib.oqp_gpu_metc_session_begin(nf, nm, nbf, 2, len(ints))
            self.assertGreater(session, 0)
            try:
                self.assertEqual(lib.oqp_gpu_metc_session_upload_d3(session, d3.ctypes.data_as(ptr_d)), 0)
                # Reuse the arena across passes; zero resets diagnostics too.
                for cur_pass in (1, 2):
                    self.assertEqual(lib.oqp_gpu_metc_session_zero_f3(session), 0)
                    ref = reference_contract(mode, cur_pass, ids, ints, np.zeros_like(f0), d3,
                                             nf, nm, nbf, .71, 1.13)
                    for thread in (0, 1):
                        for flush in range(thread + 1):
                            self.assertEqual(lib.oqp_gpu_metc_session_contract(
                                session, thread, ids.ctypes.data_as(ptr_i), ints.ctypes.data_as(ptr_d),
                                len(ints), cur_pass, .71, 1.13, int(mode == 'umrsf')), 0)
                    out = np.zeros(2 * len(f0), dtype=np.float64)
                    self.assertEqual(lib.oqp_gpu_metc_session_download_f3(session, out.ctypes.data_as(ptr_d)), 0)
                    np.testing.assert_allclose(out[:len(f0)], ref, rtol=1e-9, atol=1e-11)
                    np.testing.assert_allclose(out[len(f0):], 2*ref, rtol=1e-9, atol=1e-11)
            finally:
                self.assertEqual(lib.oqp_gpu_metc_session_end(session), 0)

if __name__ == "__main__":
    unittest.main()
