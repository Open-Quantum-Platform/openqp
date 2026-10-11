"""Native tests for the TDHF/TDDFT Davidson residual step.

Drives ``oqp_rpa_residual_selftest`` (tests/fortran/rpa_residual_selftest.F90),
which calls the production ``rparesvec`` and ``rpa_max_residual`` from
``source/tdhf_lib.F90``. A root whose squared residual norm exceeds 1 must keep
its true error, receive a preconditioned correction, and keep the Davidson
iteration unconverged.
"""

import os
import unittest

import numpy as np


NATIVE_TESTS_REQUIRED = os.getenv("OQP_REQUIRE_NATIVE_TESTS") == "1"
_RUNTIME_ERROR = ""


def _runtime_available():
    global _RUNTIME_ERROR
    try:
        # The installed oqp locates its own runtime; the source tree is not put
        # on any search path (tests/test_oqp_module_isolation.py).
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp

        if not hasattr(oqp.lib, "oqp_rpa_residual_selftest"):
            _RUNTIME_ERROR = "missing native symbol: oqp_rpa_residual_selftest"
            return False
        return True
    except Exception as exc:
        _RUNTIME_ERROR = f"{type(exc).__name__}: {exc}"
        return False


RUNTIME_AVAILABLE = _runtime_available()


class NativeRpaResidualBuildGate(unittest.TestCase):
    @unittest.skipUnless(
        NATIVE_TESTS_REQUIRED,
        "source-only run; set OQP_REQUIRE_NATIVE_TESTS=1 after building OpenQP",
    )
    def test_required_runtime_and_symbol(self):
        self.assertTrue(
            RUNTIME_AVAILABLE,
            f"compiled OpenQP runtime is required: {_RUNTIME_ERROR}",
        )


def _preconditioned(residual, energy, abd):
    # Test data keep every |energy - abd| far above the denominator floor.
    return residual / (energy - abd)


@unittest.skipUnless(RUNTIME_AVAILABLE, "compiled OpenQP runtime not available")
class RpaResidualTests(unittest.TestCase):
    TOL = 1.0e-5

    @staticmethod
    def _ptr(ffi, array):
        return ffi.cast("double *", int(array.ctypes.data))

    def run_resvec(self, w_l, w_r, v_l, v_r, ee, abd, imax, tamm_dancoff):
        import oqp

        ffi = oqp.ffi
        n, ndsr = w_r.shape
        w_l = np.asfortranarray(w_l, dtype=np.float64).copy(order="F")
        w_r = np.asfortranarray(w_r, dtype=np.float64).copy(order="F")
        v_l = np.asfortranarray(v_l, dtype=np.float64).copy(order="F")
        v_r = np.asfortranarray(v_r, dtype=np.float64).copy(order="F")
        ee = np.ascontiguousarray(ee, dtype=np.float64)
        abd = np.ascontiguousarray(abd, dtype=np.float64)
        errors = np.full(ndsr, np.nan)
        q = np.full((n, ndsr, 2), np.nan, order="F")
        mxerr = ffi.new("double *")
        converged = ffi.new("int *")
        oqp.lib.oqp_rpa_residual_selftest(
            n, ndsr, imax, bool(tamm_dancoff), self.TOL,
            self._ptr(ffi, w_l), self._ptr(ffi, w_r),
            self._ptr(ffi, v_l), self._ptr(ffi, v_r),
            self._ptr(ffi, ee), self._ptr(ffi, abd),
            self._ptr(ffi, errors), self._ptr(ffi, q),
            mxerr, converged,
        )
        return errors, q, float(mxerr[0]), int(converged[0])

    @staticmethod
    def _tda_problem():
        """Non-symmetric A; root 1 exact, root 2 ||r||^2 > 1, root 3 small."""
        rng = np.random.default_rng(20261011)
        n = 7
        abd = np.array([0.31, 0.42, 0.55, 0.71, 0.93, 1.20, 1.55])
        coupling = 0.9 * rng.standard_normal((n, n))
        a = np.diag(abd) + coupling
        evals, evecs = np.linalg.eig(a)
        k = int(np.argmin(np.abs(evals.imag) + 10.0 * (evals.real < 0)))
        e1 = float(evals[k].real)
        v1 = evecs[:, k].real / np.linalg.norm(evecs[:, k].real)
        v2 = rng.standard_normal(n)
        v2 /= np.linalg.norm(v2)
        v3 = np.zeros(n)
        v3[0] = 1.0
        a3 = np.diag(abd).copy()
        a3[1, 0] = 1.0e-2  # residual 2e-4 > tol for root 3
        v = np.column_stack([v1, v2, v3])
        w = np.column_stack([a @ v1, a @ v2, a3 @ v3])
        # Root 2 and 3 energies stay > 0.01 from every abd entry.
        ee = np.array([e1, 5.0, 0.30])
        return w, v, ee, abd

    def test_tda_large_residual_root_stays_unconverged(self):
        w, v, ee, abd = self._tda_problem()
        residual = w - v * ee
        expected = np.sum(residual**2, axis=0)
        self.assertGreater(expected[1], 1.0)
        self.assertLess(expected[0], self.TOL)

        errors, q, mxerr, converged = self.run_resvec(
            w, w, v, v, ee, abd, imax=0, tamm_dancoff=True)

        np.testing.assert_allclose(errors, expected, rtol=1e-13, atol=1e-28)
        self.assertEqual(mxerr, float(np.max(errors)))
        self.assertGreater(mxerr, 1.0)
        self.assertEqual(converged, 0)
        # Root 1 is converged: no correction. Roots 2 and 3 get the
        # preconditioned residual, including the large-residual root.
        np.testing.assert_array_equal(q[:, 0, 0], 0.0)
        for ivec in (1, 2):
            np.testing.assert_allclose(
                q[:, ivec, 0], _preconditioned(residual[:, ivec], ee[ivec], abd),
                rtol=1e-13, atol=0.0)
        self.assertGreater(np.linalg.norm(q[:, 1, 0]), 0.0)

    def test_tda_only_large_residual_root_cannot_converge(self):
        """The case the old clamp turned into a converged result."""
        w, v, ee, abd = self._tda_problem()
        w, v, ee = w[:, :2], v[:, :2], ee[:2]
        errors, q, mxerr, converged = self.run_resvec(
            w, w, v, v, ee, abd, imax=0, tamm_dancoff=True)
        self.assertGreater(errors[1], 1.0)
        self.assertEqual(converged, 0)
        self.assertTrue(np.all(np.isfinite(q[:, 1, 0])))
        self.assertGreater(np.linalg.norm(q[:, 1, 0]), 0.0)

    def test_rpa_left_residual_large_right_small(self):
        """max(left, right) is the error; q(:,:,1) is left, q(:,:,2) right."""
        rng = np.random.default_rng(7)
        n, ndsr = 6, 2
        abd = np.array([0.25, 0.37, 0.48, 0.66, 0.89, 1.30])
        v_l = rng.standard_normal((n, ndsr))
        v_r = rng.standard_normal((n, ndsr))
        ee = np.array([0.05, 0.12])
        r_l = 0.8 * rng.standard_normal((n, ndsr))  # ||r_l||^2 > 1
        r_r = 1.0e-2 * rng.standard_normal((n, ndsr))
        w_l = r_l + v_l * ee
        w_r = r_r + v_r * ee
        er_l = np.sum(r_l**2, axis=0)
        er_r = np.sum(r_r**2, axis=0)
        self.assertTrue(np.all(er_l > 1.0))
        self.assertTrue(np.all(er_r < 1.0))

        errors, q, mxerr, converged = self.run_resvec(
            w_l, w_r, v_l, v_r, ee, abd, imax=0, tamm_dancoff=False)

        np.testing.assert_allclose(errors, np.maximum(er_l, er_r), rtol=1e-12)
        self.assertEqual(converged, 0)
        for ivec in range(ndsr):
            np.testing.assert_allclose(
                q[:, ivec, 0], _preconditioned(r_l[:, ivec], ee[ivec], abd),
                rtol=1e-12, atol=1e-15)
            np.testing.assert_allclose(
                q[:, ivec, 1], _preconditioned(r_r[:, ivec], ee[ivec], abd),
                rtol=1e-12, atol=1e-15)

    def test_converged_roots_and_negative_roots_excluded(self):
        """imax roots are skipped; exact roots below tol converge."""
        n = 4
        abd = np.array([0.2, 0.4, 0.6, 0.8])
        v = np.eye(n)[:, :3]
        ee = abd[:3].copy()
        w = v * ee
        w[:, 0] += 5.0  # root 1 is a negative-eigenvalue root (imax=1)
        errors, q, mxerr, converged = self.run_resvec(
            w, w, v, v, ee, abd, imax=1, tamm_dancoff=True)
        np.testing.assert_array_equal(errors, 0.0)
        self.assertEqual(mxerr, 0.0)
        self.assertEqual(converged, 1)
        np.testing.assert_array_equal(q[:, 1:, 0], 0.0)


if __name__ == "__main__":
    unittest.main()
