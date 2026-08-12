"""Numerical tests for the deterministic native simplex-QP solver."""

import os
import struct
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
NATIVE_TESTS_REQUIRED = os.getenv("OQP_REQUIRE_NATIVE_TESTS") == "1"
_RUNTIME_ERROR = ""


def _runtime_available():
    global _RUNTIME_ERROR
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp

        missing = [
            name for name in (
                "oqp_simplex_qp_solve",
                "oqp_simplex_qp_solve_avoid",
            )
            if not hasattr(oqp.lib, name)
        ]
        if missing:
            _RUNTIME_ERROR = f"missing native symbols: {', '.join(missing)}"
            return False
        return True
    except Exception as exc:
        _RUNTIME_ERROR = f"{type(exc).__name__}: {exc}"
        return False


RUNTIME_AVAILABLE = _runtime_available()


class NativeSimplexQPBuildGate(unittest.TestCase):
    @unittest.skipUnless(
        NATIVE_TESTS_REQUIRED,
        "source-only run; set OQP_REQUIRE_NATIVE_TESTS=1 after building OpenQP",
    )
    def test_required_runtime_and_symbols(self):
        self.assertTrue(
            RUNTIME_AVAILABLE,
            f"compiled OpenQP runtime is required: {_RUNTIME_ERROR}",
        )


@unittest.skipUnless(RUNTIME_AVAILABLE, "compiled OpenQP runtime not available")
class SimplexQPTests(unittest.TestCase):
    @staticmethod
    def _pointer(ffi, array):
        return ffi.cast("double *", int(array.ctypes.data))

    def solve(self, h, g):
        import oqp

        h = np.asfortranarray(h, dtype=np.float64)
        g = np.ascontiguousarray(g, dtype=np.float64)
        self.assertEqual(h.shape, (g.size, g.size))
        x = np.empty(g.size, dtype=np.float64)
        value = oqp.ffi.new("double *")
        status = oqp.ffi.new("int *")
        oqp.lib.oqp_simplex_qp_solve(
            g.size,
            self._pointer(oqp.ffi, h),
            self._pointer(oqp.ffi, g),
            self._pointer(oqp.ffi, x),
            value,
            status,
        )
        return x, float(value[0]), int(status[0])

    def avoid(self, h, g, forbidden, forbid_vertices_before=0):
        import oqp

        h = np.asfortranarray(h, dtype=np.float64)
        g = np.ascontiguousarray(g, dtype=np.float64)
        forbidden = np.asfortranarray(forbidden, dtype=np.float64)
        self.assertEqual(h.shape, (g.size, g.size))
        self.assertEqual(forbidden.shape[0], g.size)
        x = np.empty(g.size, dtype=np.float64)
        value = oqp.ffi.new("double *")
        status = oqp.ffi.new("int *")
        oqp.lib.oqp_simplex_qp_solve_avoid(
            g.size,
            self._pointer(oqp.ffi, h),
            self._pointer(oqp.ffi, g),
            forbidden.shape[1],
            self._pointer(oqp.ffi, forbidden),
            forbid_vertices_before,
            self._pointer(oqp.ffi, x),
            value,
            status,
        )
        return x, float(value[0]), int(status[0])

    def assert_simplex(self, x):
        self.assertTrue(np.all(np.isfinite(x)))
        self.assertGreaterEqual(float(np.min(x)), -1.0e-14)
        self.assertAlmostEqual(float(np.sum(x)), 1.0, places=13)

    def test_convex_uniform_solution(self):
        x, value, status = self.solve(2.0*np.eye(3), np.zeros(3))
        self.assertEqual(status, 0)
        np.testing.assert_allclose(x, np.full(3, 1.0/3.0), atol=1.0e-13)
        self.assertAlmostEqual(value, 1.0/3.0, places=13)

    def test_convex_vertex_solution(self):
        x, _, status = self.solve(2.0*np.eye(2), [-10.0, 0.0])
        self.assertEqual(status, 0)
        np.testing.assert_allclose(x, [1.0, 0.0], atol=1.0e-13)

    def test_nonconvex_global_vertex(self):
        x, _, status = self.solve(np.diag([-2.0, -1.0]), np.zeros(2))
        self.assertEqual(status, 0)
        np.testing.assert_allclose(x, [1.0, 0.0], atol=1.0e-13)

    def test_singular_psd_face_prefers_the_newest_equivalent_state(self):
        h = np.array([[2.0, 2.0, 0.0],
                      [2.0, 2.0, 0.0],
                      [0.0, 0.0, 2.0]])
        x, _, status = self.solve(h, np.zeros(3))
        self.assertEqual(status, 0)
        # x3 is fixed at 0.5 at every minimum. The remaining 0.5 is placed on
        # the newer of the two equivalent old states by the documented tie rule.
        np.testing.assert_allclose(x, [0.0, 0.5, 0.5], atol=2.0e-13)

    def test_flat_problem_prefers_latest_state(self):
        x, _, status = self.solve(np.zeros((3, 3)), np.zeros(3))
        self.assertEqual(status, 0)
        np.testing.assert_array_equal(x, [0.0, 0.0, 1.0])

    def test_single_state(self):
        x, value, status = self.solve([[7.0]], [-2.0])
        self.assertEqual(status, 0)
        np.testing.assert_array_equal(x, [1.0])
        self.assertEqual(value, 1.5)

    def test_extreme_scaling_is_invariant(self):
        h = np.array([[3.0, -0.5, 0.0],
                      [-0.5, 2.0, 0.25],
                      [0.0, 0.25, 1.0]])
        g = np.array([-0.3, 0.2, -0.1])
        reference, _, ref_status = self.solve(h, g)
        self.assertEqual(ref_status, 0)
        for scale in (1.0e-150, 1.0e150):
            with self.subTest(scale=scale):
                x, _, status = self.solve(scale*h, scale*g)
                self.assertEqual(status, 0)
                np.testing.assert_allclose(x, reference, atol=2.0e-13)

    def test_large_common_linear_offset_does_not_change_solution(self):
        h = 2.0*np.eye(2)
        g = np.array([1.0e14 - 0.2, 1.0e14 + 0.2])
        expected_x1 = (2.0 - (g[0] - g[1]))/4.0

        x, _, status = self.solve(h, g)

        self.assertEqual(status, 0)
        np.testing.assert_allclose(x, [expected_x1, 1.0 - expected_x1], atol=1.0e-13)

    def test_nonfinite_input_returns_latest_vertex(self):
        g = np.array([0.0, np.nan, 0.0])
        x, value, status = self.solve(np.eye(3), g)
        self.assertEqual(status, 3)
        np.testing.assert_array_equal(x, [0.0, 0.0, 1.0])
        self.assertGreater(value, 1.0e300)

    def test_repeated_calls_are_bitwise_deterministic(self):
        h = np.array([[1.7, -0.3, 0.2],
                      [-0.3, 2.1, -0.4],
                      [0.2, -0.4, 0.9]])
        g = np.array([0.2, -0.1, 0.05])
        first = self.solve(h, g)
        signature = (first[0].tobytes(), struct.pack("=d", first[1]), first[2])
        for _ in range(10):
            current = self.solve(h, g)
            self.assertEqual(
                (current[0].tobytes(), struct.pack("=d", current[1]), current[2]),
                signature,
            )

    def test_ill_conditioned_thirteen_dimensional_solution_is_certified(self):
        n = 13
        i = np.arange(n)[:, None]
        k = np.arange(n)[None, :]
        q = np.sqrt(2.0 / n) * np.cos(np.pi * (i + 0.5) * k / n)
        q[:, 0] = 1.0 / np.sqrt(n)
        eigenvalues = np.geomspace(1.0e-3, 1.0, n)
        h = q @ np.diag(eigenvalues) @ q.T
        target = np.arange(1, n + 1, dtype=np.float64)
        target /= target.sum()
        g = -h @ target

        x, value, status = self.solve(h, g)
        self.assertEqual(status, 0)
        self.assert_simplex(x)
        np.testing.assert_allclose(x, target, atol=1.0e-9, rtol=0.0)

        gradient = h @ x + g
        lagrange = -float(np.mean(gradient[x > 1.0e-10]))
        kkt_residual = float(np.max(np.abs(gradient + lagrange)))
        value_star = 0.5 * target @ h @ target + g @ target
        self.assertLess(kkt_residual, 1.0e-10)
        self.assertLessEqual(value - value_star, 1.0e-12)

    def test_above_exact_dimension_limit_is_explicit_failure(self):
        n = 14
        g = np.linspace(-0.2, 0.3, n)
        x, _, status = self.solve(2.0*np.eye(n), g)
        self.assertEqual(status, 1)
        expected = np.zeros(n)
        expected[-1] = 1.0
        np.testing.assert_array_equal(x, expected)

    def test_flat_problem_uses_allowed_interior_when_preferred_is_forbidden(self):
        forbidden = np.array([[0.0], [1.0]])
        x, _, status = self.avoid(
            np.zeros((2, 2)), np.zeros(2), forbidden
        )
        self.assertEqual(status, 0)
        self.assert_simplex(x)
        self.assertGreaterEqual(np.linalg.norm(x - forbidden[:, 0]), 1.0e-4)

    def test_flat_singleton_reports_no_allowed_candidate(self):
        x, _, status = self.avoid(
            np.zeros((1, 1)), np.zeros(1), np.ones((1, 1))
        )
        self.assertEqual(status, 5)
        np.testing.assert_array_equal(x, [1.0])

    def test_flat_problem_uses_mixture_when_all_old_vertices_are_forbidden(self):
        x, _, status = self.avoid(
            np.zeros((2, 2)), np.zeros(2), np.empty((2, 0)),
            forbid_vertices_before=2,
        )
        self.assertEqual(status, 0)
        self.assert_simplex(x)
        self.assertLess(float(np.max(x)), 1.0 - 1.0e-4)

    def test_forbidden_solution_and_old_vertex_are_avoided(self):
        uniform = np.full(3, 1.0/3.0)
        x, _, status = self.avoid(
            2.0*np.eye(3), np.zeros(3), uniform.reshape(3, 1)
        )
        self.assertEqual(status, 0)
        self.assert_simplex(x)
        self.assertGreaterEqual(np.linalg.norm(x - uniform), 1.0e-4)

        old_vertex = np.array([1.0, 0.0, 0.0])
        x, _, status = self.avoid(
            np.zeros((3, 3)), [-10.0, 0.0, 0.0],
            np.empty((3, 0)), forbid_vertices_before=2,
        )
        self.assertEqual(status, 0)
        self.assert_simplex(x)
        self.assertGreaterEqual(np.linalg.norm(x - old_vertex), 1.0e-4)


if __name__ == "__main__":
    unittest.main()
