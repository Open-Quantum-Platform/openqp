"""Guard singular state-overlap minors against zero-pivot division."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class StateOverlapDeterminantTests(unittest.TestCase):
    def test_singular_minor_returns_zero_before_elimination(self):
        source = (
            ROOT / "source" / "modules" / "get_states_overlap.F90"
        ).read_text()
        body = source.split("function comp_det(array, n)", 1)[1].split(
            "end function comp_det", 1
        )[0]
        pivot_guard = body.index("abs(array(k,k)) <= tiny(1.0_dp)")
        division = body.index("array(m,k)/array(k,k)")
        self.assertLess(pivot_guard, division)
        self.assertIn("det = 0.0_dp", body[pivot_guard:division])
        self.assertIn("return", body[pivot_guard:division])


if __name__ == "__main__":
    unittest.main()
