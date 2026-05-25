import re
import unittest
from pathlib import Path


SOURCE = Path(__file__).resolve().parents[1] / "source" / "integrals" / "int2.F90"


class Int2OpenMPWorkshareTest(unittest.TestCase):
    def test_shell_pair_workshare_has_barrier_between_dynamic_do_instances(self):
        text = SOURCE.read_text()
        match = re.search(
            r"subroutine int2_twoei\(this, int2_consumer\)(.*?)end subroutine int2_twoei",
            text,
            re.IGNORECASE | re.DOTALL,
        )
        if match is None:
            self.fail("Could not find int2_twoei")
        body = match.group(1)
        self.assertIn("!$omp do schedule(dynamic,2)", body)
        self.assertNotIn("!$omp end do nowait", body)
        self.assertIn("!$omp end do", body)


if __name__ == "__main__":
    unittest.main()
