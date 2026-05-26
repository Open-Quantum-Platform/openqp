import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INT2 = ROOT / "source" / "integrals" / "int2.F90"


class Int2OpenMPWorkshareTests(unittest.TestCase):
    def _int2_twoei_source(self):
        text = INT2.read_text()
        match = re.search(
            r"subroutine\s+int2_twoei\b(?P<body>.*?)end\s+subroutine\s+int2_twoei",
            text,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if match is None:
            self.fail("int2_twoei subroutine not found")
        return match.group("body").lower()

    def test_shell_pair_parallelism_uses_single_flat_pair_workshare(self):
        body = self._int2_twoei_source()

        self.assertIn("npairs = nshell*(nshell+1)/2", body)
        self.assertRegex(
            body,
            r"!\$omp\s+do\s+schedule\(dynamic,1\)\s*\n\s*do\s+ij_pair\s*=\s*1\s*,\s*npairs",
        )
        self.assertIn("call int2_decode_shell_pair(ij_pair, nshell, i, j)", body)
        self.assertNotIn("!$omp end do nowait", body)
        self.assertNotIn("!$omp do schedule(dynamic,2)", body)

        flat_workshare = re.search(
            r"do\s+ij_pair\s*=\s*1\s*,\s*npairs(?P<body>.*?)!\$omp\s+end\s+do",
            body,
            flags=re.DOTALL,
        )
        if flat_workshare is None:
            self.fail("flat shell-pair workshare not found")
        self.assertNotRegex(flat_workshare.group("body"), r"do\s+j\s*=\s*1\s*,\s*i")
        self.assertNotIn("!$omp do", flat_workshare.group("body"))
        self.assertNotIn("!$omp barrier", flat_workshare.group("body"))

    def test_shell_pair_decoder_preserves_descending_i_then_ascending_j_order(self):
        body = self._int2_twoei_source()

        self.assertIn("subroutine int2_decode_shell_pair(ij_pair, nshell, i, j)", body)
        self.assertIn("remaining = ij_pair", body)
        self.assertIn("do i = nshell, 1, -1", body)
        self.assertIn("if (remaining <= i) then", body)
        self.assertIn("j = remaining", body)


if __name__ == "__main__":
    unittest.main()
