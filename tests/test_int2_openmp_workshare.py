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
        self.assertIn("call int2_build_shell_pair_map(nshell, pair_i, pair_j)", body)
        self.assertIn("i = pair_i(ij_pair)", body)
        self.assertIn("j = pair_j(ij_pair)", body)
        self.assertNotIn("!$omp end do nowait", body)
        self.assertNotIn("!$omp do schedule(dynamic,2)", body)

        flat_workshare = re.search(
            r"do\s+ij_pair\s*=\s*1\s*,\s*npairs(?P<body>.*?)!\$omp\s+end\s+do",
            body,
            flags=re.DOTALL,
        )
        if flat_workshare is None:
            self.fail("flat shell-pair workshare not found")
        flat_body = flat_workshare.group("body")
        self.assertNotRegex(flat_body, r"do\s+j\s*=\s*1\s*,\s*i")
        self.assertNotIn("!$omp do", flat_body)
        self.assertNotIn("!$omp barrier", flat_body)
        self.assertNotIn("call int2_decode_shell_pair", flat_body)

    def test_precomputed_shell_pair_map_uses_stable_cost_classes(self):
        body = self._int2_twoei_source()

        build_map = re.search(
            r"subroutine\s+int2_build_shell_pair_map\b(?P<body>.*?)end\s+subroutine\s+int2_build_shell_pair_map",
            body,
            flags=re.DOTALL,
        )
        if build_map is None:
            self.fail("int2_build_shell_pair_map subroutine not found")
        map_body = build_map.group("body")

        self.assertIn("cnt = 0", map_body)
        self.assertIn("cost = sw * real(nbfi*nbfj,dp) * real(i,dp)", map_body)
        self.assertIn("off(0) = 0", map_body)
        self.assertRegex(
            map_body,
            r"do\s+c\s*=\s*1\s*,\s*nclass\s*\n\s*off\(c\)\s*=\s*off\(c-1\)\s*\+\s*cnt\(c-1\)",
        )
        self.assertRegex(
            map_body,
            r"do\s+i\s*=\s*nshell\s*,\s*1\s*,\s*-1\s*\n[\s\S]*?do\s+j\s*=\s*1\s*,\s*i",
        )
        self.assertRegex(
            map_body,
            r"off\(cls\)\s*=\s*off\(cls\)\s*\+\s*1\s*\n\s*ij_pair\s*=\s*off\(cls\)\s*\n\s*shell_pair_i\(ij_pair\)\s*=\s*i\s*\n\s*shell_pair_j\(ij_pair\)\s*=\s*j",
        )


if __name__ == "__main__":
    unittest.main()
