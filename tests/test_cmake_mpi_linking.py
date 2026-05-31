import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class CMakeMPILinkingTests(unittest.TestCase):
    def test_oqp_target_links_mpi_fortran_when_mpi_enabled(self):
        cmake = (ROOT / "source" / "CMakeLists.txt").read_text()
        match = re.search(r"if\(ENABLE_MPI\)(.*?)endif\(\)", cmake, flags=re.S)

        self.assertIsNotNone(match, "source/CMakeLists.txt must contain an ENABLE_MPI block for oqp")
        if match is None:
            return
        block = match.group(1)
        self.assertIn("target_compile_definitions(oqp PRIVATE ENABLE_MPI)", block)
        self.assertRegex(
            block,
            r"target_link_libraries\s*\(\s*oqp\s+MPI::MPI_Fortran\s*\)",
            "ENABLE_MPI must link the oqp target against MPI::MPI_Fortran, not only define ENABLE_MPI",
        )


if __name__ == "__main__":
    unittest.main()
