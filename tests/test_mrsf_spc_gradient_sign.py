from pathlib import Path
import re
import unittest


class MrsfSpcGradientSignTests(unittest.TestCase):
    def test_mrsf_ovov_intra_square_uses_same_sign_as_coco(self):
        """Eq. 3.24 gives CO-CO and OV-OV intra-square SPC terms the same sign."""
        src = Path(__file__).resolve().parents[1] / "source/modules/tdhf_mrsf_gradient.F90"
        text = src.read_text()

        coco = re.search(r"df1\s*=\s*df1\s*([+-])\s*sgnk\s*\*\s*qfspcp1\s*\*\s*db1", text)
        ovov = re.search(r"df1\s*=\s*df1\s*([+-])\s*sgnk\s*\*\s*qfspcp2\s*\*\s*db2", text)

        self.assertIsNotNone(coco, "CO-CO SPC gradient term not found")
        self.assertIsNotNone(ovov, "OV-OV SPC gradient term not found")
        assert coco is not None
        assert ovov is not None
        self.assertEqual(coco.group(1), "+")
        self.assertEqual(ovov.group(1), coco.group(1))


if __name__ == "__main__":
    unittest.main()
