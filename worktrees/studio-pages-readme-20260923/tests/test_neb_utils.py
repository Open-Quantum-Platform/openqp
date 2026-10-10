import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class TestNEBEndpointInterpolation(unittest.TestCase):
    def setUp(self):
        self.neb_utils = load_module(
            "neb_utils_under_test",
            "pyoqp/oqp/library/neb_utils.py",
        )

    def test_interpolate_preserves_endpoints_and_linearly_spaces_images(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            reactant = Path(tmpdir) / "reactant.xyz"
            product = Path(tmpdir) / "product.xyz"
            reactant.write_text(
                "2\nreactant\nH 0.0 0.0 0.0\nH 0.0 0.0 1.0\n"
            )
            product.write_text(
                "2\nproduct\nH 0.0 0.0 0.2\nH 0.0 0.0 1.4\n"
            )

            images = self.neb_utils.interpolate_xyz_endpoints(
                reactant,
                product,
                nimage=3,
            )

        self.assertEqual([image.symbols for image in images], [["H", "H"]] * 3)
        self.assertEqual(images[0].coordinates_angstrom, [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        self.assertEqual(images[1].coordinates_angstrom, [[0.0, 0.0, 0.1], [0.0, 0.0, 1.2]])
        self.assertEqual(images[2].coordinates_angstrom, [[0.0, 0.0, 0.2], [0.0, 0.0, 1.4]])


if __name__ == "__main__":
    unittest.main()
