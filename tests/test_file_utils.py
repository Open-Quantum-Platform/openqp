import importlib.util
import sys
import types
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_file_utils():
    stub_names = (
        "oqp",
        "oqp.molden",
        "oqp.molden.moldenwriter",
        "oqp.periodic_table",
        "oqp.runtime",
        "oqp.utils",
        "oqp.utils.constants",
        "oqp.utils.mpi_utils",
    )
    saved_modules = {name: sys.modules.get(name) for name in stub_names}

    try:
        oqp = types.ModuleType("oqp")
        oqp.__path__ = []
        sys.modules["oqp"] = oqp

        molden = types.ModuleType("oqp.molden")
        molden.__path__ = []
        sys.modules["oqp.molden"] = molden

        moldenwriter = types.ModuleType("oqp.molden.moldenwriter")
        setattr(moldenwriter, "write_frequency", lambda *args, **kwargs: "")
        sys.modules["oqp.molden.moldenwriter"] = moldenwriter

        periodic_table = types.ModuleType("oqp.periodic_table")
        setattr(periodic_table, "SYMBOL_MAP", {1: 1, 8: 8, "1": 1, "8": 8})
        elements = [""] * 9
        elements[1] = "H"
        elements[8] = "O"
        setattr(periodic_table, "ELEMENTS_NAME", elements)
        sys.modules["oqp.periodic_table"] = periodic_table

        runtime = types.ModuleType("oqp.runtime")
        setattr(runtime, "basis_search_paths", lambda: [])
        sys.modules["oqp.runtime"] = runtime

        utils = types.ModuleType("oqp.utils")
        utils.__path__ = []
        sys.modules["oqp.utils"] = utils

        constants = types.ModuleType("oqp.utils.constants")
        setattr(constants, "ANGSTROM_TO_BOHR", 0.529177210903)
        sys.modules["oqp.utils.constants"] = constants

        mpi_utils = types.ModuleType("oqp.utils.mpi_utils")
        setattr(mpi_utils, "mpi_dump", lambda func: func)
        sys.modules["oqp.utils.mpi_utils"] = mpi_utils

        spec = importlib.util.spec_from_file_location(
            "file_utils_under_test",
            ROOT / "pyoqp/oqp/utils/file_utils.py",
        )
        assert spec is not None
        assert spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        sys.modules["file_utils_under_test"] = module
        spec.loader.exec_module(module)
        return module
    finally:
        for name, module in saved_modules.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


class TestWriteXYZ(unittest.TestCase):
    def test_accepts_numpy_2d_atomic_number_arrays(self):
        file_utils = load_file_utils()
        atoms = np.array([[8], [1], [1]], dtype=np.int64)
        coord = np.array(
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.8897261254576558,
                0.0,
                1.8897261254576558,
                0.0,
            ],
            dtype=float,
        )

        xyz = file_utils.write_xyz(atoms, coord, (3, -76.0))

        lines = xyz.splitlines()
        self.assertEqual(lines[0], "3")
        self.assertEqual(lines[1], "Geom 3 -76.0")
        self.assertTrue(lines[2].startswith("O"), lines[2])
        self.assertTrue(lines[3].startswith("H"), lines[3])
        self.assertTrue(lines[4].startswith("H"), lines[4])


if __name__ == "__main__":
    unittest.main()
