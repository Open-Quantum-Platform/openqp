"""The analytic Hessian's L<=3 limit for custom ``file:`` bases.

A ``file:`` basis is read relative to the input directory, exactly as
set_basis does.  ``hess.type=auto`` then picks the analytic Hessian when the
custom basis stays within L<=3, and an explicit ``hess.type=analytical`` is
rejected before the native run when the basis has g functions or cannot be
read at all.
"""
import importlib.util
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_checker():
    # Stub the runtime package only while the module loads, then restore
    # sys.modules so later tests in the same process import the real oqp.
    saved = {key: sys.modules.get(key) for key in ("oqp", "oqp.utils", "oqp.utils.mpi_utils")}
    try:
        return _load_checker()
    finally:
        for key, value in saved.items():
            if value is None:
                sys.modules.pop(key, None)
            else:
                sys.modules[key] = value


def _load_checker():
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        size = 1
        use_mpi = False

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils
    name = "input_checker_hessian_custom_basis_under_test"
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "pyoqp/oqp/utils/input_checker.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def shell(l, exps, coefs):
    return {"function_type": "gto", "region": "", "angular_momentum": [l],
            "exponents": [str(e) for e in exps], "coefficients": [[str(c) for c in coefs]]}


def write_basis(directory, name, max_l_h):
    shells_h = [shell(0, [1.0], [1.0])] + [shell(l, [0.8], [1.0]) for l in range(1, max_l_h + 1)]
    data = {"molssi_bse_schema": {"schema_type": "complete", "schema_version": "0.1"},
            "elements": {"1": {"electron_shells": shells_h}}}
    (Path(directory) / name).write_text(json.dumps(data))


def config(basis, hess_type):
    return {
        "input": {"runtype": "hess", "method": "hf", "basis": basis,
                  "system": "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74"},
        "guess": {}, "scf": {"type": "rhf", "multiplicity": 1},
        "tdhf": {}, "properties": {}, "hess": {"type": hess_type, "state": 0},
    }


def basis_errors(report):
    return [d for d in report.diagnostics
            if d.severity == "ERROR" and d.path == "input.basis"]


class CustomBasisHessianCheck(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.checker = load_checker()

    def test_auto_reads_a_custom_basis_within_the_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_basis(tmp, "spd.json", 2)
            cfg = config("file:spd.json", "auto")
            self.assertEqual(self.checker.resolve_hessian_type(cfg, tmp)[0], "analytical")
            # without the input directory the file cannot be found: numerical
            self.assertEqual(self.checker.resolve_hessian_type(cfg)[0], "numerical")

    def test_explicit_analytical_rejects_g_functions_in_a_custom_basis(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_basis(tmp, "spdfg.json", 4)
            report = self.checker.check_input_values(
                config("file:spdfg.json", "analytical"),
                raise_error=False, emit=False, input_dir=tmp)
            self.assertTrue(any("L=4" in str(d.value) for d in basis_errors(report)),
                            report.to_text())

    def test_explicit_analytical_rejects_an_unreadable_custom_basis(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = self.checker.check_input_values(
                config("file:missing.json", "analytical"),
                raise_error=False, emit=False, input_dir=tmp)
            self.assertTrue(basis_errors(report), report.to_text())

    def test_explicit_analytical_accepts_a_readable_custom_basis(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_basis(tmp, "spd.json", 2)
            report = self.checker.check_input_values(
                config("file:spd.json", "analytical"),
                raise_error=False, emit=False, input_dir=tmp)
            self.assertFalse(basis_errors(report), report.to_text())


if __name__ == "__main__":
    unittest.main()
