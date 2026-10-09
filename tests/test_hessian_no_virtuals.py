"""hess.type=auto with an empty virtual space or a meta-GGA functional.

He2/STO-3G fills its basis with occupied orbitals, so the analytic kernel has
no response space and stores no Hessian.  The default (auto) must therefore
resolve to the numerical Hessian for it, while an ordinary molecule keeps the
analytic one.  A meta-GGA (LibXC tau-dependent family, here TPSS) also resolves
to numerical, since the analytic Hessian has no tau channel; B3LYP5 does not.

Skipped unless the compiled OpenQP runtime is importable.
"""

import os
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

INPUT = """[input]
system=
{system}
charge=0
runtype=energy
basis=sto-3g
method=hf
{dft}
[scf]
type=rhf
multiplicity=1
conv=1.0e-10
"""

HE2 = """   2  0.0 0.0 0.0
   2  0.0 0.0 3.0"""
H2O = """   8   0.000000000   0.000000000   0.117300000
   1   0.000000000   0.757200000  -0.469200000
   1   0.000000000  -0.757200000  -0.469200000"""


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class AutoHessianWithoutVirtuals(unittest.TestCase):
    def _hessian_type(self, tmp, name, system, dft=""):
        from oqp.pyoqp import Runner
        from oqp.library.single_point import Hessian
        inp = Path(tmp) / f"{name}.inp"
        inp.write_text(INPUT.format(system=system, dft=dft))
        runner = Runner(project=name, input_file=str(inp),
                        log=str(Path(tmp) / f"{name}.log"), silent=1, usempi=False)
        runner.run()
        runner.mol.config["hess"]["type"] = "auto"
        hess = Hessian(runner.mol)
        self.last_needs_tau = hess._functional_needs_tau()
        return hess.hess_type, hess.hess_type_reason

    def test_empty_virtual_space_uses_the_numerical_hessian(self):
        with tempfile.TemporaryDirectory() as tmp:
            kind, reason = self._hessian_type(tmp, "he2", HE2)
            self.assertEqual(kind, "numerical")
            self.assertIn("no virtual orbitals", reason)
            kind, _ = self._hessian_type(tmp, "h2o", H2O)
            self.assertEqual(kind, "analytical")

    def test_meta_gga_uses_the_numerical_hessian(self):
        import oqp
        if not hasattr(oqp, "oqp_functional_needs_tau"):
            self.skipTest("runtime predates oqp_functional_needs_tau")
        with tempfile.TemporaryDirectory() as tmp:
            kind, reason = self._hessian_type(tmp, "h2o_tpss", H2O, "functional=tpss")
            self.assertEqual(kind, "numerical")
            self.assertIn("meta-GGA", reason)
            # the LibXC family query itself (independent of the name set)
            self.assertTrue(self.last_needs_tau)
            kind, _ = self._hessian_type(tmp, "h2o_b3lyp", H2O, "functional=b3lyp5")
            self.assertEqual(kind, "analytical")
            self.assertFalse(self.last_needs_tau)


if __name__ == "__main__":
    unittest.main()
