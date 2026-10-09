"""hess.type=auto with an empty virtual space.

He2/STO-3G fills its basis with occupied orbitals, so the analytic kernel has
no response space and stores no Hessian.  The default (auto) must therefore
resolve to the numerical Hessian for it, while an ordinary molecule keeps the
analytic one.

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
    def _hessian_type(self, tmp, name, system):
        from oqp.pyoqp import Runner
        from oqp.library.single_point import Hessian
        inp = Path(tmp) / f"{name}.inp"
        inp.write_text(INPUT.format(system=system))
        runner = Runner(project=name, input_file=str(inp),
                        log=str(Path(tmp) / f"{name}.log"), silent=1, usempi=False)
        runner.run()
        runner.mol.config["hess"]["type"] = "auto"
        hess = Hessian(runner.mol)
        return hess.hess_type, hess.hess_type_reason

    def test_empty_virtual_space_uses_the_numerical_hessian(self):
        with tempfile.TemporaryDirectory() as tmp:
            kind, reason = self._hessian_type(tmp, "he2", HE2)
            self.assertEqual(kind, "numerical")
            self.assertIn("no virtual orbitals", reason)
            kind, _ = self._hessian_type(tmp, "h2o", H2O)
            self.assertEqual(kind, "analytical")


if __name__ == "__main__":
    unittest.main()
