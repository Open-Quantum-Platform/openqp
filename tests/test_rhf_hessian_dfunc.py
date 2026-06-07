"""Regression test: closed-shell RHF analytic Hessian with d functions.

The closed-shell HF/DFT analytic Hessian (hf_hessian_mod::hf_hessian) brings the
derivative integrals dS/dT/dV into the OpenQP normalized (bfnrm) convention before
the CPHF RHS and response contractions.  That scaling must be applied EXACTLY ONCE:
for basis functions with bfnrm /= 1 (d/f shells, e.g. cc-pVDZ) a double application
scales the integrals by bfnrm**2 and corrupts the Hessian, while staying invisible
for s/p-only bases (bfnrm == 1, e.g. STO-3G / 6-31G).

This guards specifically against that failure mode (a double-bfnrm regression that
once slipped in via a semantic merge): H2O / RHF / cc-pVDZ analytic vs the OpenQP
numerical Hessian must agree to the finite-difference truncation level, and the
analytic matrix must be exactly symmetric.

Skipped unless the compiled OpenQP runtime is importable.
"""

import os
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

INPUT_TMPL = """[input]
system=
   8   0.000000000   0.000000000   0.117790000
   1   0.000000000   0.755450000  -0.471160000
   1   0.000000000  -0.755450000  -0.471160000
charge=0
runtype=hess
basis=cc-pvdz
method=hf
[guess]
type=huckel
[scf]
multiplicity=1
type=rhf
[hess]
type={hess_type}
state=0
dx=0.005
nproc=1
"""


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class RhfDFunctionAnalyticHessian(unittest.TestCase):
    def _hessian(self, hess_type):
        from oqp.pyoqp import Runner

        workdir = Path(f"/tmp/oqp_rhf_dfunc_hess_{hess_type}")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "h2o.inp"
        inp.write_text(INPUT_TMPL.format(hess_type=hess_type))
        runner = Runner(project=f"h2o_{hess_type}", input_file=str(inp),
                        log=str(workdir / "h2o.log"))
        runner.run()
        return np.array(runner.mol.hessian)

    def test_rhf_dfunction_analytic_matches_numerical(self):
        han = self._hessian("analytical")
        hnum = self._hessian("numerical")
        self.assertEqual(han.shape, hnum.shape)
        # analytic Hessian must be (essentially) exactly symmetric
        self.assertLess(np.max(np.abs(han - han.T)), 1.0e-9)
        # and match the numerical reference at the FD-truncation level; a
        # double-bfnrm regression drives this to O(1).
        self.assertLess(np.max(np.abs(han - hnum)), 5.0e-4)


if __name__ == "__main__":
    unittest.main()
