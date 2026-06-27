"""Validate the dftd4 dispersion contribution to the analytic Hessian.

The native electronic analytic Hessian (oqp.hf_hessian) does not include the
empirical dftd4 dispersion term.  The numerical Hessian includes it implicitly,
because each displaced gradient is dispersion-corrected.  To keep the two paths
consistent, the analytic Hessian driver adds d2 E_disp / dR2, obtained by central
finite difference of the (cheap, analytic) dftd4 dispersion gradient at the same
step the numerical Hessian uses (Hessian._dispersion_hessian).

This test runs H2O / RKS-PBE / cc-pVDZ with d4=true and checks that:
  * the analytic Hessian (electronic + dispersion) matches the OpenQP numerical
    Hessian (which also includes dispersion) to the finite-difference level;
  * the dispersion term is actually non-zero (guards against silently adding 0).

`pbe` is used because both OpenQP and dftd4 recognize that functional name (OpenQP
names like `b3lyp5` are not known to dftd4's DampingParam, which would fail the
whole d4 path -- energy, gradient and Hessian alike).

Skipped unless both the compiled OpenQP runtime and the dftd4 package are available.
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
functional=pbe
d4={d4}
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
        # dftd4 is linked natively into liboqp (no Python dftd4 package needed).
        return hasattr(oqp.lib, "oqp_dftd4_disp")
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(),
                     "compiled OpenQP runtime with native dftd4 not available")
class D4AnalyticHessian(unittest.TestCase):
    def _hessian(self, hess_type, d4):
        from oqp.pyoqp import Runner

        workdir = Path(f"/tmp/oqp_d4_hess_{hess_type}_{d4}")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "h2o.inp"
        inp.write_text(INPUT_TMPL.format(hess_type=hess_type, d4=d4))
        runner = Runner(project=f"h2o_{hess_type}_{d4}", input_file=str(inp),
                        log=str(workdir / "h2o.log"))
        runner.run()
        return np.array(runner.mol.hessian)

    def test_analytic_d4_matches_numerical_d4(self):
        han = self._hessian("analytical", "true")
        hnum = self._hessian("numerical", "true")
        self.assertEqual(han.shape, hnum.shape)
        # analytic Hessian (electronic + symmetrized dispersion FD) stays symmetric
        self.assertLess(np.max(np.abs(han - han.T)), 1.0e-9)
        # analytic + D4 must match numerical + D4 at the FD-truncation level
        self.assertLess(np.max(np.abs(han - hnum)), 5.0e-4)

    def test_dispersion_term_is_nonzero(self):
        han_no = self._hessian("analytical", "false")
        han_d4 = self._hessian("analytical", "true")
        # the dispersion second derivative must actually change the Hessian
        self.assertGreater(np.max(np.abs(han_d4 - han_no)), 1.0e-6)


if __name__ == "__main__":
    unittest.main()
