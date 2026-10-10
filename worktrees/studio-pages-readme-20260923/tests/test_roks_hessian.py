"""Validate the native open-shell (ROKS) DFT analytic Hessian.

ROKS evaluates the HF response semi-numerically (resp_grad: central FD of the
analytic open-shell gradient along the CPHF-relaxed orbital path).  For DFT the
XC is folded into that same finite difference, with NO separate analytic XC
energy-weighting term (which would double-count):
  * the spin Fock used to build the energy-weighted density W' is the full KS
    Fock (Hcore + J - c_x K + Vxc), Vxc from the open-shell dftexcor at the
    displaced geometry/orbitals;
  * the explicit open-shell XC gradient (derexc_blk, urohf) at the displaced
    density is added to the gradient.
The CPHF right-hand side carries the XC skeleton dVxc/dR + f_xc[d0] (FD of the
spin XC Fock), and the CPHF operator already includes the open-shell f_xc kernel.

Checked against the OpenQP numerical Hessian (central FD of the analytic ROKS
gradient): they must agree to the finite-difference truncation level and the
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
   8   0.000000000   0.000000000   0.000000000
   1   0.000000000   0.000000000   0.970000000
charge=0
runtype=hess
basis=6-31g
functional=bhhlyp
method=hf
[guess]
type=huckel
[scf]
multiplicity=2
type=rohf
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
class RoksAnalyticHessian(unittest.TestCase):
    def _hessian(self, hess_type):
        from oqp.pyoqp import Runner

        workdir = Path(f"/tmp/oqp_roks_hess_{hess_type}")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "oh.inp"
        inp.write_text(INPUT_TMPL.format(hess_type=hess_type))
        runner = Runner(project=f"ohroks_{hess_type}", input_file=str(inp),
                        log=str(workdir / "oh.log"))
        runner.run()
        return np.array(runner.mol.hessian)

    def test_roks_analytic_matches_numerical(self):
        han = self._hessian("analytical")
        hnum = self._hessian("numerical")
        self.assertEqual(han.shape, hnum.shape)
        self.assertLess(np.max(np.abs(han - han.T)), 1.0e-9)
        self.assertLess(np.max(np.abs(han - hnum)), 5.0e-4)


if __name__ == "__main__":
    unittest.main()
