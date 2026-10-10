"""Validate the native open-shell (UKS) DFT analytic Hessian.

The UKS Hessian extends the validated UHF analytic assembly with the
exchange-correlation contribution, mirroring the closed-shell RKS structure for
two spins:
  * RHS:  central FD of the spin XC Fock (open-shell dftexcor) along the combined
          geometry + reorthonormalization path -> the XC skeleton dVxc/dR and
          f_xc[d0], subtracted from the CPHF right-hand side;
  * dHse: XC skeleton + density response, central FD of the analytic open-shell
          XC gradient (derexc_blk, urohf) along R +/- h, P^s +/- h dP^s;
  * dHt3: -sum_s sum_kl s1oo^s,x_kl (vxc^s,y + fxc[dP^s,y])_kl, the XC part of the
          energy-weighted term (coefficient -1 per spin, the UHF occupation).
The HF-exchange fraction is carried by hfscale in the Coulomb/exchange terms; the
CPHF operator already includes the open-shell f_xc kernel (utddft_fxc).

Checked against the OpenQP numerical Hessian (central FD of the analytic UKS
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
type=uhf
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
class UksAnalyticHessian(unittest.TestCase):
    def _hessian(self, hess_type):
        from oqp.pyoqp import Runner

        workdir = Path(f"/tmp/oqp_uks_hess_{hess_type}")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "oh.inp"
        inp.write_text(INPUT_TMPL.format(hess_type=hess_type))
        runner = Runner(project=f"ohuks_{hess_type}", input_file=str(inp),
                        log=str(workdir / "oh.log"))
        runner.run()
        return np.array(runner.mol.hessian)

    def test_uks_analytic_matches_numerical(self):
        han = self._hessian("analytical")
        hnum = self._hessian("numerical")
        self.assertEqual(han.shape, hnum.shape)
        self.assertLess(np.max(np.abs(han - han.T)), 1.0e-9)
        self.assertLess(np.max(np.abs(han - hnum)), 5.0e-4)


if __name__ == "__main__":
    unittest.main()
