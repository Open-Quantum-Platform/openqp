"""Validate the native analytic Hessian for range-separated (CAM/LC) functionals.

CAM-type functionals split the 1/r12 operator into a long-range part (full
Coulomb J + exchange K[1/r] scaled by cam_alpha) and a short-range part
(erfc-attenuated exchange K[erf(mu r)/r] scaled by cam_beta).  The OpenQP 2e
derivative-integral engine is erfc-attenuation capable, so every Fock build in
the analytic Hessian runs the same two-pass split keyed off infos%dft%cam_flag:
  * the 2e second-derivative skeleton -- grd2_hess_driver (two-pass wrapper);
  * the CPHF response-Fock derivative -- fock_deriv_contract -> grd2_driver
    (which already auto-detects cam_flag);
  * the reorthonormalization / relaxed-density Fock and the CPHF A-matrix --
    fock_jk / cphf_apbx (already CAM-aware via int2_driver).
The DFT exchange-correlation skeleton/kernel terms are the functional's own
(handled by dftexcor / the f_xc kernel) and are unchanged by the range split.

Checked against the OpenQP numerical Hessian (central FD of the analytic CAM
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
charge={charge}
runtype=hess
basis=6-31g
functional=cam-b3lyp
method=hf
[guess]
type=huckel
[scf]
multiplicity={mult}
type={scf}
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
class CamAnalyticHessian(unittest.TestCase):
    def _hessian(self, tag, charge, mult, scf, hess_type):
        from oqp.pyoqp import Runner

        workdir = Path(f"/tmp/oqp_cam_hess_{tag}_{hess_type}")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "oh.inp"
        inp.write_text(INPUT_TMPL.format(charge=charge, mult=mult, scf=scf,
                                         hess_type=hess_type))
        runner = Runner(project=f"oh_{tag}_{hess_type}", input_file=str(inp),
                        log=str(workdir / "oh.log"))
        runner.run()
        return np.array(runner.mol.hessian)

    def _check(self, tag, charge, mult, scf):
        han = self._hessian(tag, charge, mult, scf, "analytical")
        hnum = self._hessian(tag, charge, mult, scf, "numerical")
        self.assertEqual(han.shape, hnum.shape)
        self.assertLess(np.max(np.abs(han - han.T)), 1.0e-9)
        self.assertLess(np.max(np.abs(han - hnum)), 5.0e-4)

    def test_rks_cam_analytic_matches_numerical(self):
        # closed-shell OH+ (RKS/CAM-B3LYP), analytic CPHF response path
        self._check("rks", 1, 1, "rhf")

    def test_uks_cam_analytic_matches_numerical(self):
        # open-shell OH (UKS/CAM-B3LYP), per-spin response
        self._check("uks", 0, 2, "uhf")

    def test_roks_cam_analytic_matches_numerical(self):
        # open-shell OH (ROKS/CAM-B3LYP), semi-numerical response
        self._check("roks", 0, 2, "rohf")


if __name__ == "__main__":
    unittest.main()
