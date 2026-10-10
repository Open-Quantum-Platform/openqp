"""Validate the native analytic Hessian with effective core potentials (ECP).

ECP second derivatives are supplied by libecpint (derivative order 2):
  * the ECP skeleton d^2 V_ECP/dR^2 (fixed density) is contracted analytically
    via ecp_tool::add_ecphess into the Cartesian Hessian, packed by the libecpint
    atom-pair convention {AA, AB, AC, ..., BB, ...};
  * the ECP first-derivative integrals enter the core-Hamiltonian derivative
    dHcore/dR (ecp_tool::ecp_deriv_ints added into the nuclear-attraction
    derivative tensor), so the ECP also drives the CPHF right-hand side and the
    orbital-relaxation response, exactly as point-charge nuclear attraction does;
  * the point-charge nuclear attraction uses the ECP-screened charge
    (zn - ecp_zn_num) consistently in BOTH the skeleton (hess_en/hess_nn) and the
    response (der_nucattr), the latter being a latent bug exposed by ECP.
For ROHF/ROKS (semi-numerical response) the ECP gradient (add_ecpder) is folded
into resp_grad with the ECP centre moved in lockstep with its atom.

Checked against the OpenQP numerical Hessian (central FD of the analytic ECP
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
   35   0.000000000   0.000000000   0.000000000
   1    0.000000000   0.000000000   1.407611000
charge={charge}
runtype=hess
basis=LANL2DZ
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
class EcpAnalyticHessian(unittest.TestCase):
    def _hessian(self, tag, charge, mult, scf, hess_type):
        from oqp.pyoqp import Runner

        workdir = Path(f"/tmp/oqp_ecp_hess_{tag}_{hess_type}")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "hbr.inp"
        inp.write_text(INPUT_TMPL.format(charge=charge, mult=mult, scf=scf,
                                         hess_type=hess_type))
        runner = Runner(project=f"hbr_{tag}_{hess_type}", input_file=str(inp),
                        log=str(workdir / "hbr.log"))
        runner.run()
        return np.array(runner.mol.hessian)

    def _check(self, tag, charge, mult, scf):
        han = self._hessian(tag, charge, mult, scf, "analytical")
        hnum = self._hessian(tag, charge, mult, scf, "numerical")
        self.assertEqual(han.shape, hnum.shape)
        self.assertLess(np.max(np.abs(han - han.T)), 1.0e-9)
        self.assertLess(np.max(np.abs(han - hnum)), 5.0e-4)

    def test_rhf_ecp_analytic_matches_numerical(self):
        # closed-shell HBr, analytic CPHF response path
        self._check("rhf", 0, 1, "rhf")

    def test_rohf_ecp_analytic_matches_numerical(self):
        # open-shell HBr+ doublet, semi-numerical (resp_grad) response path
        self._check("rohf", 1, 2, "rohf")


if __name__ == "__main__":
    unittest.main()
