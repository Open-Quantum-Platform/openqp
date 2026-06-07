"""Validate the native open-shell (UHF) HF analytic Hessian.

The open-shell analytic Hessian (hf_hessian_mod::hf_hessian_uhf) assembles, per
spin and summed over s in {alpha, beta} with single occupation factors:

  * the fixed-density skeleton (1e total density + open-shell Lagrangian W, 2e
    via grd2_uhf_compute_data_t) -- already validated by hess_skel_open_selftest;
  * the CPHF orbital-relaxation response, built from the validated UHF CPHF
    solver (cphf_mod::cphf_solve_uhf) and the open-shell two-density
    derivative-Fock contraction (fock_deriv_mod::fock_deriv_contract_os, itself
    validated by fockx_os_selftest).

The unambiguous correctness check is the OpenQP numerical Hessian, i.e. a central
finite difference of the validated analytic UHF gradient.  For OH / STO-3G the
two must agree to the finite-difference truncation level (dx = 0.005 -> a few
1e-5), and the analytic matrix must be exactly symmetric.

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
basis=sto-3g
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
class UhfAnalyticHessian(unittest.TestCase):
    def _hessian(self, hess_type):
        from oqp.pyoqp import Runner

        workdir = Path(f"/tmp/oqp_uhf_hess_{hess_type}")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "oh.inp"
        inp.write_text(INPUT_TMPL.format(hess_type=hess_type))
        runner = Runner(project=f"oh_{hess_type}", input_file=str(inp),
                        log=str(workdir / "oh.log"))
        runner.run()
        return np.array(runner.mol.hessian)

    def test_uhf_analytic_matches_numerical(self):
        han = self._hessian("analytical")
        hnum = self._hessian("numerical")
        self.assertEqual(han.shape, hnum.shape)
        # analytic Hessian must be (essentially) exactly symmetric
        self.assertLess(np.max(np.abs(han - han.T)), 1.0e-9)
        # and match the numerical reference at the FD-truncation level
        self.assertLess(np.max(np.abs(han - hnum)), 5.0e-4)


if __name__ == "__main__":
    unittest.main()
