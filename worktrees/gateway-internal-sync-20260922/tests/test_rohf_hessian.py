"""Validate the native open-shell (ROHF) HF analytic Hessian.

The ROHF analytic Hessian (hf_hessian_mod::hf_hessian_rohf) solves the orbital
response over the docc/socc/virt rotation space with the exact ROHF orbital
Hessian:
  * operator (cphf_mod::cphf_solve_rohf): the Fock-transform part is the full
    commutator [F^s_MO, K]_vo (not the canonical Fvv K - K Foo), required because
    the raw spin-Fock vir-occ blocks are nonzero for non-canonical ROHF orbitals;
  * CPHF right-hand side: the non-canonical Pulay form with the occupied-projected
    anticommutator  sum_j ( S^x_aj F^s_ji + F^s_aj S^x_ji ),  reducing to
    eps_i S^x_ai in the canonical limit.
The orbital-relaxation response is then the central finite difference, over BOTH
geometry and the CPHF-relaxed orbital path, of the analytic open-shell electronic
gradient (nuclear repulsion added analytically).

The unambiguous correctness check is the OpenQP numerical Hessian (central finite
difference of the validated analytic ROHF gradient): the two must agree to the
finite-difference truncation level, and the analytic matrix must be exactly
symmetric.  The per-perturbation relaxed density derivatives dPa/dPb match a
re-SCF reference to ~2e-5.

Skipped unless the compiled OpenQP runtime is importable.
"""

import os
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

INPUT_TMPL = """[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=1
runtype=hess
basis=sto-3g
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
class RohfAnalyticHessian(unittest.TestCase):
    def _hessian(self, hess_type):
        from oqp.pyoqp import Runner

        workdir = Path(f"/tmp/oqp_rohf_hess_{hess_type}")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "h2op.inp"
        inp.write_text(INPUT_TMPL.format(hess_type=hess_type))
        runner = Runner(project=f"h2op_{hess_type}", input_file=str(inp),
                        log=str(workdir / "h2op.log"))
        runner.run()
        return np.array(runner.mol.hessian)

    def test_rohf_analytic_matches_numerical(self):
        han = self._hessian("analytical")
        hnum = self._hessian("numerical")
        self.assertEqual(han.shape, hnum.shape)
        # analytic ROHF Hessian must be (essentially) exactly symmetric
        self.assertLess(np.max(np.abs(han - han.T)), 1.0e-9)
        # and match the numerical reference at the FD-truncation level
        self.assertLess(np.max(np.abs(han - hnum)), 5.0e-4)


if __name__ == "__main__":
    unittest.main()
