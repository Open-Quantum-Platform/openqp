"""Numeric regression: the exact formula replica and the closed-form
gamma^formula kernel (nac-lagrangian campaign, 2026-07-31).

Asserts, on H2O/BHHLYP/6-31G* (compiled code required; skips otherwise):
  R1. The Python replica of compute_states_overlap (exact tlf=0 minors +
      8-term contraction + column normalization) reproduces the Fortran
      state overlap at a finite MO rotation to 1e-12.
  R2. The closed-form gamma^formula (cofactor sensitivities x multilinear
      contraction partials) equals a spot-check of generator-sweep
      derivatives of the replica to 1e-10.

These freeze the two central derivation results: the formula is exactly
understood, and its orbital-response kernel is available in one
linear-algebra pass.
"""
import os
import sys
import tempfile
import unittest

try:
    import oqp                      # noqa: F401 (ILP64 before numpy)
    import numpy as np
    from scipy.linalg import expm
    from oqp.pyoqp import Runner
    HAVE_OQP = True
except Exception:                   # pragma: no cover
    HAVE_OQP = False

H2O_INP = """[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g*
functional=bhhlyp
method=tdhf

[guess]
type=huckel

[scf]
multiplicity=3
type=rohf

[tdhf]
type=mrsf
nstate=3
tlf=0
"""

GATE_DIR = os.path.join(os.path.dirname(__file__), '..', 'tools',
                        'nac_lagrangian')


@unittest.skipUnless(HAVE_OQP, "compiled oqp package not importable")
class FormulaKernelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, os.path.abspath(GATE_DIR))
        cls.tmp = tempfile.TemporaryDirectory(prefix="nac_formula_test_")
        inp = os.path.join(cls.tmp.name, "h2o.inp")
        with open(inp, "w") as f:
            f.write(H2O_INP)
        runner = Runner(input_file=inp, log=inp.replace(".inp", ".log"))
        runner.run()
        cls.mol = runner.mol
        import nac_formula_kernel as K
        cls.K = K
        cls.ctx = K.build_context(cls.mol)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_replica_matches_fortran(self):
        mol, ctx, K = self.mol, self.ctx, self.K
        rng = np.random.default_rng(23)
        nbf = ctx['nbf']
        Krot = rng.standard_normal((nbf, nbf))
        Krot = Krot - Krot.T
        th = 1e-3
        SF = K.fortran_S(mol, ctx, th, Krot)
        SP = K.replica_S(ctx, expm(th * Krot))
        self.assertLessEqual(np.abs(SP - SF).max(), 1e-12)

    def test_closed_form_gamma_matches_sweep_spotcheck(self):
        ctx, K = self.ctx, self.K
        gam = K.gamma_closed(ctx)
        nbf = ctx['nbf']
        rng = np.random.default_rng(5)
        hh = 1e-4
        # spot-check 8 random generators against Richardson FD of the replica
        for _ in range(8):
            p = int(rng.integers(1, nbf))
            q = int(rng.integers(0, p))
            Kg = np.zeros((nbf, nbf))
            Kg[p, q] = 1.0
            Kg[q, p] = -1.0
            Sp1 = K.replica_S(ctx, expm(hh * Kg))
            Sm1 = K.replica_S(ctx, expm(-hh * Kg))
            Sp2 = K.replica_S(ctx, expm(2 * hh * Kg))
            Sm2 = K.replica_S(ctx, expm(-2 * hh * Kg))
            dS = (8.0 * (Sp1 - Sm1) - (Sp2 - Sm2)) / (12.0 * hh)
            pred = 2.0 * gam[:, :, p, q]      # slot convention: dS = 2*gam[p,q]
            self.assertLessEqual(np.abs(pred - dS).max(), 1e-10,
                                 msg=f"generator ({p},{q})")


if __name__ == "__main__":
    unittest.main()
