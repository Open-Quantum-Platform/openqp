"""Numeric regression: the exact formula replica and the closed-form
gamma^formula kernel (nac-lagrangian campaign, 2026-07-31).

Asserts, on H2O/BHHLYP/6-31G* (compiled code required; skips otherwise):
  R1. The Python replica of compute_states_overlap (exact tlf=0 minors +
      8-term contraction + column normalization) reproduces the Fortran
      state overlap at a finite MO rotation to 1e-12.
  R2. The closed-form gamma^formula (cofactor sensitivities x multilinear
      contraction partials) equals a spot-check of generator-sweep
      derivatives of the replica to 1e-10.
  R3. The resident Fortran sparse-cofactor kernel equals gamma^formula and
      fourth-order exact-overlap generator derivatives in every orbital block.

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

        oqp.mrsf_nac_metric_data(cls.mol)
        nbf = cls.ctx['nbf']
        nstate = cls.ctx['nstate']
        flat = np.array(
            cls.mol.data['OQP::nac_gamma_tlf'], copy=True
        ).reshape(-1)
        cls.resident_gamma = flat.reshape(
            (nbf, nbf, nstate, nstate), order='F'
        ).transpose(2, 3, 0, 1)

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
        # Production TLF deliberately retains the raw projection because a
        # finite state window need not be closed under a nuclear step.  The
        # formula replica below tests the historical column-normalized
        # cofactor expression, so apply that normalization only to this oracle.
        SF /= np.linalg.norm(SF, axis=0, keepdims=True)
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

    def test_resident_fortran_gamma_matches_closed_form(self):
        reference = self.K.gamma_closed(self.ctx)
        resident_error = np.abs(self.resident_gamma - reference).max()
        if os.environ.get('NAC_GATE_REPORT'):
            print(f'resident-vs-cofactor maxdiff={resident_error:.12e}')
        self.assertLessEqual(resident_error, 1e-10)
        self.assertLessEqual(
            np.abs(
                self.resident_gamma
                + self.resident_gamma.transpose(0, 1, 3, 2)
            ).max(),
            1e-14,
        )

        # H2O freezes both facts that invalidated the retired approximation:
        # same-space response is material and state labels are not antisymmetry
        # partners after one-sided column normalization.
        nocb, noca = self.ctx['nocb'], self.ctx['noca']
        same_space = self.resident_gamma[0, 1, nocb:noca, nocb:noca]
        self.assertGreater(np.abs(same_space).max(), 0.5)
        state_antisymmetry_residual = np.abs(
            self.resident_gamma
            + self.resident_gamma.transpose(1, 0, 2, 3)
        ).max()
        self.assertGreater(state_antisymmetry_residual, 0.5)

    def test_resident_gamma_matches_all_orbital_blocks(self):
        ctx, K = self.ctx, self.K
        nbf, nocb, noca = ctx['nbf'], ctx['nocb'], ctx['noca']
        spaces = {
            'd': slice(0, nocb),
            's': slice(nocb, noca),
            'v': slice(noca, nbf),
        }
        rng = np.random.default_rng(37)
        hh = 1e-4

        for block_name in ('dd', 'ds', 'dv', 'ss', 'sv', 'vv'):
            left = spaces[block_name[0]]
            right = spaces[block_name[1]]
            generator = np.zeros((nbf, nbf))
            shape = (left.stop - left.start, right.stop - right.start)
            block = rng.standard_normal(shape)
            if block_name[0] == block_name[1]:
                generator[left, right] = block - block.T
            else:
                generator[left, right] = block
                generator[right, left] = -block.T

            plus_one = K.replica_S(ctx, expm(hh * generator))
            minus_one = K.replica_S(ctx, expm(-hh * generator))
            plus_two = K.replica_S(ctx, expm(2 * hh * generator))
            minus_two = K.replica_S(ctx, expm(-2 * hh * generator))
            derivative = (
                8.0 * (plus_one - minus_one) - (plus_two - minus_two)
            ) / (12.0 * hh)
            prediction = np.einsum(
                'ijpq,pq->ij', self.resident_gamma, generator
            )
            block_error = np.abs(prediction - derivative).max()
            if os.environ.get('NAC_GATE_REPORT'):
                print(f'block {block_name} maxdiff={block_error:.12e}')
            self.assertLessEqual(
                block_error,
                2e-9,
                msg=f"orbital block {block_name}",
            )


if __name__ == "__main__":
    unittest.main()
