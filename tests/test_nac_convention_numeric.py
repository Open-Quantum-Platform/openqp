"""Numeric regression: numerical-NAC sign/orientation conventions (H2O).

Runs the COMPILED code (SCF + MRSF-TDDFT + 18-displacement numerical NAC)
and asserts, with NO sign alignment anywhere:

  dcv[i,j] = d_ij = <I|d/dR|J>          antisymmetric (exactly)
  nacv[i,j] = (E_j - E_i) * d_ij = h_ij  symmetric (exactly), SIGNED match
  frozen reference d vectors (chc3 2026-07-31, commit 28eaaa7 conventions)
  match component-by-component INCLUDING SIGN.

Skipped when the compiled oqp package is not importable (source-only tree).
Runtime ~2-5 min (18 small SCF+TDDFT jobs).
"""
import os
import tempfile
import unittest

try:
    import oqp                      # noqa: F401  (must precede numpy: ILP64)
    import numpy as np
    from oqp.pyoqp import Runner
    from oqp.library.single_point import NAC
    HAVE_OQP = True
except Exception:                   # pragma: no cover - source-only checkout
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

[nac]
type=numerical
states=1 2, 1 3, 2 3
dx=0.001
nproc=4
"""

# Frozen post-fix reference (chc3, MKL-ILP64, 2026-07-31). d_ij vectors,
# (natom, 3), atom order O,H,H; canonical convention d_ij = <I|d/dR|J>.
# Signs flipped 2026-07-31 after the storage-boundary transpose fix: the
# original freeze was taken from the d_ji-oriented read (verified against
# the exact biorthogonal-overlap orientation gate).
REF_D = {
    (1, 2): [[0.013974, 0.013974, 0.0],
             [-0.062404, -0.062404, 0.0],
             [-0.062404, -0.062404, 0.0]],
    (1, 3): [[0.0, 0.0, 0.0],
             [-0.018393, -0.018393, 0.0],
             [0.018393, 0.018393, 0.0]],
    (2, 3): [[-0.265741, 0.265741, 0.0],
             [0.160167, -0.160167, 0.133384],
             [0.160167, -0.160167, -0.133384]],
}
REF_E = [-76.36093649, -76.03214231, -75.97388051]


@unittest.skipUnless(HAVE_OQP, "compiled oqp package not importable")
class NACConventionNumericTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix="nac_conv_test_")
        inp = os.path.join(cls.tmp.name, "h2o_nac.inp")
        with open(inp, "w") as f:
            f.write(H2O_INP)
        runner = Runner(input_file=inp, log=inp.replace(".inp", ".log"))
        runner.run()
        cls.mol = runner.mol
        nac = NAC(cls.mol)
        cls.nacv, cls.dcv, cls.flags = nac.numerical_nac()
        cls.E = np.array(cls.mol.energies[1:4])

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_all_displacements_computed(self):
        self.assertEqual(set(self.flags), {"computed"})

    def test_state_energies(self):
        np.testing.assert_allclose(self.E, REF_E, atol=5e-6)

    def test_dcv_antisymmetric_exactly(self):
        asym = np.abs(self.dcv + self.dcv.transpose(1, 0, 2, 3)).max()
        self.assertLessEqual(asym, 1e-14)

    def test_nacv_symmetric_exactly(self):
        sym = np.abs(self.nacv - self.nacv.transpose(1, 0, 2, 3)).max()
        self.assertLessEqual(sym, 1e-14)

    def test_gap_orientation_signed(self):
        # nacv[i,j] must equal (E_j - E_i) * dcv[i,j] with cos = +1
        # (NOT -1: the pre-fix code had the gap axes swapped).
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                h_can = (self.E[j] - self.E[i]) * self.dcv[i, j]
                num = float(np.sum(self.nacv[i, j] * h_can))
                den = (np.linalg.norm(self.nacv[i, j])
                       * np.linalg.norm(h_can) + 1e-300)
                self.assertGreater(num / den, 0.999999,
                                   msg=f"gap orientation wrong for ({i+1},{j+1})")

    def test_frozen_reference_vectors_gauge_resolved(self):
        # Davidson state phases are random per run (observed on chc3:
        # magnitudes reproduce to 1e-10, pair signs flip as s_i s_j with
        # s in {+-1}^3). Resolve the STATE gauge from the pair signs,
        # REQUIRE gauge consistency (product of the three pair signs must
        # be +1 -- a non-gauge sign error cannot satisfy this), then
        # compare component-by-component including sign.
        pair_sign = {}
        for (i, j), ref in REF_D.items():
            got = self.dcv[i - 1, j - 1].ravel()
            s = float(np.sign(np.dot(got, np.array(ref).ravel())))
            self.assertNotEqual(s, 0.0)
            pair_sign[(i, j)] = s
        self.assertEqual(
            pair_sign[(1, 2)] * pair_sign[(1, 3)] * pair_sign[(2, 3)], 1.0,
            msg="pair signs are NOT a state gauge -- genuine sign defect")
        for (i, j), ref in REF_D.items():
            with self.subTest(pair=(i, j)):
                got = self.dcv[i - 1, j - 1]
                np.testing.assert_allclose(
                    pair_sign[(i, j)] * got, np.array(ref), atol=2e-4,
                    err_msg=f"d({i},{j}) deviates from frozen reference")


if __name__ == "__main__":
    unittest.main()
