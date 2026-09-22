"""Numeric regression: the derived interstate transition density gamma^IJ.

Runs the COMPILED code once (SCF + MRSF-TDDFT, H2O/BHHLYP/6-31G*), then
asserts two machine-precision identities established 2026-07-31:

  A. The closed-form gamma (dgemm over unfolded amplitudes, operator phase
     convention |i->a> = a+_{a beta} a_{i alpha}|ref>) equals the exact
     Slater-Condon 1-TDM over the 90 SF determinants, both spins, every
     state pair. This replaces the rejected sign-scanned TLF kernel.
  B. Self-consistency (code-independent): the finite-difference of the
     EXACT biorthogonal determinant overlap <Psi_I(C)|Psi_J(C e^{tK})>
     equals sum_pq gamma_pq K_pq on every ROHF rotation block, in
     particular gamma_dv = 0 exactly and the socc-socc response is REAL
     (not gauge) for S0-touching pairs.

Skipped when the compiled oqp package is not importable. Runtime ~1 min.
"""
import os
import tempfile
import unittest

try:
    import oqp                      # noqa: F401  (must precede numpy: ILP64)
    import numpy as np
    from scipy.linalg import expm
    from oqp.pyoqp import Runner
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
"""


@unittest.skipUnless(HAVE_OQP, "compiled oqp package not importable")
class GammaExactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix="nac_gamma_test_")
        inp = os.path.join(cls.tmp.name, "h2o.inp")
        with open(inp, "w") as f:
            f.write(H2O_INP)
        runner = Runner(input_file=inp, log=inp.replace(".inp", ".log"))
        runner.run()
        mol = runner.mol

        cls.nstate = mol.config["tdhf"]["nstate"]
        noca = int(np.asarray(mol.data["nelec_A"]).ravel()[0])
        nocb = noca - 2
        C = np.array(mol.data["OQP::VEC_MO_A"], copy=True)
        nbf = C.shape[0]
        nvirb = nbf - nocb
        nij = noca * nvirb
        rs = 1.0 / np.sqrt(2.0)
        X0 = np.array(mol.data["OQP::td_bvec_mo"], copy=True
                      ).reshape(-1).reshape((cls.nstate, nij)).T.copy()

        def unfold(bv, st):
            ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
            ijlr2 = (noca - nocb - 1) * noca + noca
            x = np.zeros((noca, nvirb))
            for i in range(1, noca + 1):
                for jj in range(nocb + 1, nbf + 1):
                    ij = (jj - nocb - 1) * noca + i
                    if ij == ijlr1:
                        x[i - 1, jj - nocb - 1] = bv[ijlr1 - 1, st - 1] * rs
                    elif ij == ijlr2:
                        x[i - 1, jj - nocb - 1] = -bv[ijlr1 - 1, st - 1] * rs
                    else:
                        x[i - 1, jj - nocb - 1] = bv[ij - 1, st - 1]
            return x

        cls.noca, cls.nocb, cls.nbf, cls.nvirb = noca, nocb, nbf, nvirb
        cls.Xt = [unfold(X0, s + 1) for s in range(cls.nstate)]

        ref_a = list(range(noca))
        ref_b = list(range(nocb))
        cls.dets, cls.amp_index = [], {}
        for i in range(noca):
            for a in range(nvirb):
                aocc = tuple(sorted(set(ref_a) - {i}))
                bocc = tuple(sorted(set(ref_b) | {nocb + a}))
                cls.amp_index[(i, a)] = len(cls.dets)
                cls.dets.append((aocc, bocc))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    # -- helpers ----------------------------------------------------------
    @classmethod
    def coefs(cls, x):
        c = np.zeros(len(cls.dets))
        for (i, a), idx in cls.amp_index.items():
            c[idx] = x[i, a] * ((-1.0) ** (cls.noca - 1 - i))
        return c

    @classmethod
    def tdm_exact(cls, I, J):
        nbf = cls.nbf
        ga = np.zeros((nbf, nbf))
        gb = np.zeros((nbf, nbf))
        cb, ck = cls.coefs(cls.Xt[I]), cls.coefs(cls.Xt[J])

        def one_sector(ob, ok):
            sb, sk = set(ob), set(ok)
            db = sorted(sb - sk)
            dk = sorted(sk - sb)
            if not db and not dk:
                return "same", None
            if len(db) == 1 and len(dk) == 1:
                p, q = db[0], dk[0]
                sign = ((-1) ** (len(ok) - 1 - ok.index(q))
                        * (-1) ** (len(ob) - 1 - ob.index(p)))
                return "single", (p, q, sign)
            return "far", None

        for m, (am, bm) in enumerate(cls.dets):
            if cb[m] == 0.0:
                continue
            for k, (an, bn) in enumerate(cls.dets):
                w = cb[m] * ck[k]
                if w == 0.0:
                    continue
                ka, da = one_sector(am, an)
                if ka == "far":
                    continue
                kb, db = one_sector(bm, bn)
                if kb == "far":
                    continue
                if ka == "same" and kb == "same":
                    for p in am:
                        ga[p, p] += w
                    for p in bm:
                        gb[p, p] += w
                elif ka == "single" and kb == "same":
                    p, q, s = da
                    ga[p, q] += s * w
                elif ka == "same" and kb == "single":
                    p, q, s = db
                    gb[p, q] += s * w
        return ga, gb

    @classmethod
    def gamma_closed(cls, I, J):
        """The DERIVED closed form (replaces the sign-scanned TLF kernel)."""
        nbf, noca, nocb = cls.nbf, cls.noca, cls.nocb
        xi, xj = cls.Xt[I], cls.Xt[J]
        ga = np.zeros((nbf, nbf))
        gb = np.zeros((nbf, nbf))
        ov = float(np.sum(xi * xj))
        for p in range(noca):
            ga[p, p] += ov
        for p in range(nocb):
            gb[p, p] += ov
        ga[:noca, :noca] -= xj @ xi.T
        gb[nocb:, nocb:] += xi.T @ xj
        return ga, gb

    # -- tests ------------------------------------------------------------
    def test_closed_form_equals_slater_condon(self):
        for I in range(self.nstate):
            for J in range(self.nstate):
                if I >= J:
                    continue
                ga_e, gb_e = self.tdm_exact(I, J)
                ga_c, gb_c = self.gamma_closed(I, J)
                self.assertLessEqual(np.abs(ga_e - ga_c).max(), 1e-12,
                                     msg=f"alpha gamma mismatch ({I+1},{J+1})")
                self.assertLessEqual(np.abs(gb_e - gb_c).max(), 1e-12,
                                     msg=f"beta gamma mismatch ({I+1},{J+1})")

    def test_gamma_matches_exact_overlap_response(self):
        nbf, noca, nocb = self.nbf, self.noca, self.nocb
        th = 1e-5
        rng = np.random.default_rng(11)
        cvecs = [self.coefs(self.Xt[s]) for s in range(self.nstate)]

        def exact_overlap(M):
            S = np.zeros((self.nstate, self.nstate))
            for m, (am, bm) in enumerate(self.dets):
                wa = np.array([cv[m] for cv in cvecs])
                if np.all(wa == 0.0):
                    continue
                for k, (an, bn) in enumerate(self.dets):
                    wb = np.array([cv[k] for cv in cvecs])
                    if np.all(wb == 0.0):
                        continue
                    ov = (np.linalg.det(M[np.ix_(am, an)])
                          * np.linalg.det(M[np.ix_(bm, bn)]))
                    S += np.outer(wa, wb) * ov
            return S

        gam = {}
        for I in range(self.nstate):
            for J in range(self.nstate):
                if I == J:
                    continue
                ga, gb = self.tdm_exact(min(I, J), max(I, J))
                g = ga + gb
                if I > J:
                    g = g.T
                gam[(I, J)] = g

        blocks = {"ds": (slice(0, nocb), slice(nocb, noca)),
                  "dv": (slice(0, nocb), slice(noca, nbf)),
                  "sv": (slice(nocb, noca), slice(noca, nbf)),
                  "ss": (slice(nocb, noca), slice(nocb, noca))}
        for bname, (lo, hi) in blocks.items():
            K = np.zeros((nbf, nbf))
            blk = rng.standard_normal((hi.stop - hi.start, lo.stop - lo.start))
            K[hi, lo] = blk
            K[lo, hi] = -blk.T
            dS = (exact_overlap(expm(th * K))
                  - exact_overlap(expm(-th * K))) / (2 * th)
            for I in range(self.nstate):
                for J in range(self.nstate):
                    if I >= J:
                        continue
                    an = float(np.sum(gam[(I, J)] * K))
                    self.assertLessEqual(
                        abs(dS[I, J] - an), 5e-7,
                        msg=f"block {bname} pair ({I+1},{J+1}): "
                            f"FD={dS[I, J]:+.8f} vs gamma.K={an:+.8f}")

    def test_gamma_dv_block_is_zero(self):
        nbf, noca, nocb = self.nbf, self.noca, self.nocb
        for I in range(self.nstate):
            for J in range(self.nstate):
                if I >= J:
                    continue
                ga, gb = self.tdm_exact(I, J)
                g = ga + gb
                self.assertLessEqual(np.abs(g[:nocb, noca:]).max(), 1e-14)
                self.assertLessEqual(np.abs(g[noca:, :nocb]).max(), 1e-14)


if __name__ == "__main__":
    unittest.main()
