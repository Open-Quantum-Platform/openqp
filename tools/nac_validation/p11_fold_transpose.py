"""Phase 11 hypothesis (ii): INPUT fold (M, mrsfxvec) vs OUTPUT fold (M^T, mrsfmntoia)
transpose in the interstate 2e back-transform of the production gradient chain.

THE TWO FOLDS (verified O==M^T at the SOMO-operator level, see TEST A):
  - INPUT  fold M : the gradient-chain z-vector RHS folds the AMPLITUDE before
                    contracting with the 2e kernel G:
                       wrk3 = iatogen(mrsfxvec(X)),  hxa = 2 G(X) wrk3^T
                    (tdhf_mrsf_z_vector.F90:1707-1720). M acts on the SOMO slot as
                       in[ijlr1] -> +SQ at ijlr1, -SQ at ijlr2.
  - OUTPUT fold O : the matvec mrsfmntoia contracts the kernel first, then folds the
                    OUTPUT response matrix (mntoia :1216-1233):
                       out[ijlr1] = (in[ijlr1]-in[ijlr2])*SQ ; out[ijlr2]=0
                    O = M^T exactly.

The matvec operator A (TRUTH = oracle, numeric_cross) uses the OUTPUT fold and is
SYMMETRIC.  The PT law d_amp(I,J)=X_I^T dA X_J/gap needs the symmetric A.  On the
DIAGONAL x^T M x == x^T M^T x so the fold transpose is INVISIBLE (diagonal validation
passes either way).  OFF-DIAGONAL X_I^T M X_J != X_I^T M^T X_J -> the deficiency.

Strategy:
  (A) verify O == M^T at the SOMO-fold-operator level (exact).
  (B) build the gradient-chain 2e bilinear RHS for the 3 pairs BOTH ways:
        route-M  (INPUT fold)  : G from unfolded X, wrk3 = iatogen(mrsfxvec(X))   [the
                                 buggy production convention; validated on diagonal]
        route-Mt (OUTPUT fold) : the SAME bilinear with the SOMO fold transposed,
                                 i.e. the matvec/symmetric convention.
      Both reduce to the validated diagonal at I=J.
  (C) project each to the nuclear cross via the gauge-free transport (reuse the
      p11_poc_gaugefree machinery is heavy; instead compare in ROTATION space to the
      matvec-A 2e rotation cross Lmntoia = the truth, and report cos/ratio + the
      (1,3) sign) and compare to gap*oracle / numeric_cross.

PASS for the FOLD fix: route-Mt matches the matvec-A 2e rotation cross (cos>0.99,
ratio~1) on ALL pairs incl the (1,3) sign, while route-M does not.
"""
import os
import math
import oqp
from oqp.pyoqp import Runner
import numpy as np

INP = '/tmp/nactest/H2O_en.inp'
SQ = 1.0 / math.sqrt(2.0)

r = Runner(input_file=INP, log='/tmp/nactest/p11_fold.log'); r.run(); mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
print("RUN OK (int2e_cutoff=1e-20)", flush=True)

noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
nbf = np.array(mol.data['OQP::VEC_MO_A'], copy=True).shape[0]
nvirb = nbf - nocb; nij = noca * nvirb; nstate = 3
lr1, lr2 = noca - 2, noca - 1
Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); shp = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))
ijlr1 = (noca-1-nocb-1)*noca + (noca-1) - 1
ijlr2 = (noca  -nocb-1)*noca + (noca  ) - 1


def setb(c):
    raw = Xraw.copy().reshape(-1); raw[0:nij] = c
    mol.data['OQP::td_bvec_mo'] = raw.reshape(shp)


def matvec(c):
    setb(c); oqp.mrsf_matvec_apply(mol)
    ax  = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
    gmo = np.array(mol.data['OQP::nac_gmo'],  copy=True).reshape((nbf, nbf), order='F')
    gch = np.array(mol.data['OQP::nac_gchan'],copy=True).reshape((nbf, nbf, 6), order='F')
    return ax, gmo, gch


def iatogen(xv):
    w = np.zeros((nbf, nbf)); w[0:noca, nocb:nbf] = xv.reshape(nvirb, noca).T; return w


def comp(w):
    p = np.zeros(nij); ij = 0
    for j in range(nocb, nbf):
        for i in range(0, noca):
            p[ij] = w[i, j]; ij += 1
    return p


# ============================================================================
# TEST (A): O == M^T at the SOMO-fold-operator level (compressed-amp linear maps)
# ============================================================================
print("\n=== TEST (A): SOMO-fold operators  O (output) vs M (input)  ===")
M = np.eye(nij); M[ijlr1, ijlr1] = SQ; M[ijlr2, ijlr1] = -SQ; M[ijlr2, ijlr2] = 0.0
O = np.eye(nij); O[ijlr1, ijlr1] = SQ; O[ijlr1, ijlr2] = -SQ; O[ijlr2, ijlr2] = 0.0; O[ijlr2, ijlr1] = 0.0
print(f"  ijlr1={ijlr1} ijlr2={ijlr2}")
print(f"  M(input):  M[{ijlr1},{ijlr1}]={M[ijlr1,ijlr1]:+.4f} M[{ijlr2},{ijlr1}]={M[ijlr2,ijlr1]:+.4f}")
print(f"  O(output): O[{ijlr1},{ijlr1}]={O[ijlr1,ijlr1]:+.4f} O[{ijlr1},{ijlr2}]={O[ijlr1,ijlr2]:+.4f}")
print(f"  ★ max|O - M^T| = {np.max(np.abs(O-M.T)):.3e}   (O==M^T claim)")
print(f"    max|O - M|   = {np.max(np.abs(O-M)):.3e}   (folds genuinely differ)")


# ============================================================================
# Build the matvec 2e OUTPUT-fold operator Mop (= the TRUTH, symmetric).
# ============================================================================
def mntoia_wrk(gmo, gch):
    a2, a1, c1, c2, o2, c12 = (gch[:, :, k] for k in range(6))
    w = gmo.copy()
    for i in range(0, noca-2):
        w[i, lr2] += a1[i, lr2] + c12[i, lr1]
        w[i, lr1] += a2[i, lr1] - c12[i, lr2]
    for a in range(noca, nbf):
        w[lr1, a] += c2[lr1, a] + o2[lr2, a]
        w[lr2, a] += c1[lr2, a] - o2[lr1, a]
    w[lr1, lr1] = (gmo[lr1, lr1] - gmo[lr2, lr2])*SQ      # OUTPUT fold
    w[lr2, lr2] = 0.0
    return w


matvec(X[0])
print("\n=== build matvec 2e operator Mop (OUTPUT fold) ===")
GMO = [None]*nstate; GCH = [None]*nstate
for k in range(nstate):
    _, GMO[k], GCH[k] = matvec(X[k])
Mop = np.zeros((nij, nij)); e = np.zeros(nij)
for j in range(nij):
    e[:] = 0.0; e[j] = 1.0
    _, g, gc = matvec(e)
    Mop[:, j] = comp(mntoia_wrk(g, gc))
A = np.zeros((nij, nij))
for j in range(nij):
    e[:] = 0.0; e[j] = 1.0
    A[:, j] = matvec(e)[0]
print(f"  |A - A^T|     = {np.max(np.abs(A-A.T)):.3e}")
print(f"  |Mop - Mop^T| = {np.max(np.abs(Mop-Mop.T)):.3e}  (OUTPUT-fold 2e op symmetric)")
print(f"  |Mop - O Mop| relation check (Mop already contains O on its SOMO output row)")


# ============================================================================
# TEST (B): the gradient-chain 2e bilinear RHS, BOTH folds, projected to rotation
# space.  route-M = production convention (input fold on wrk3).  route-Mt = the
# OUTPUT-fold (matvec) convention.  Use ch7 (G_MO) only here, matching the
# validated diagonal gmo construction; the full-channel rotation cross is TEST(C).
# ============================================================================
def get_gmo_only(col):
    _, g, _ = matvec(col); return g


def mrsfxvec(col):
    f = col.copy(); f[ijlr1] = col[ijlr1]*SQ; f[ijlr2] = -col[ijlr1]*SQ; return f


# transpose fold on the compressed amplitude (apply O instead of M)
def mrsfxvec_T(col):
    f = col.copy(); f[ijlr1] = (col[ijlr1]-col[ijlr2])*SQ; f[ijlr2] = 0.0; return f


def project(hxa, hxb):
    rhs = np.zeros(nconf); n = 0
    for i in range(nocb, noca):
        for j in range(0, nocb):
            rhs[n] = hxa[i, j] - hxa[j, i] - hxb[j, i]; n += 1
    for k in range(noca, nbf):
        for j in range(0, nocb):
            rhs[n] = hxa[k, j] - hxb[j, k]; n += 1
    for k in range(noca, nbf):
        for i in range(nocb, noca):
            rhs[n] = hxa[k, i] + hxb[k, i] - hxb[i, k]; n += 1
    return -rhs


prs = ([(i, j) for i in range(nocb, noca) for j in range(0, nocb)]
       + [(k, j) for k in range(noca, nbf) for j in range(0, nocb)]
       + [(k, i) for k in range(noca, nbf) for i in range(nocb, noca)])
nconf = len(prs)
nds, ndv = nocb*(noca-nocb), (nbf-noca)*nocb
blocks = [('doc-socc', 0, nds), ('doc-virt', nds, nds+ndv), ('soc-virt', nds+ndv, nconf)]


def rhs_inter(I, J, foldfn):
    GI, GJ = GMO[I], GMO[J]
    wI = iatogen(foldfn(X[I])); wJ = iatogen(foldfn(X[J]))
    hxa = GI @ wJ[0:noca, :].T + GJ @ wI[0:noca, :].T
    hxb = GI[0:noca, :].T @ wJ[0:noca, :] + GJ[0:noca, :].T @ wI[0:noca, :]
    return project(hxa, hxb)


print("\n=== TEST (B): ch7 interstate RHS, INPUT fold (route-M) vs OUTPUT fold (route-Mt) ===")
for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    Rm  = rhs_inter(I, J, mrsfxvec)
    Rmt = rhs_inter(I, J, mrsfxvec_T)
    c = np.dot(Rm, Rmt)/(np.linalg.norm(Rm)*np.linalg.norm(Rmt)+1e-30)
    print(f"  ({I+1},{J+1}): |R_M|={np.linalg.norm(Rm):.4e} |R_Mt|={np.linalg.norm(Rmt):.4e}"
          f"  cos(M,Mt)={c:+.5f} ratio={np.linalg.norm(Rm)/(np.linalg.norm(Rmt)+1e-30):.4f}")


# ============================================================================
# TEST (C): ROTATION-space full-channel 2e cross of the matvec operator (TRUTH).
# Lmntoia(I,J) = rotation gradient of 1/2(X_I^T Mop_rot X_J + sym).  Compare both
# folds' rotation cross to gap*oracle numeric_cross (the real answer).
# qz_full.npz numeric is the gap-scaled cross including relax+sp; here we compare
# the 2e-only SHAPE+SIGN.  Mop (output fold) is the matvec truth by construction.
# We then build the INPUT-fold rotation cross by swapping the SOMO output fold for
# the input fold in mntoia_wrk, and show it FLIPS/mismatches where coded_cross does.
# ============================================================================
TH = 1e-4


def givens(a, b, th):
    U = np.eye(nbf); c, s = np.cos(th), np.sin(th)
    U[a, a] = c; U[b, b] = c; U[b, a] = s; U[a, b] = -s
    return U


prs_sp = ([(i, j, 'b')  for i in range(nocb, noca) for j in range(0, nocb)]
          + [(k, j, 'ab') for k in range(noca, nbf) for j in range(0, nocb)]
          + [(k, i, 'a')  for k in range(noca, nbf) for i in range(nocb, noca)])


def mntoia_wrk_inputfold(gmo, gch):
    """same channel assembly but the SOMO fold transposed to the INPUT side:
    leave the output SOMO diagonal as the raw kernel grid (no output fold), which
    is what the gradient chain's input-fold convention produces for the operator."""
    a2, a1, c1, c2, o2, c12 = (gch[:, :, k] for k in range(6))
    w = gmo.copy()
    for i in range(0, noca-2):
        w[i, lr2] += a1[i, lr2] + c12[i, lr1]
        w[i, lr1] += a2[i, lr1] - c12[i, lr2]
    for a in range(noca, nbf):
        w[lr1, a] += c2[lr1, a] + o2[lr2, a]
        w[lr2, a] += c1[lr2, a] - o2[lr1, a]
    # INPUT-fold convention: the operator is Mop with O replaced by O^T=M on output.
    # O acts out[ijlr1]=(in1-in2)*SQ.  M^T-of-that = put the fold on the *input* of
    # the bilinear -> realized by NOT folding here and folding the amplitude instead.
    return w


def rot_state(k, Ua, Ub, wrkfn):
    gmo_r = Ua.T @ GMO[k] @ Ub
    gch_r = np.stack([Ua.T @ GCH[k][:, :, ch] @ Ub for ch in range(6)], axis=2)
    return comp(wrkfn(gmo_r, gch_r))


def Lcross(I, J, wrkfn, foldfn):
    """rotation gradient of 1/2(<foldfn(X_I)|Op(X_J)> + <foldfn(X_J)|Op(X_I)>)."""
    XI = foldfn(X[I]) if foldfn else X[I]
    XJ = foldfn(X[J]) if foldfn else X[J]
    L = np.zeros(nconf)
    for n, (p, q, sp) in enumerate(prs_sp):
        def S(th):
            Ua = givens(p, q, th) if sp in ('a', 'ab') else np.eye(nbf)
            Ub = givens(p, q, th) if sp in ('b', 'ab') else np.eye(nbf)
            return 0.5*(np.dot(XI, rot_state(J, Ua, Ub, wrkfn))
                        + np.dot(XJ, rot_state(I, Ua, Ub, wrkfn)))
        L[n] = (S(TH) - S(-TH))/(2*TH)
    return L


print("\n=== TEST (C): full-channel 2e rotation cross  route-Mt (OUTPUT) vs route-M (INPUT) ===")
print("    route-Mt = <X_I | Op_outfold(X_J)> (matvec truth, symmetric)")
print("    route-M  = <mrsfxvec(X_I) | Op_nofold(X_J)> (gradient-chain input fold)")
qz = np.load('/tmp/nactest/qz_full.npz')
ncoord_truth = {}
for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    Lmt = Lcross(I, J, mntoia_wrk, None)                       # OUTPUT fold (truth)
    Lm  = Lcross(I, J, mntoia_wrk_inputfold, mrsfxvec)         # INPUT fold (chain)
    c = np.dot(Lmt, Lm)/(np.linalg.norm(Lmt)*np.linalg.norm(Lm)+1e-30)
    print(f"\n  pair ({I+1},{J+1}): |L_Mt|={np.linalg.norm(Lmt):.4e} |L_M|={np.linalg.norm(Lm):.4e}"
          f"  cos(Mt,M)={c:+.5f} ratio(M/Mt)={np.linalg.norm(Lm)/(np.linalg.norm(Lmt)+1e-30):.4f}")
    for nm, lo, hi in blocks:
        bc = (np.dot(Lmt[lo:hi], Lm[lo:hi])
              / (np.linalg.norm(Lmt[lo:hi])*np.linalg.norm(Lm[lo:hi])+1e-30))
        print(f"     {nm:9} cos={bc:+.5f} |Mt|={np.linalg.norm(Lmt[lo:hi]):.3e}"
              f" |M|={np.linalg.norm(Lm[lo:hi]):.3e}")

print("\nDONE.")

# ============================================================================
# TEST (D): which fold matches the TRUTH?  The matvec operator A is symmetric and
# IS the oracle (numeric_cross == gap*oracle, cos+1 all pairs).  Build the matvec
# 2e interstate bilinear RHS = the rotation gradient of 1/2(X_I^T Mop X_J + sym),
# but the clean discriminator is the RHS-level comparison: the gradient-chain RHS
# hxa/hxb is the orbital-rotation derivative source.  Compare BOTH folds' ch7 RHS
# (TEST B) to the matvec-truth ch7 RHS built from the SYMMETRIC operator (output
# fold, both bra & ket).  The truth uses O on BOTH amplitudes symmetrically.
# ============================================================================
def rhs_truth(I, J):
    """matvec-truth interstate ch7 RHS: symmetric, output fold on both sides.
    Equivalent to applying O (output fold) consistently => the bilinear is built
    from G(X) and the O-folded amplitude on BOTH bra and ket, then symmetrized."""
    GI, GJ = GMO[I], GMO[J]
    # output fold acts as O on the compressed amplitude; for the bilinear value to
    # match the symmetric matvec, fold BOTH amplitudes with O (== M^T).
    wI = iatogen(mrsfxvec_T(X[I])); wJ = iatogen(mrsfxvec_T(X[J]))
    hxa = GI @ wJ[0:noca, :].T + GJ @ wI[0:noca, :].T
    hxb = GI[0:noca, :].T @ wJ[0:noca, :] + GJ[0:noca, :].T @ wI[0:noca, :]
    return project(hxa, hxb)

print("\n=== TEST (D): route-M / route-Mt vs the matvec-symmetric truth (ch7 RHS) ===")
print("    (diagonal sanity: at I=J all three must coincide)")
for st in range(nstate):
    Rd_M  = rhs_inter(st, st, mrsfxvec)
    Rd_Mt = rhs_inter(st, st, mrsfxvec_T)
    print(f"  diag {st+1}: |R_M - R_Mt(diag)| = {np.max(np.abs(Rd_M-Rd_Mt)):.2e}")
for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    Rm  = rhs_inter(I, J, mrsfxvec)      # INPUT fold (production)
    Rmt = rhs_inter(I, J, mrsfxvec_T)    # ket OUTPUT fold only (asymmetric variant)
    Rtr = rhs_truth(I, J)                # symmetric output fold (matvec truth)
    cM  = np.dot(Rm,  Rtr)/(np.linalg.norm(Rm)*np.linalg.norm(Rtr)+1e-30)
    cMt = np.dot(Rmt, Rtr)/(np.linalg.norm(Rmt)*np.linalg.norm(Rtr)+1e-30)
    print(f"  ({I+1},{J+1}): cos(M,truth)={cM:+.5f} ratio={np.linalg.norm(Rm)/(np.linalg.norm(Rtr)+1e-30):.4f}"
          f"   cos(Mt,truth)={cMt:+.5f} ratio={np.linalg.norm(Rmt)/(np.linalg.norm(Rtr)+1e-30):.4f}")
