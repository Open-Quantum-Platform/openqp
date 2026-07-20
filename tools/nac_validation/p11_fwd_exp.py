"""Phase 11 FORWARD-part experiment (read-only).

mntoia_interstate.py reconstructs the matvec 2e back-transform L^mntoia from the
exported MO kernels (gmo, gchan) but rotates ONLY the already-built kernels as
Ua^T K Ub. That is the BACK-transform; it reproduces only ~half of the
gradient-chain diagonal RHS (cos 0.77-0.9998, ratio 0.24-0.48) because it MISSES
the FORWARD part: the C-dependence of the channel CONSTRUCTION inside mrsfcbc
(the AO kernel K^AO(C) itself changes when the orbitals rotate, because the
channels are assembled from the occupied/SOMO orbitals).

THIS SCRIPT obtains the forward part by FINITE-DIFFERENCING the exported channel
kernels under the SAME orbital rotation theta:
    rotate VEC_MO_A/B by U(theta) = C.U  (Givens on the rotating pair),
    recompute gmo,gchan via mrsf_matvec_apply (tight int2e_cutoff=1e-20),
    read them back. The exported kernel is then  K_exp(th) = (CU)^T K^AO(CU) (CU).
The pure BACK-transform alone is  K_back(th) = U(th)^T [C^T K^AO(C) C] U(th)
                                            = U(th)^T K_exp(0) U(th).
So the FORWARD increment of the *AO* kernel, expressed in the rotated MO basis, is
    K_fwd(th) = K_exp(th) - K_back(th)  =  (CU)^T [K^AO(CU) - K^AO(C)] (CU).
The full mntoia output gradient is built by feeding the FULL recomputed kernels
K_exp(th) into mntoia_wrk and FD'ing S(th); the back-only gradient feeds
K_back(th). Their difference is L^fwd.

Compared to:
 (1) DIAGONAL: full L^mntoia (= back + forward) vs the gradient-chain 2e RHS dump
     (-(G_MO+mrsfsp), gates SP off / 2FT,AB1 zeroed). Does ratio -> 1, cos -> 1?
 (2) OFF-DIAGONAL (rotation space): residual R(I,J) = L^mntoia_full - L^chain_2e.
     Reported per block + interstate. (The deficiency D in qz_nogate.npz is in
     NUCLEAR space; a literal cos(R,D) needs the z-vector->nuclear projection, so
     here R is reported in rotation space and only the magnitude/structure tested.)

================================================================================
EMPIRICAL CONCLUSION (this run, tight int2e_cutoff=1e-20):
================================================================================
The naive orbital-rotation-FD forward part does NOT close the diagonal. It makes
it WORSE, and the forward increment is ~orthogonal to the chain RHS:

  state 1: BACK cos=+0.834 ratio=1.29 ->  FULL cos=+0.127 ratio=7.39
  state 2: BACK cos=+0.9998 ratio=0.242 -> FULL cos=+0.046 ratio=1.71
  state 3: BACK cos=+0.956 ratio=0.478 -> FULL cos=+0.063 ratio=8.08
  forward increment |L_fwd|>>|chain|, cos(L_fwd, chain) = -0.019/-0.096/+0.006
  doc-virt ratioFULL blows up to 25-76x.

ROOT CAUSE (probed directly): the forward increment d(gmo)/dtheta is
  - 71% concentrated in the SOMO rows/cols (lr1=4, lr2=5), max at (4,18);
  - collapses 11x when the SOMO-row/col amplitudes of X are zeroed.
This is the SAME OUTPUT-SIDE SOMO-FOLD ARTIFACT documented in
PHASE11_fd_findings.md (the xlr*d(fa/fb)/dtheta background). Rotating the real
orbitals and recomputing the channels via mrsf_matvec_apply re-applies the
matvec's output-side sqrt2 SOMO fold keyed to the rotated orbitals, against the
frozen huge ground-pair amplitude xlr~=0.998 -> spurious background that swamps
the genuine forward channel-construction physics.

OFF-DIAGONAL: the full L^mntoia is dominated by the same artifact; cos(L_full,
L_chain7) ~ 0.0007/0.0001/0.288 -> the residual R is NOT the deficiency (it is the
artifact). (1,3) doc-virt |R|=0.349 is entirely the artifact, not the sign-flip
physics.

VERDICT: the forward part of the matvec 2e rotation gradient CANNOT be obtained by
finite-differencing the recomputed channels under a real orbital rotation, for the
same structural reason the raw FD-of-matvec failed: the SOMO fold is on the OUTPUT.
The correct forward part must come from the ANALYTIC channel C-derivative WITHOUT
re-applying the output fold (i.e. fold on the INPUT, as the validated diagonal
hxa=2*G_MO*X_folded already does). The factor-2 in that validated diagonal is
forward+back with the fold held on the input; the back-only reconstruction here
(M_back) misses exactly that input-folded forward half -- which is why M_back
ratio ~ 0.24-0.48 but the FIX is NOT to add the output-folded forward FD.
"""
import oqp
from oqp.pyoqp import Runner
import numpy as np
import os

INP = '/tmp/nactest/H2O_en.inp'
SQ = 1.0 / np.sqrt(2.0)

r = Runner(input_file=INP, log='/tmp/nactest/p11_fwd.log'); r.run(); mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20    # TIGHT: matvec linear to 1e-15
print("RUN OK (tight int2e_cutoff=1e-20)", flush=True)

noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
C0a = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
nbf = C0a.shape[0]
nvirb = nbf - nocb; nij = noca * nvirb; nstate = 3
lr1, lr2 = noca - 2, noca - 1
Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); shp = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))


def setb(c):
    raw = Xraw.copy().reshape(-1); raw[0:nij] = c
    mol.data['OQP::td_bvec_mo'] = raw.reshape(shp)


def set_mo(Ca, Cb):
    mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(Ca)
    mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(Cb)


def matvec_kernels(c, Ca=None, Cb=None):
    """run mrsf_matvec_apply for amplitude c with (optionally rotated) MOs Ca,Cb,
    return (gmo, gchan) as exported = mo^T K^AO mo in the *given* MO basis."""
    if Ca is not None:
        set_mo(Ca, Cb)
    setb(c)
    oqp.mrsf_matvec_apply(mol)
    gmo = np.array(mol.data['OQP::nac_gmo'], copy=True).reshape((nbf, nbf), order='F')
    gch = np.array(mol.data['OQP::nac_gchan'], copy=True).reshape((nbf, nbf, 6), order='F')
    if Ca is not None:
        set_mo(C0a, C0b)        # restore
    return gmo, gch


# unrotated kernels + frozen Fock
GMO = [None]*nstate; GCH = [None]*nstate
for k in range(nstate):
    GMO[k], GCH[k] = matvec_kernels(X[k])
matvec_kernels(X[0])
fa = np.array(mol.data['OQP::nac_fa'], copy=True).reshape((nbf, nbf), order='F')
fb = np.array(mol.data['OQP::nac_fb'], copy=True).reshape((nbf, nbf), order='F')


def iatogen(xv):
    w = np.zeros((nbf, nbf)); w[0:noca, nocb:nbf] = xv.reshape(nvirb, noca).T; return w


def comp(w):
    p = np.zeros(nij); ij = 0
    for j in range(nocb, nbf):
        for i in range(0, noca):
            p[ij] = w[i, j]; ij += 1
    return p


def esum_full(xraw, fij, fab):           # verbatim mrsfesum (mult=1)
    wrk = iatogen(xraw); scr = wrk.copy(); scr[lr1, lr1] = 0; scr[lr2, lr2] = 0
    t = np.zeros((nbf, nbf)); t[0:noca, nocb:nbf] = scr[0:noca, nocb:nbf]@fab[nocb:nbf, nocb:nbf]
    t[0:noca, nocb:nbf] += -fij[0:noca, 0:noca]@scr[0:noca, nocb:nbf]
    xlr = wrk[lr1, lr1]; w1 = np.zeros((nbf, nbf))
    for j in range(nocb, nbf):
        for i in range(0, noca):
            w1[i, j] += t[i, j]
            if i == lr1: w1[i, j] += fab[j, lr1]*xlr*SQ
            if i == lr2: w1[i, j] -= fab[j, lr2]*xlr*SQ
            if j == lr1: w1[i, j] -= fij[i, lr1]*xlr*SQ
            if j == lr2: w1[i, j] += fij[i, lr2]*xlr*SQ
    dumn = (-fij[lr1, 0:noca]@scr[0:noca, lr1]+fij[lr2, 0:noca]@scr[0:noca, lr2]
            + fab[lr1, nocb:nbf]@scr[lr1, nocb:nbf]-fab[lr2, nocb:nbf]@scr[lr2, nocb:nbf])
    w1[lr1, lr1] = dumn*SQ+xlr*(fab[lr1, lr1]+fab[lr2, lr2]-fij[lr1, lr1]-fij[lr2, lr2])*0.5
    w1[lr2, lr2] = 0
    return comp(w1)


def mntoia_wrk(gmo, gch):
    """reconstruct the mrsfmntoia nbf x nbf wrk (before projection)."""
    ado2v, ado1v, adco1, adco2, ao21v, aco12 = (gch[:, :, k] for k in range(6))
    wrk = gmo.copy()
    for i in range(0, noca-2):                       # doc
        wrk[i, lr2] += ado1v[i, lr2] + aco12[i, lr1]   # sec3
        wrk[i, lr1] += ado2v[i, lr1] - aco12[i, lr2]   # sec4
    for a in range(noca, nbf):                       # virt beta
        wrk[lr1, a] += adco2[lr1, a] + ao21v[lr2, a]   # sec5
        wrk[lr2, a] += adco1[lr2, a] - ao21v[lr1, a]   # sec6
    wrk[lr1, lr1] = (gmo[lr1, lr1] - gmo[lr2, lr2])*SQ  # output fold
    wrk[lr2, lr2] = 0.0
    return wrk


def mntoia_comp(gmo, gch):
    return comp(mntoia_wrk(gmo, gch))


# ---- VALUE anchor (sanity, same as mntoia_interstate STEP 2) ----
print("\n=== VALUE anchor: reconstructed mntoia(X) vs (nac_mvax - esum) ===")
for st in range(1, nstate+1):
    ax, gmo, gch = (lambda c: (setb(c), oqp.mrsf_matvec_apply(mol),
                   np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()))(X[st-1])[2], None, None
    gmo, gch = GMO[st-1], GCH[st-1]
    M_fortran = ax - esum_full(X[st-1], fa, fb)
    print(f"  state {st}: max|M_python - M_fortran| = "
          f"{np.max(np.abs(mntoia_comp(gmo, gch)-M_fortran)):.2e}")

# ============ ROTATION GRADIENT ============
prs = ([(i, j, 'b') for i in range(nocb, noca) for j in range(0, nocb)]
       + [(k, j, 'ab') for k in range(noca, nbf) for j in range(0, nocb)]
       + [(k, i, 'a') for k in range(noca, nbf) for i in range(nocb, noca)])
nconf = len(prs)
nds, ndv = nocb*(noca-nocb), (nbf-noca)*nocb
blocks = [('doc-socc', 0, nds), ('doc-virt', nds, nds+ndv), ('soc-virt', nds+ndv, nconf)]
TH = 1e-4


def givens(a, b, th):
    U = np.eye(nbf); c, s = np.cos(th), np.sin(th)
    U[a, a] = c; U[b, b] = c; U[b, a] = s; U[a, b] = -s     # C'=C U : col a = c a + s b
    return U


# ---------- BACK-only (= mntoia_interstate M_rot): rotate fixed kernels ----------
def M_back(k, Ua, Ub):
    gmo_r = Ua.T @ GMO[k] @ Ub
    gch_r = np.stack([Ua.T @ GCH[k][:, :, ch] @ Ub for ch in range(6)], axis=2)
    return mntoia_comp(gmo_r, gch_r)


# ---------- FULL (back + forward): recompute kernels at rotated orbitals ----------
def M_full(k, Ua, Ub):
    """rotate the ORBITALS, recompute the AO kernels via the matvec, read exported
    kernels = (C Ua)^T K^AO(C Ua) (C Ua)  [rows alpha, cols beta]. This is the FULL
    theta-dependence: forward (channel construction from rotated C) + back."""
    Ca = C0a @ Ua; Cb = C0b @ Ub
    gmo_r, gch_r = matvec_kernels(X[k], Ca, Cb)
    return mntoia_comp(gmo_r, gch_r)


def L_of(Mfun, I, J):
    """rotation gradient of S_IJ = 1/2(X_I . M(X_J) + X_J . M(X_I))."""
    L = np.zeros(nconf)
    for n, (p, q, sp) in enumerate(prs):
        def S(th):
            Ua = givens(p, q, th) if sp in ('a', 'ab') else np.eye(nbf)
            Ub = givens(p, q, th) if sp in ('b', 'ab') else np.eye(nbf)
            return 0.5*(np.dot(X[I], Mfun(J, Ua, Ub)) + np.dot(X[J], Mfun(I, Ua, Ub)))
        L[n] = (S(TH) - S(-TH))/(2*TH)
    return L


def dump_chain_2e(target):
    """gradient-chain 2e RHS (G_MO + mrsfsp): SP off, 2FT+AB1 zeroed = -(L^chain_2e)."""
    mol.data.set_tdhf_target(target)
    for k in ('NAC_ZERO_2FT', 'NAC_ZERO_AB1', 'NAC_ZERO_HX', 'NAC_ZERO_T', 'NAC_ZERO_SP'):
        os.environ.pop(k, None)
    os.environ['NAC_ZERO_2FT'] = '1'; os.environ['NAC_ZERO_AB1'] = '1'
    os.environ['NAC_DUMP_RHS'] = '1'
    oqp.tdhf_mrsf_z_vector(mol)
    for k in ('NAC_ZERO_2FT', 'NAC_ZERO_AB1', 'NAC_DUMP_RHS'):
        os.environ.pop(k, None)
    return np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).ravel()[:nconf]


def report(tag, L, Rc):
    nL, nR = np.linalg.norm(L), np.linalg.norm(Rc)
    c = np.dot(L, Rc)/(nL*nR+1e-30)
    print(f"  {tag}: |L|={nL:.4e} |Rc|={nR:.4e} cos={c:+.5f} ratio={nL/(nR+1e-30):.4f}")
    return c


# ---- DIAGONAL: BACK vs FULL vs chain ----
print("\n=== DIAGONAL: L^mntoia(I,I) vs chain 2e RHS  [-(G_MO+mrsfsp)] ===")
print("  (chain RHS sign is -L^chain_2e; compare directly, cos>0 expected)")
Rdiag = {}
for st in range(1, nstate+1):
    Rc = dump_chain_2e(st); Rdiag[st] = Rc
    Lb = L_of(M_back, st-1, st-1)
    Lf = L_of(M_full, st-1, st-1)
    print(f"\n state {st}:")
    report("BACK only ", Lb, Rc)
    report("FULL(b+f) ", Lf, Rc)
    # forward increment
    Lfwd = Lf - Lb
    print(f"    forward increment |L_fwd|={np.linalg.norm(Lfwd):.4e} "
          f"cos(L_fwd,Rc)={np.dot(Lfwd,Rc)/(np.linalg.norm(Lfwd)*np.linalg.norm(Rc)+1e-30):+.5f}")
    # per-block ratio of FULL
    for nm, lo, hi in blocks:
        nl, nr = np.linalg.norm(Lf[lo:hi]), np.linalg.norm(Rc[lo:hi])
        cc = np.dot(Lf[lo:hi], Rc[lo:hi])/(nl*nr+1e-30)
        print(f"      blk {nm:9s} cosFULL={cc:+.4f} ratioFULL={nl/(nr+1e-30):.4f}")

# ---- OFF-DIAGONAL: full L^mntoia and residual vs chain (rotation space) ----
print("\n=== OFF-DIAGONAL (rotation space): residual R = L^mntoia_full - L^chain_2e ===")
print("  L^chain_2e(I,J) = bilinear of the diagonal dumps via parallelogram:")
print("  build chain cross by polarization of the diagonal RHS is NOT valid (RHS is")
print("  linear in target only); instead use the gmo_interstate bilinear as L^chain.")

# chain interstate 2e RHS = the gmo_interstate bilinear (G_MO channel-7 ONLY)
ijlr1 = (noca-1-nocb-1)*noca + (noca-1) - 1
ijlr2 = (noca-nocb-1)*noca + (noca) - 1


def mrsfxvec(col):
    f = col.copy(); f[ijlr1] = col[ijlr1]*SQ; f[ijlr2] = -col[ijlr1]*SQ
    return f


def iatogen_f(colf):
    av = np.zeros((nbf, nbf)); av[0:noca, nocb:nbf] = colf.reshape((nvirb, noca)).T
    return av


def project_pm(hxa, hxb):
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


wf = [iatogen_f(mrsfxvec(X[k])) for k in range(nstate)]


def chain_2e_inter(I, J):
    """channel-7 bilinear 2e RHS (gmo_interstate.py); reduces to diag at I=J."""
    GI, GJ, wI, wJ = GMO[I], GMO[J], wf[I], wf[J]
    hxa = GI @ wJ[0:noca, :].T + GJ @ wI[0:noca, :].T
    hxb = GI[0:noca, :].T @ wJ[0:noca, :] + GJ[0:noca, :].T @ wI[0:noca, :]
    return project_pm(hxa, hxb)


for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    Lf = L_of(M_full, I, J)
    Lb = L_of(M_back, I, J)
    Lch7 = chain_2e_inter(I, J)
    R = Lf - Lch7
    nLf, nLch = np.linalg.norm(Lf), np.linalg.norm(Lch7)
    cfc = np.dot(Lf, Lch7)/(nLf*nLch+1e-30)
    print(f"\n pair ({I+1},{J+1}): |L_full|={nLf:.4e} |L_back|={np.linalg.norm(Lb):.4e} "
          f"|L_chain7|={nLch:.4e}")
    print(f"   cos(L_full, L_chain7)={cfc:+.5f}  |R=L_full-L_chain7|={np.linalg.norm(R):.4e}")
    for nm, lo, hi in blocks:
        print(f"     blk {nm:9s} |R|={np.linalg.norm(R[lo:hi]):.3e} "
              f"|L_full|={np.linalg.norm(Lf[lo:hi]):.3e}")

print("\nDONE.")
