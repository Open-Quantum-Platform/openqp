"""Phase 11 Route A: reconstruct the matvec 2e back-transform mrsfmntoia from the
exported channel MO kernels (OQP::nac_gchan + nac_gmo), validate the VALUE against
the Fortran (nac_mvax - esum), then build its interstate ROTATION GRADIENT L^mntoia
and compare to the gradient chain [G_MO + mrsfsp] -> residual R = the candidate
off-diagonal fix (O==M^T output-fold transpose hypothesis).

mrsfmntoia (tdhf_mrsf_lib.F90:996-1260), mult=1, mapped to MO kernels (ROHF mo_a==mo_b):
  base   wrk = gmo (= mo_a^T agdlr mo_a, channel 7)
  sec3   wrk[i,lr2] += ado1v[i,lr2] + aco12[i,lr1]   (i in doc=0..noca-3)
  sec4   wrk[i,lr1] += ado2v[i,lr1] - aco12[i,lr2]   (i in doc)
  sec5   wrk[lr1,a] += adco2[lr1,a] + ao21v[lr2,a]   (a in virt_b=noca..nbf-1)
  sec6   wrk[lr2,a] += adco1[lr2,a] - ao21v[lr1,a]   (a in virt_b)
  fold   wrk[lr1,lr1] = (gmo[lr1,lr1]-gmo[lr2,lr2])*sq ; wrk[lr2,lr2] = 0
  proj   M[i,j] for i in 0..noca-1, j in nocb..nbf-1
channel index (gchan[:,:,k], k=0..5) = fmrst ch1..6 = ado2v,ado1v,adco1,adco2,ao21v,aco12.
"""
import oqp
from oqp.pyoqp import Runner
import numpy as np

INP = '/tmp/nactest/H2O_en.inp'
SQ = 1.0/np.sqrt(2.0)

r = Runner(input_file=INP, log='/tmp/nactest/mntoia.log'); r.run(); mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20    # TIGHT screening: the default 5e-11
# ball-based density screening corrupts the small mixed channels 5,6 (o21v/co12),
# making the matvec artifactually nonlinear. 1e-20 -> matvec linear to 1e-15.
print("RUN OK (tight int2e_cutoff=1e-20)", flush=True)
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca-2
nbf = np.array(mol.data['OQP::VEC_MO_A'], copy=True).shape[0]
nvirb = nbf-nocb; nij = noca*nvirb; nstate = 3
lr1, lr2 = noca-2, noca-1
Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); shp = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))


def setb(c):
    raw = Xraw.copy().reshape(-1); raw[0:nij] = c
    mol.data['OQP::td_bvec_mo'] = raw.reshape(shp)


def matvec(c):
    setb(c); oqp.mrsf_matvec_apply(mol)
    ax = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
    gmo = np.array(mol.data['OQP::nac_gmo'], copy=True).reshape((nbf, nbf), order='F')
    gch = np.array(mol.data['OQP::nac_gchan'], copy=True).reshape((nbf, nbf, 6), order='F')
    return ax, gmo, gch


matvec(X[0])      # populate nac_fa/fb tags
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
    """reconstruct the mrsfmntoia nbf x nbf wrk (before projection) from MO kernels."""
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


# cache per-state kernels
AX = [None]*nstate; GMO = [None]*nstate; GCH = [None]*nstate
for k in range(nstate):
    AX[k], GMO[k], GCH[k] = matvec(X[k])


def M_of(k):                       # reconstructed pure-2e operator output for state k
    return mntoia_comp(GMO[k], GCH[k])


# ---- STEP 2: VALUE anchor: reconstructed mntoia == nac_mvax - esum ----
print("\n=== STEP 2: reconstructed mntoia(X) vs Fortran (nac_mvax - esum) ===")
for st in range(1, nstate+1):
    M_fortran = AX[st-1] - esum_full(X[st-1], fa, fb)
    print(f"  state {st}: max|M_python - M_fortran| = {np.max(np.abs(M_of(st-1)-M_fortran)):.2e}"
          f"   |M|={np.linalg.norm(M_fortran):.4e}")

# ---- STEP 1b: operator M linearity (must hold: M = A - esum, both linear) ----
print("\n=== STEP 1b: M linearity  M(aX1+bX3) vs aM(X1)+bM(X3) ===")
a, b = 0.7, -1.3
axm = matvec(a*X[0]+b*X[2])[0]                       # A.(aX1+bX3)
M_comb = axm - esum_full(a*X[0]+b*X[2], fa, fb)
print(f"  max|M(aX1+bX3) - (aM1+bM3)| = {np.max(np.abs(M_comb-(a*M_of(0)+b*M_of(2)))):.2e}")

# ============ ROTATION GRADIENT of the mntoia 2e bilinear ============
# Spin-resolved ROHF rotation (sfrogen): doc-socc=beta, doc-virt=both, soc-virt=alpha.
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


def M_rot(k, Ua, Ub):
    """reconstructed mntoia output for state k with rotated orbitals (channels held
    fixed in AO => MO kernels rotate as Ua^T K Ub; rows=alpha, cols=beta)."""
    gmo_r = Ua.T @ GMO[k] @ Ub
    gch_r = np.stack([Ua.T @ GCH[k][:, :, ch] @ Ub for ch in range(6)], axis=2)
    return mntoia_comp(gmo_r, gch_r)


def Lmntoia(I, J):
    """rotation gradient of S_IJ = 1/2(X_I . M(X_J) + X_J . M(X_I))."""
    L = np.zeros(nconf)
    for n, (p, q, sp) in enumerate(prs):
        def S(th):
            Ua = givens(p, q, th) if sp in ('a', 'ab') else np.eye(nbf)
            Ub = givens(p, q, th) if sp in ('b', 'ab') else np.eye(nbf)
            return 0.5*(np.dot(X[I], M_rot(J, Ua, Ub)) + np.dot(X[J], M_rot(I, Ua, Ub)))
        L[n] = (S(TH) - S(-TH))/(2*TH)
    return L


import os


def dump_chain_2e(target):
    """gradient-chain 2e RHS (G_MO + mrsfsp), no 2FT, no AB1: gates SP OFF (keep
    mrsfsp), 2FT+AB1 ON->zeroed. = -(L^chain_2e)."""
    mol.data.set_tdhf_target(target)
    for k in ('NAC_ZERO_2FT', 'NAC_ZERO_AB1', 'NAC_ZERO_HX', 'NAC_ZERO_T', 'NAC_ZERO_SP'):
        os.environ.pop(k, None)
    os.environ['NAC_ZERO_2FT'] = '1'; os.environ['NAC_ZERO_AB1'] = '1'
    os.environ['NAC_DUMP_RHS'] = '1'
    oqp.tdhf_mrsf_z_vector(mol)
    for k in ('NAC_ZERO_2FT', 'NAC_ZERO_AB1', 'NAC_DUMP_RHS'):
        os.environ.pop(k, None)
    return np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).ravel()[:nconf]


# ---- STEP 3: DIAGONAL anchor: Lmntoia(I,I) vs chain 2e RHS (resolve sign/factor) ----
print("\n=== STEP 3: diagonal Lmntoia(I,I) vs chain 2e RHS dump (-(G_MO+mrsfsp)) ===")
for st in range(1, nstate+1):
    L = Lmntoia(st-1, st-1)
    Rc = dump_chain_2e(st)
    for sgn in (+1, -1):
        c = np.dot(L, sgn*Rc)/(np.linalg.norm(L)*np.linalg.norm(Rc)+1e-30)
        if c > 0.5:
            print(f"  state {st}: |L|={np.linalg.norm(L):.4e} |Rc|={np.linalg.norm(Rc):.4e} "
                  f"cos={c:+.5f} ratio(L/Rc)={np.linalg.norm(L)/(np.linalg.norm(Rc)+1e-30):.4f} sign={sgn:+d}")
