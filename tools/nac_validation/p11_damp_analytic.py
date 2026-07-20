"""Phase 11 milestone test: replace the SEMI-NUMERICAL amplitude term
   d_amp(I,J) = X_I^T (dX_J/dR)
with the matvec-A form
   d_amp(I,J) = X_I^T (dA/dR) X_J / (Omega_J - Omega_I).
Oracle = the benchmark's semi-numerical d_amp (cos=1.0 total NAC), saved in
/tmp/nactest/p11_damp_oracle.npz.

Approach: at each +/- nuclear displacement, run SCF (relax orbitals), build the
matvec operator A from the EXPORTED per-geometry kernels (nac_gmo+nac_gchan) +
frozen-Fock fa/fb, evaluate the bilinear  B_IJ(R) = X_I(R)^T A(R) X_J(R)  with
the reference amplitudes parallel-transported (same Procrustes gauge as the
benchmark), FD it, divide by the gap.  At R0, A X_K = Omega_K X_K so
B_IJ(R0)=0 (orthogonal); the FD picks up X_I^T dA X_J + cross terms.

We test BOTH:
  (a) full bilinear FD  X_I^T dA X_J / gap            (the perturbation-theory form)
  (b) the explicit one-sided  X_I^T (A X_J)'/gap projected.
"""
import math
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint
import numpy as np

INP = '/tmp/nactest/H2O_tight_dx0.001.inp'
DELTA = 1e-3
SQ = 1.0 / math.sqrt(2.0)
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=INP, log='/tmp/nactest/p11damp.log'); r.run(); mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = 3
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True); nbf = C0raw.shape[0]
nvirb = nbf - nocb; nij = noca * nvirb
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
natom = mol.data['natom']
xyz0 = np.array(mol.get_system(), copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))   # (nstate, nij)
E = list(mol.energies); Om = [E[k+1]-E[0] for k in range(nstate)]
lr1, lr2 = noca-2, noca-1
ijlr1 = (noca-1-nocb-1)*noca + (noca-1) - 1
ijlr2 = (noca-nocb-1)*noca + (noca) - 1


def iatogen(xv):
    w = np.zeros((nbf, nbf)); w[0:noca, nocb:nbf] = xv.reshape(nvirb, noca).T; return w


def comp(w):
    p = np.zeros(nij); ij = 0
    for j in range(nocb, nbf):
        for i in range(0, noca):
            p[ij] = w[i, j]; ij += 1
    return p


def esum_full(xraw, fij, fab):
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


def mntoia_comp(gmo, gch):
    ado2v, ado1v, adco1, adco2, ao21v, aco12 = (gch[:, :, k] for k in range(6))
    wrk = gmo.copy()
    for i in range(0, noca-2):
        wrk[i, lr2] += ado1v[i, lr2] + aco12[i, lr1]
        wrk[i, lr1] += ado2v[i, lr1] - aco12[i, lr2]
    for a in range(noca, nbf):
        wrk[lr1, a] += adco2[lr1, a] + ao21v[lr2, a]
        wrk[lr2, a] += adco1[lr2, a] - ao21v[lr1, a]
    wrk[lr1, lr1] = (gmo[lr1, lr1] - gmo[lr2, lr2])*SQ
    wrk[lr2, lr2] = 0.0
    return comp(wrk)


def set_bvec(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def matvec_A(col):
    """A(col) via the exported kernels at the CURRENT geometry/orbitals."""
    set_bvec(col); oqp.mrsf_matvec_apply(mol)
    gmo = np.array(mol.data['OQP::nac_gmo'], copy=True).reshape((nbf, nbf), order='F')
    gch = np.array(mol.data['OQP::nac_gchan'], copy=True).reshape((nbf, nbf, 6), order='F')
    fa = np.array(mol.data['OQP::nac_fa'], copy=True).reshape((nbf, nbf), order='F')
    fb = np.array(mol.data['OQP::nac_fb'], copy=True).reshape((nbf, nbf), order='F')
    return mntoia_comp(gmo, gch) + esum_full(col, fa, fb)


# --- check A at R0 reproduces eigen relation: A X_K = Om_K X_K ---
print("=== R0 eigen check: |A X_K - Om_K X_K| ===")
for k in range(nstate):
    ax = matvec_A(X0[k])
    print(f"  state {k+1}: |A X - Om X|={np.linalg.norm(ax-Om[k]*X0[k]):.2e}  Om={Om[k]:.5f} <X|AX>={np.dot(X0[k],ax):.5f}")

# --- displaced amplitude transport (same gauge as benchmark) ---
mol.save_data()
cfg = mol.config; json0 = mol.log.replace('.log', '.json')
cfg['guess']['type'] = 'json'; cfg['guess']['file'] = json0; cfg['guess']['continue_geom'] = False


def unfold(bv_st):
    x = np.zeros((noca, nvirb))
    for i in range(1, noca+1):
        for jj in range(nocb+1, nbf+1):
            ij = (jj-nocb-1)*noca + i
            if ij == ijlr1+1:
                x[i-1, jj-nocb-1] = bv_st[ijlr1]*SQ
            elif ij == ijlr2+1:
                x[i-1, jj-nocb-1] = -bv_st[ijlr1]*SQ
            else:
                x[i-1, jj-nocb-1] = bv_st[ij-1]
    return x


def refold(cp):
    f = cp.copy(); f[nocb, 0] = math.sqrt(2.0)*cp[nocb, 0]; f[nocb+1, 1] = 0.0
    return f.reshape(-1, order='F')[: nij] if False else f.T.reshape(-1)


def displaced(coord):
    """relax SCF at coord; return (transported reference amps in displaced basis,
    the displaced-basis matvec evaluator state, Om at displaced)."""
    mol.update_system(coord); oqp.library.ints_1e(mol); oqp.library.guess(mol)
    SinglePoint(mol).energy()
    mol.data['OQP::xyz_old'] = xyz0.reshape((3, natom))
    mol.data['OQP::VEC_MO_A_old'] = C0raw; mol.data['OQP::E_MO_A_old'] = e0
    mol.data['OQP::VEC_MO_B_old'] = C0b; mol.data['OQP::E_MO_B_old'] = e0b
    oqp.get_structures_ao_overlap(mol)
    M = np.array(mol.data[OVTAG], copy=True).reshape(-1).reshape((nbf, nbf)).T
    Q = np.zeros((nbf, nbf))
    for lo, hi in ((0, nocb), (nocb, noca), (noca, nbf)):
        sub = M[lo:hi, lo:hi]
        W, s, Vt = np.linalg.svd(sub)
        Q[lo:hi, lo:hi] = Vt.T @ W.T          # nearest-orthogonal Procrustes rotation
    # transport the REFERENCE amplitudes X0 into the displaced MO basis
    Xt = np.zeros((nstate, nij))
    for st in range(nstate):
        c = unfold(X0[st]); cp = Q[:noca, :noca].T @ c @ Q[nocb:, nocb:]
        Xt[st] = refold(cp)
        if np.dot(X0[st], Xt[st]) < 0:       # global phase alignment to reference
            Xt[st] *= -1.0
    return Xt


def Bmat(Xset, col_I, col_J):
    """X_I^T A X_J at the current geometry, with transported amps."""
    return np.dot(col_I, matvec_A(col_J))


print("\n=== analytic d_amp(I,J) = X_I^T dA X_J / gap  vs semi-numeric oracle ===")
oracle = np.load('/tmp/nactest/p11_damp_oracle.npz')
PAIRS = [(0, 1), (0, 2), (1, 2)]
ncoord = 3*natom


def AX_all(Xt):
    """matvec A applied to each transported amplitude at the current geometry."""
    return [matvec_A(Xt[st]) for st in range(nstate)]


# B[k][sign] = nstate x nstate bilinear  Xt_I^T A Xt_J  at displacement k
Bp_all = np.zeros((ncoord, nstate, nstate)); Bm_all = np.zeros((ncoord, nstate, nstate))
for k in range(ncoord):
    Xt = displaced(xyz0 + DELTA*np.eye(ncoord)[k]); ax = AX_all(Xt)
    for I in range(nstate):
        for J in range(nstate):
            Bp_all[k, I, J] = np.dot(Xt[I], ax[J])
    Xt = displaced(xyz0 - DELTA*np.eye(ncoord)[k]); ax = AX_all(Xt)
    for I in range(nstate):
        for J in range(nstate):
            Bm_all[k, I, J] = np.dot(Xt[I], ax[J])

for (I, J) in PAIRS:
    damp = ((Bp_all[:, I, J] - Bm_all[:, I, J])/(2*DELTA))/(Om[J]-Om[I])
    o = oracle[f'{I+1}{J+1}_damp']
    c = np.dot(damp, o)/(np.linalg.norm(damp)*np.linalg.norm(o)+1e-30)
    rr = np.linalg.norm(damp)/(np.linalg.norm(o)+1e-30)
    print(f"  ({I+1},{J+1}): cos={c:+.5f} ratio={rr:.4f} |analytic|={np.linalg.norm(damp):.5f} |oracle|={np.linalg.norm(o):.5f}")
    print(f"        analytic={np.array2string(damp,precision=4,suppress_small=True)}")
    print(f"        oracle  ={np.array2string(o,precision=4,suppress_small=True)}")
mol.update_system(xyz0)
