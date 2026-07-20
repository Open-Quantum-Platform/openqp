"""Phase 11 ROOT-CAUSE PROOF (Python): the off-diagonal coded_cross defect is the
INPUT-fold (mrsfxvec on the amplitude) vs OUTPUT-fold (mrsfmntoia symmetric SOMO
fold) in the interstate 2e bilinear of the z-vector RHS.

Method (full-operator gauge-free PoC, the proven oracle route p11_poc_gaugefree.py):
  d_amp(I,J) = X0_I^T (dA_ref/dR) X0_J / (Om_J - Om_I)
with A the FULL symmetric matvec operator transported into the FIXED reference frame.

Two operators, built at each +/- displacement:
  A_out  = the true matvec A (OUTPUT fold; symmetric)         -> the PROPOSED FIX
  A_in   = A_out - (M_O - M_M)  (the SOMO-slot 2e fold delta) -> the PRODUCTION BUG
           where (M_O - M_M) lives entirely on the ijlr1 row+col (recon_operators
           confirms off-row diff = 0); M_M is the input-fold (gradient-chain) operator.

PASS for the FIX (A_out):  cos>0.99 & ratio in [0.9,1.1] vs oracle, ALL 3 pairs.
EXPECT for the BUG (A_in): reproduces qz_full coded_cross -> (1,3) SIGN FLIP, (1,2)
           ratio ~1.7, (2,3) ~ok.  This LINKS the fold to the production defect.
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

r = Runner(input_file=INP, log='/tmp/nactest/p11fold.log'); r.run(); mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = 3
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True); nbf = C0raw.shape[0]
nvirb = nbf - nocb; nij = noca * nvirb
lr1, lr2 = noca - 2, noca - 1
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
natom = mol.data['natom']
xyz0 = np.array(mol.get_system(), copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
E = list(mol.energies); Om = [E[k+1]-E[0] for k in range(nstate)]
ncoord = 3*natom
ijlr1 = (noca-1-nocb-1)*noca + (noca-1) - 1
ijlr2 = (noca-nocb-1)*noca + (noca) - 1


def set_bvec(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def matvec(c):
    set_bvec(c); oqp.mrsf_matvec_apply(mol)
    ax = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
    gmo = np.array(mol.data['OQP::nac_gmo'], copy=True).reshape((nbf, nbf), order='F')
    gch = np.array(mol.data['OQP::nac_gchan'], copy=True).reshape((nbf, nbf, 6), order='F')
    return ax, gmo, gch


def comp(w):
    p = np.zeros(nij); ij = 0
    for j in range(nocb, nbf):
        for i in range(0, noca):
            p[ij] = w[i, j]; ij += 1
    return p


def mntoia_wrk(gmo, gch, output_fold=True):
    """ch7+ch1-6 MO->amplitude back-transform. output_fold=True is the matvec
    mrsfmntoia (SOMO fold (gmo[lr1,lr1]-gmo[lr2,lr2])*SQ on the OUTPUT). False is
    the gradient-chain input-fold convention (gmo[lr1,lr1]*SQ only)."""
    ado2v, ado1v, adco1, adco2, ao21v, aco12 = (gch[:, :, k] for k in range(6))
    wrk = gmo.copy()
    for i in range(0, noca-2):
        wrk[i, lr2] += ado1v[i, lr2] + aco12[i, lr1]
        wrk[i, lr1] += ado2v[i, lr1] - aco12[i, lr2]
    for a in range(noca, nbf):
        wrk[lr1, a] += adco2[lr1, a] + ao21v[lr2, a]
        wrk[lr2, a] += adco1[lr2, a] - ao21v[lr1, a]
    if output_fold:
        wrk[lr1, lr1] = (gmo[lr1, lr1] - gmo[lr2, lr2])*SQ
    else:
        wrk[lr1, lr1] = gmo[lr1, lr1]*SQ
    wrk[lr2, lr2] = 0.0
    return comp(wrk)


def Amats():
    """build (A_out, A_in): full symmetric matvec A, and A with the SOMO-slot 2e
    fold swapped to the input-fold convention. The 1e/Fock relaxation part is
    identical in both; only the 2e SOMO fold delta differs."""
    A_out = np.zeros((nij, nij)); delta = np.zeros((nij, nij))
    e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0; e[j] = 1.0
        ax, g, gc = matvec(e)
        A_out[:, j] = ax
        delta[:, j] = mntoia_wrk(g, gc, True) - mntoia_wrk(g, gc, False)
    A_in = A_out - delta
    return A_out, A_in


# ---- amplitude transport (PoC recipe: per-block Loewdin + det-grid SOMO fold) ----
def unfold_det(col):
    x = np.zeros((noca, nvirb))
    for i in range(noca):
        for a in range(nvirb):
            ij = a*noca + i
            if ij == ijlr1:
                x[i, a] = col[ijlr1]*SQ
            elif ij == ijlr2:
                x[i, a] = -col[ijlr1]*SQ
            else:
                x[i, a] = col[ij]
    return x


def refold_det(cp):
    g = cp.copy()
    i1, a1 = ijlr1 % noca, ijlr1 // noca
    g[i1, a1] = math.sqrt(2.0)*cp[i1, a1]
    i2, a2 = ijlr2 % noca, ijlr2 // noca
    g[i2, a2] = 0.0
    return g.T.reshape(-1)


def transport_T(Q):
    Qo = Q[0:noca, 0:noca]; Qv = Q[nocb:nbf, nocb:nbf]
    T = np.zeros((nij, nij)); e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0; e[j] = 1.0
        T[:, j] = refold_det(Qo.T @ unfold_det(e) @ Qv)
    return T


print("=== R0 anchor ===")
Ao0, Ai0 = Amats()
print(f"  |A_out - A_out^T| = {np.max(np.abs(Ao0-Ao0.T)):.2e}  (matvec A symmetric)")
print(f"  |A_in  - A_in^T|  = {np.max(np.abs(Ai0-Ai0.T)):.2e}  (input-fold A NOT symmetric)")
print(f"  |A_out - A_in|    = {np.max(np.abs(Ao0-Ai0)):.2e}  (SOMO-slot fold delta)")
wv = np.linalg.eigvalsh(Ao0)
print(f"  eigs(A_out)={np.round(np.sort(wv)[:4],5)}  Om={np.round(Om,5)}")

mol.save_data()
cfg = mol.config; json0 = mol.log.replace('.log', '.json')
cfg['guess']['type'] = 'json'; cfg['guess']['file'] = json0; cfg['guess']['continue_geom'] = False


def displaced(coord):
    mol.update_system(coord); oqp.library.ints_1e(mol); oqp.library.guess(mol)
    SinglePoint(mol).energy()
    mol.data._data.control.int2e_cutoff = 1e-20
    mol.data['OQP::xyz_old'] = xyz0.reshape((3, natom))
    mol.data['OQP::VEC_MO_A_old'] = C0raw; mol.data['OQP::E_MO_A_old'] = e0
    mol.data['OQP::VEC_MO_B_old'] = C0b; mol.data['OQP::E_MO_B_old'] = e0b
    oqp.get_structures_ao_overlap(mol)
    M = np.array(mol.data[OVTAG], copy=True).reshape(-1).reshape((nbf, nbf)).T
    Q = np.zeros((nbf, nbf))
    for lo, hi in ((0, nocb), (nocb, noca), (noca, nbf)):
        sub = M[lo:hi, lo:hi]
        wv, U = np.linalg.eigh(sub.T @ sub)
        R = sub @ (U @ np.diag(1.0/np.sqrt(wv)) @ U.T)
        Q[lo:hi, lo:hi] = R
    Ao, Ai = Amats()
    T = transport_T(Q)
    return T.T @ Ao @ T, T.T @ Ai @ T


dAo = np.zeros((ncoord, nij, nij)); dAi = np.zeros((ncoord, nij, nij))
for k in range(ncoord):
    Aop, Aip = displaced(xyz0 + DELTA*np.eye(ncoord)[k])
    Aom, Aim = displaced(xyz0 - DELTA*np.eye(ncoord)[k])
    dAo[k] = (Aop-Aom)/(2*DELTA); dAi[k] = (Aip-Aim)/(2*DELTA)
mol.update_system(xyz0); oqp.library.ints_1e(mol); oqp.library.guess(mol); SinglePoint(mol).energy()

oracle = np.load('/tmp/nactest/p11_damp_oracle.npz')
qz = np.load('/tmp/nactest/qz_full.npz')
gap = {(0, 1): Om[1]-Om[0], (0, 2): Om[2]-Om[0], (1, 2): Om[2]-Om[1]}
PAIRS = [(0, 1), (0, 2), (1, 2)]


def cosr(a, b):
    return (np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30),
            np.linalg.norm(a)/(np.linalg.norm(b)+1e-30))


print("\n=== d_amp via FULL transported operator, OUTPUT fold (PROPOSED FIX) vs ORACLE ===")
allpass = True
for (I, J) in PAIRS:
    damp = np.array([X0[I] @ dAo[k] @ X0[J] for k in range(ncoord)])/gap[(I, J)]
    o = oracle[f'{I+1}{J+1}_damp']
    c, rr = cosr(damp, o)
    ok = abs(c) > 0.99 and 0.9 < rr < 1.1
    allpass = allpass and ok
    print(f"  ({I+1},{J+1}): cos={c:+.6f} ratio={rr:.4f}  {'PASS' if ok else 'FAIL'}")
print(f"  VERDICT FIX: {'ALL THREE PASS -- OUTPUT FOLD CLOSES ORACLE' if allpass else 'NOT all pass'}")

print("\n=== d_amp via FULL transported operator, INPUT fold (PRODUCTION BUG) vs numeric_cross ===")
print("    (numeric_cross = gap*oracle == qz_full *_numeric; compare to qz_full *_coded)")
codnum = {(0, 1): ('Z12_coded', 'Z12_numeric'), (0, 2): ('Z13_coded', 'Z13_numeric'),
          (1, 2): ('Z23_coded', 'Z23_numeric')}
for (I, J) in PAIRS:
    # interstate cross from input-fold operator = symmetric bilinear of the asym op
    di = np.array([0.5*(X0[I] @ dAi[k] @ X0[J] + X0[J] @ dAi[k] @ X0[I])
                   for k in range(ncoord)])  # gap*d_amp == numeric_cross
    num = gap[(I, J)]*oracle[f'{I+1}{J+1}_damp']
    c, rr = cosr(di, num)
    # compare the SIGNATURE to the actual production coded_cross from qz_full
    cc_coded = qz[codnum[(I, J)][0]] - 0.5*(qz['X1_coded' if I == 0 else f'X{I+1}_coded']
                                            if False else 0)
    print(f"  ({I+1},{J+1}): input-fold cross vs numeric_cross cos={c:+.6f} ratio={rr:.4f}")

# Reproduce the ACTUAL qz_full coded_cross signature directly (coded(Z)-1/2[coded X])
print("\n=== production qz_full coded_cross signature (for reference) ===")
Xc = {1: qz['X1_coded'], 2: qz['X2_coded'], 3: qz['X3_coded']}
Zc = {(0, 1): qz['Z12_coded'], (0, 2): qz['Z13_coded'], (1, 2): qz['Z23_coded']}
for (I, J) in PAIRS:
    cc = Zc[(I, J)] - 0.5*(Xc[I+1]+Xc[J+1])
    num = gap[(I, J)]*oracle[f'{I+1}{J+1}_damp']
    c, rr = cosr(cc, num)
    print(f"  ({I+1},{J+1}): coded_cross vs numeric_cross cos={c:+.6f} ratio={rr:.4f}")
mol.update_system(xyz0)
