"""Phase 11 ROOT-CAUSE PROOF 2: reproduce the EXACT production coded_cross
signature from the input-fold operator using the production polarization
  coded_cross(I,J) = Q(Z_IJ) - 1/2[Q(X_I)+Q(X_J)],  Z_IJ=(X_I+X_J)/sqrt2
where Q(c) is the production single-amplitude quadratic form
  Q(c) = d/dR [ c_folded^T  Op_2e  c_input ]   (input fold on ONE leg, NO output fold)
        + the 1e/Fock relaxation part (same in both folds).

Production builds the 2e RHS as hxa = 2 G[mrsfxvec(c)]  iatogen(mrsfxvec(c))^T, i.e.
the amplitude operator acting on c is  Op_in(c) where the SOMO slot uses the INPUT
fold (mrsfxvec) on the contracted leg.  Crucially Q is a QUADRATIC FORM in c with a
NON-SYMMETRIC kernel A_in, so
  Q(c) = c^T A_in c   (as a quadratic form the antisymmetric part drops on the DIAGONAL
                       Q(X_K), but the polarization cross PICKS IT UP):
  coded_cross(I,J) = 1/2[ X_I^T A_in X_J + X_J^T A_in X_I ]  ... (symmetric, same as proof1)
        WRONG -> proof1 already showed this does NOT reproduce the (1,3) sign flip.

The real production path: Q(c) = c^T A_in c uses the SAME c on BOTH legs, BUT the two
legs fold DIFFERENTLY: leg1 = mrsfxvec(c) (folded), leg2 = the kernel G[mrsfxvec(c)]
contracted, and the PROJECTION sfrorhs adds hxa-hxa^T-hxb terms.  The net production
quadratic form is  Q(c) = c^T (A_in) c with A_in the FULL input-fold matvec-like op,
and the polarization cross is 1/2[X_I A_in X_J + X_J A_in X_I].  If proof1's symmetric
bilinear did not match, the asymmetry must enter through the gradient FD of the
NON-SYMMETRIC A_in transported frame (the transport T does not commute with the
antisymmetric part).  TEST: contract the TRANSPORTED-FRAME FD of A_in WITHOUT
re-symmetrizing -- i.e. exactly  X_I^T dA_in X_J  (the production runs the chain on Z
then subtracts, which is the un-symmetrized polarization of a non-symmetric form)."""
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

r = Runner(input_file=INP, log='/tmp/nactest/p11fold2.log'); r.run(); mol = r.mol
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
    A_out = np.zeros((nij, nij)); delta = np.zeros((nij, nij))
    e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0; e[j] = 1.0
        ax, g, gc = matvec(e)
        A_out[:, j] = ax
        delta[:, j] = mntoia_wrk(g, gc, True) - mntoia_wrk(g, gc, False)
    A_in = A_out - delta
    return A_out, A_in


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
gap = {(0, 1): Om[1]-Om[0], (0, 2): Om[2]-Om[0], (1, 2): Om[2]-Om[1]}
PAIRS = [(0, 1), (0, 2), (1, 2)]


def cosr(a, b):
    return (np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30),
            np.linalg.norm(a)/(np.linalg.norm(b)+1e-30))


def Qform(dA, c):
    """production single-amplitude quadratic form derivative c^T dA c (per coord)."""
    return np.array([c @ dA[k] @ c for k in range(ncoord)])


print("=== Reproduce production coded_cross via POLARIZATION of the single-amplitude")
print("    quadratic form Q(c)=c^T dA c (Z_IJ then subtract), each fold ===")
for tag, dA in (('OUTPUT-fold A (FIX)', dAo), ('INPUT-fold A (production BUG)', dAi)):
    print(f"\n  -- {tag} --")
    for (I, J) in PAIRS:
        Z = (X0[I]+X0[J])*SQ
        QZ = Qform(dA, Z)
        QI = Qform(dA, X0[I]); QJ = Qform(dA, X0[J])
        cc = QZ - 0.5*(QI+QJ)               # production coded_cross polarization
        num = gap[(I, J)]*oracle[f'{I+1}{J+1}_damp']  # numeric_cross
        c, rr = cosr(cc, num)
        ok = abs(c) > 0.99 and 0.9 < rr < 1.1
        print(f"     ({I+1},{J+1}): cos={c:+.6f} ratio={rr:.4f}  {'PASS' if ok else ''}")
mol.update_system(xyz0)
