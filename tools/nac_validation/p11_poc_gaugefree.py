"""Phase 11 PoC (gauge-FREE): the analytic amplitude term d_amp(I,J)=X_I^T dX_J as
   d_amp(I,J) = X0_I^T (dA_ref/dR) X0_J / (Omega_J - Omega_I),
where dA_ref is the derivative of the matvec operator A TRANSPORTED into the FIXED
reference MO frame (so the reference eigenvectors X0_I, X0_J never rotate -> no
quadratic gauge contamination).

A is SYMMETRIC in the 90-dim amplitude space (|A-A^T|=3e-13, verified), so the
perturbation-theory identity  X_I^T dX_J = X_I^T dA X_J / (Omega_J-Omega_I)  is
exact.  We build A_disp at each +/- displacement column-by-column from the matvec
(A is linear at int2e_cutoff=1e-20), transport it to the reference frame via the
amplitude-space orbital transport T (built from the per-block Procrustes Q of the
ref x displaced MO overlap), FD, and contract with the FIXED reference X0.

Oracle = benchmark semi-numerical d_amp (cos=1.0 total NAC), p11_damp_oracle.npz.
PASS: cos>0.99 and ratio in [0.9,1.1] for ALL of (1,2),(1,3),(2,3).
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

r = Runner(input_file=INP, log='/tmp/nactest/p11gf.log'); r.run(); mol = r.mol
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
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
E = list(mol.energies); Om = [E[k+1]-E[0] for k in range(nstate)]
ncoord = 3*natom
ijlr1 = (noca-1-nocb-1)*noca + (noca-1) - 1
ijlr2 = (noca-nocb-1)*noca + (noca) - 1


def set_bvec(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def Acol(c):
    set_bvec(c); oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()


def Amat():
    """full nij x nij matvec A at the current geometry/orbitals (A is linear)."""
    A = np.zeros((nij, nij)); e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0; e[j] = 1.0; A[:, j] = Acol(e)
    return A


# ---- amplitude-space orbital transport T from the per-block Procrustes Q ----
# A compressed amplitude col (nij) maps to the square M[i,a]=iatogen(col); under an
# orbital rotation Q (occ x occ, vir x vir) the amplitude transforms M -> Qo^T M Qv,
# i.e. col -> comp(Qo^T iatogen(col) Qv).  Build T (nij x nij) by columns.
def iatogen(col):
    w = np.zeros((nbf, nbf)); w[0:noca, nocb:nbf] = col.reshape(nvirb, noca).T; return w


def comp_block(w):
    return w[0:noca, nocb:nbf].T.reshape(-1)


def unfold_det(col):
    """compressed amplitude -> (noca,nvirb) DETERMINANT grid (split ijlr1 into the
    two SOMO-SOMO slots with +/-1/sqrt2; this is what transports like a grid)."""
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
    """determinant grid -> compressed amplitude (refold the singlet SOMO slot)."""
    g = cp.copy()
    i1, a1 = ijlr1 % noca, ijlr1 // noca       # SOMO1->SOMO1 grid position
    g[i1, a1] = math.sqrt(2.0)*cp[i1, a1]
    i2, a2 = ijlr2 % noca, ijlr2 // noca       # SOMO2->SOMO2 grid position
    g[i2, a2] = 0.0                            # singlet: ijlr2 slot stays empty
    out = g.T.reshape(-1)
    return out


def transport_T(Q):
    # transport in the DETERMINANT grid (the spin-adapted ijlr slot does NOT
    # transform like a plain grid entry -- unfold -> transport -> refold).
    Qo = Q[0:noca, 0:noca]
    Qv = Q[nocb:nbf, nocb:nbf]
    T = np.zeros((nij, nij)); e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0; e[j] = 1.0
        c = unfold_det(e)
        ct = Qo.T @ c @ Qv
        T[:, j] = refold_det(ct)
    return T


print("=== R0 anchor: A symmetric, eigs == Om ===")
A0 = Amat()
print(f"  |A-A^T|={np.max(np.abs(A0-A0.T)):.2e}")
wv = np.linalg.eigvalsh(A0)
print(f"  lowest eigs={np.round(np.sort(wv)[:4],5)}  Om={np.round(Om,5)}")
for k in range(nstate):
    print(f"  X{k+1}^T A X{k+1}={X0[k]@A0@X0[k]:+.6f}")

mol.save_data()
cfg = mol.config; json0 = mol.log.replace('.log', '.json')
cfg['guess']['type'] = 'json'; cfg['guess']['file'] = json0; cfg['guess']['continue_geom'] = False


def displaced_Aref(coord):
    """SCF at coord, build A_disp, transport to the reference frame: A_ref=T^T A_disp T."""
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
        # Loewdin symmetric orthogonalization: R = sub (sub^T sub)^{-1/2}.
        # This is the unique nearest-orthogonal matrix to sub with NO SVD sign
        # ambiguity, hence sign-CONTINUOUS across displacement directions
        # (the per-block Procrustes SVD sign flip was corrupting (1,2)/(1,3)).
        wv, U = np.linalg.eigh(sub.T @ sub)
        R = sub @ (U @ np.diag(1.0/np.sqrt(wv)) @ U.T)
        Q[lo:hi, lo:hi] = R
    A_disp = Amat()
    T = transport_T(Q)
    # A_disp acts on displaced-frame amps; reference amp x maps to displaced T x.
    # so A in reference frame = T^T A_disp T (gives <x_ref|A_ref|y_ref>)
    return T.T @ A_disp @ T


print("\n=== gauge-free analytic d_amp(I,J) = X0_I^T dA_ref X0_J / gap ===")
oracle = np.load('/tmp/nactest/p11_damp_oracle.npz')
dA = np.zeros((ncoord, nij, nij))
for k in range(ncoord):
    Ap = displaced_Aref(xyz0 + DELTA*np.eye(ncoord)[k])
    Am = displaced_Aref(xyz0 - DELTA*np.eye(ncoord)[k])
    dA[k] = (Ap - Am)/(2*DELTA)
mol.update_system(xyz0); oqp.library.ints_1e(mol); oqp.library.guess(mol); SinglePoint(mol).energy()

PAIRS = [(0, 1), (0, 2), (1, 2)]
allpass = True
for (I, J) in PAIRS:
    damp = np.array([X0[I] @ dA[k] @ X0[J] for k in range(ncoord)])/(Om[J]-Om[I])
    o = oracle[f'{I+1}{J+1}_damp']
    c = np.dot(damp, o)/(np.linalg.norm(damp)*np.linalg.norm(o)+1e-30)
    rr = np.linalg.norm(damp)/(np.linalg.norm(o)+1e-30)
    # global-sign fix for the residual (overall NAC sign is a convention)
    s = 1.0 if c >= 0 else -1.0
    resid = np.linalg.norm(s*damp - o); relresid = resid/(np.linalg.norm(o)+1e-30)
    ok = abs(c) > 0.99 and 0.9 < rr < 1.1
    allpass = allpass and ok
    print(f"  ({I+1},{J+1}): cos={c:+.8f} ratio={rr:.6f} resid={resid:.2e} ({100*relresid:.3f}%)  {'PASS' if ok else 'FAIL'}")
    print(f"        analytic={np.array2string(s*damp,precision=5,suppress_small=True)}")
    print(f"        oracle  ={np.array2string(o,precision=5,suppress_small=True)}")
print(f"\nVERDICT: {'ALL THREE PAIRS PASS -- ORACLE CLOSED' if allpass else 'NOT all pass'}")
mol.update_system(xyz0)
