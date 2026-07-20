"""HELD-OUT STATE-4 TEST for the Phase-11 gauge-free closure.

Extends BHHLYP/H2O to FOUR MRSF states and tests the NEW pairs involving state 4:
(1,4),(2,4),(3,4) -- pairs the implementation never saw.

ORACLE (LHS):  d_amp_orc(I,J) = X0_I^T dX_J  built EXACTLY like benchmark_full_nac.py
   (the producer of the validated p11_damp_oracle.npz): SVD-Procrustes per-block Q,
   transport the DISPLACED eigenvector amplitudes into the reference frame in the
   determinant grid (unfold/transport/refold), per-state sign-align, compressed dot.
   Validated below: the three OLD pairs reproduce the saved oracle to cos=+/-1.0.

ANALYTIC (RHS): d_amp_ana(I,J) = X0_I^T dA_ref X0_J / gap, the p11_poc_gaugefree route,
   with the Loewdin det-grid transport UNCHANGED from the 3-state PoC.

If the closure is GENUINE the two must agree for the held-out state-4 pairs to the FD
floor.  We ALSO re-verify the three old pairs in the 4-state run (cross-check that
nothing broke, and that the saved oracle is reproduced).
"""
import math
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint
import numpy as np

INP = '/tmp/nactest/H2O_tight_dx0.001_4st.inp'
DELTA = 1e-3
SQ = 1.0 / math.sqrt(2.0)
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=INP, log='/tmp/nactest/p11gf4.log'); r.run(); mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = 4
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True); nbf = C0raw.shape[0]
nvirb = nbf - nocb; nij = noca * nvirb
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
natom = mol.data['natom']
xyz0 = np.array(mol.get_system(), copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))            # compressed reference amps
E = list(mol.energies); Om = [E[k + 1] - E[0] for k in range(nstate)]
ncoord = 3 * natom
ijlr1 = (noca - 1 - nocb - 1) * noca + (noca - 1) - 1
ijlr2 = (noca - nocb - 1) * noca + (noca) - 1


def set_bvec(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def Acol(c):
    set_bvec(c); oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()


def Amat():
    A = np.zeros((nij, nij)); e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0; e[j] = 1.0; A[:, j] = Acol(e)
    return A


def unfold_det(col):
    x = np.zeros((noca, nvirb))
    for i in range(noca):
        for a in range(nvirb):
            ij = a * noca + i
            if ij == ijlr1:
                x[i, a] = col[ijlr1] * SQ
            elif ij == ijlr2:
                x[i, a] = -col[ijlr1] * SQ
            else:
                x[i, a] = col[ij]
    return x


def refold_det(cp):
    g = cp.copy()
    i1, a1 = ijlr1 % noca, ijlr1 // noca
    g[i1, a1] = math.sqrt(2.0) * cp[i1, a1]
    i2, a2 = ijlr2 % noca, ijlr2 // noca
    g[i2, a2] = 0.0
    return g.T.reshape(-1)


def transport_T(Q):
    """amplitude-space transport from per-block orbital rotation Q (p11 PoC, UNCHANGED)."""
    Qo = Q[0:noca, 0:noca]; Qv = Q[nocb:nbf, nocb:nbf]
    T = np.zeros((nij, nij)); e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0; e[j] = 1.0
        c = unfold_det(e); ct = Qo.T @ c @ Qv; T[:, j] = refold_det(ct)
    return T


def loewdin_Q(M):
    Q = np.zeros((nbf, nbf))
    for lo, hi in ((0, nocb), (nocb, noca), (noca, nbf)):
        sub = M[lo:hi, lo:hi]
        wv, U = np.linalg.eigh(sub.T @ sub)
        Q[lo:hi, lo:hi] = sub @ (U @ np.diag(1.0 / np.sqrt(wv)) @ U.T)
    return Q


def svd_Q(M):
    Q = np.zeros((nbf, nbf))
    for lo, hi in ((0, nocb), (nocb, noca), (noca, nbf)):
        W, sv, Vt = np.linalg.svd(M[lo:hi, lo:hi]); Q[lo:hi, lo:hi] = Vt.T @ W.T
    return Q


print("=== R0 anchor (4-state) ===")
A0 = Amat()
print(f"  |A-A^T|={np.max(np.abs(A0 - A0.T)):.2e}")
print(f"  lowest eigs={np.round(np.sort(np.linalg.eigvalsh(A0))[:6], 5)}")
print(f"  Om        ={np.round(Om, 5)}")
for k in range(nstate):
    print(f"  X{k + 1}^T A X{k + 1}={X0[k] @ A0 @ X0[k]:+.6f}")

mol.save_data()
cfg = mol.config; json0 = mol.log.replace('.log', '.json')
cfg['guess']['type'] = 'json'; cfg['guess']['file'] = json0; cfg['guess']['continue_geom'] = False


def get_M(coord):
    mol.update_system(coord); oqp.library.ints_1e(mol); oqp.library.guess(mol)
    SinglePoint(mol).energy()
    mol.data['OQP::xyz_old'] = xyz0.reshape((3, natom))
    mol.data['OQP::VEC_MO_A_old'] = C0raw; mol.data['OQP::E_MO_A_old'] = e0
    mol.data['OQP::VEC_MO_B_old'] = C0b; mol.data['OQP::E_MO_B_old'] = e0b
    oqp.get_structures_ao_overlap(mol)
    return np.array(mol.data[OVTAG], copy=True).reshape(-1).reshape((nbf, nbf)).T


def displaced_Aref(coord):
    """RHS analytic side: A_ref = T^T A_disp T with the LOEWDIN det-grid transport (p11)."""
    M = get_M(coord)
    T = transport_T(loewdin_Q(M))
    A_disp = Amat()
    return T.T @ A_disp @ T


def displaced_amps_oracle(coord):
    """LHS oracle side: transport the DISPLACED eigenvector amplitudes into the ref
    frame with the SVD-Procrustes det-grid transport (EXACT benchmark_full_nac route)."""
    M = get_M(coord)
    Q = svd_Q(M)
    Xd = np.array(mol.data['OQP::td_bvec_mo'], copy=True).reshape(-1).reshape((nstate, nij))
    Qo = Q[:noca, :noca]; Qv = Q[nocb:, nocb:]
    out = np.zeros((nstate, nij))
    for st in range(nstate):
        c = unfold_det(Xd[st]); cp = Qo.T @ c @ Qv
        out[st] = refold_det(cp)
        if np.dot(X0[st], out[st]) < 0:
            out[st] *= -1.0
    return out


print("\n=== building FD of A_ref (analytic RHS) and of transported amps (oracle LHS) ===")
dA = np.zeros((ncoord, nij, nij))
dX = np.zeros((ncoord, nstate, nij))
for k in range(ncoord):
    dA[k] = (displaced_Aref(xyz0 + DELTA * np.eye(ncoord)[k]) -
             displaced_Aref(xyz0 - DELTA * np.eye(ncoord)[k])) / (2 * DELTA)
    Xp = displaced_amps_oracle(xyz0 + DELTA * np.eye(ncoord)[k])
    Xm = displaced_amps_oracle(xyz0 - DELTA * np.eye(ncoord)[k])
    dX[k] = (Xp - Xm) / (2 * DELTA)
mol.update_system(xyz0); oqp.library.ints_1e(mol); oqp.library.guess(mol); SinglePoint(mol).energy()


def report(I, J, saved=None):
    rhs = np.array([X0[I] @ dA[k] @ X0[J] for k in range(ncoord)]) / (Om[J] - Om[I])
    lhs = np.array([X0[I] @ dX[k, J] for k in range(ncoord)])      # = X0_I^T dX_J oracle
    c = np.dot(rhs, lhs) / (np.linalg.norm(rhs) * np.linalg.norm(lhs) + 1e-30)
    rr = np.linalg.norm(rhs) / (np.linalg.norm(lhs) + 1e-30)
    s = 1.0 if c >= 0 else -1.0
    resid = np.linalg.norm(s * rhs - lhs); rel = resid / (np.linalg.norm(lhs) + 1e-30)
    ok = abs(c) > 0.99 and 0.9 < rr < 1.1
    line = (f"  ({I + 1},{J + 1}): cos={c:+.8f} ratio={rr:.6f} resid={resid:.2e} "
            f"({100 * rel:.3f}%) |oracle|={np.linalg.norm(lhs):.5f}  {'PASS' if ok else 'FAIL'}")
    print(line)
    print(f"        analytic={np.array2string(s * rhs, precision=5, suppress_small=True)}")
    print(f"        oracle  ={np.array2string(lhs, precision=5, suppress_small=True)}")
    oks = ok
    if saved is not None:
        cs = np.dot(lhs, saved) / (np.linalg.norm(lhs) * np.linalg.norm(saved) + 1e-30)
        rs = np.linalg.norm(lhs) / (np.linalg.norm(saved) + 1e-30)
        print(f"        [oracle reproduces SAVED 3-state oracle: cos={cs:+.6f} ratio={rs:.4f}]")
        oks = oks and abs(cs) > 0.99 and 0.9 < rs < 1.1
    return oks, c, rr, resid, np.linalg.norm(lhs)


print("\n=== HELD-OUT pairs involving state 4 (never seen) ===")
results = {}
allheld = True
for (I, J) in [(0, 3), (1, 3), (2, 3)]:
    ok, c, rr, rs, nl = report(I, J)
    results[(I + 1, J + 1)] = (c, rr, rs, nl); allheld = allheld and ok

print("\n=== OLD pairs (cross-check + reproduce saved oracle) ===")
oracle3 = np.load('/tmp/nactest/p11_damp_oracle.npz')
allold = True
for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    ok, c, rr, rs, nl = report(I, J, saved=oracle3[f'{I + 1}{J + 1}_damp'])
    results[(I + 1, J + 1)] = (c, rr, rs, nl); allold = allold and ok

print(f"\nHELD-OUT (1,4),(2,4),(3,4): {'ALL PASS' if allheld else 'SOME FAIL'}")
print(f"OLD pairs (reproduce oracle):{'ALL PASS' if allold else 'SOME FAIL'}")
print("RESULTS:", {k: (round(v[0], 7), round(v[1], 6), float(f'{v[2]:.2e}'), round(v[3], 5)) for k, v in results.items()})
mol.update_system(xyz0)
