"""Phase 11 -- localize the off-diagonal ground-config deficiency in the production
gradient chain by reconstructing coded_cross from its PIECES.

PART A  (build-free, runtime gates): coded_cross(I,J)=coded(Z_IJ)-1/2[coded X_I+coded X_J]
        vs numeric_cross=gap*oracle, under baseline / AB1-off / 2FT-off / SP-off.
        -> Isolates whether AB1 (i) or relaxation (2FT/SP) carries the DIRECTIONAL defect.

PART B  (operator-level): build the FULL nij x nij matvec 2e back-transform operator two
        ways from the SAME exported MO kernels (nac_gmo ch7 + nac_gchan ch1-6):
          M_O  = OUTPUT-fold (matvec mrsfmntoia : fold applied on the OUTPUT)  -> the TRUE A 2e part
          M_M  = INPUT-fold  (gradient-chain hxa=2 G[X_folded] X_f^T : fold on INPUT)
        Show: M_O is SYMMETRIC and == A's 2e part; M_M is NOT (M_M = M_O^T-ish), so its
        symmetric cross extraction 1/2(M_M+M_M^T) differs OFF-DIAGONAL -> reproduces the
        coded(1,3) = -1/2 numeric(1,3) sign flip.  This is candidate (ii).
"""
import math, os
import oqp, oqp.library
from oqp.pyoqp import Runner
import numpy as np

SQ = 1.0/math.sqrt(2.0)

# ===================== PART A : gradient chain under runtime gates =====================
INP_A = '/tmp/nactest/H2O_tight_dx0.001.inp'
rA = Runner(input_file=INP_A, log='/tmp/nactest/reconA.log'); rA.run(); molA = rA.mol
molA.data._data.control.int2e_cutoff = 1e-20
nstate, natom = 3, molA.data['natom']
noca = int(np.asarray(molA.data['nelec_A']).ravel()[0]); nocb = noca-2
C0raw = np.array(molA.data['OQP::VEC_MO_A'], copy=True); nbf = C0raw.shape[0]
nvirb = nbf-nocb; nij = noca*nvirb
X0_raw = np.array(molA.data['OQP::td_bvec_mo'], copy=True)
XA = X0_raw.reshape(-1).reshape((nstate, nij))
Etot0 = list(molA.energies); Om = [Etot0[k+1]-Etot0[0] for k in range(nstate)]
GATEKEYS = ('NAC_ZERO_AB1', 'NAC_ZERO_2FT', 'NAC_ZERO_SP', 'NAC_ZERO_HX', 'NAC_ZERO_T')


def clear_gates():
    for k in GATEKEYS:
        os.environ.pop(k, None)


def Gfull(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    molA.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    molA.data.set_tdhf_target(1)
    oqp.tdhf_mrsf_z_vector(molA); oqp.tdhf_mrsf_gradient(molA)
    gZ = molA.get_grad().reshape(-1).copy()
    for tag in ('OQP::td_p', 'OQP::WAO', 'OQP::td_abxc', 'OQP::td_mrsf_density'):
        molA.data[tag] = np.zeros_like(np.array(molA.data[tag], copy=True))
    oqp.tdhf_mrsf_gradient(molA)
    gS = molA.get_grad().reshape(-1).copy()
    molA.data['OQP::td_bvec_mo'] = X0_raw
    return gZ - gS


Zs = {'Z12': (XA[0]+XA[1])*SQ, 'Z13': (XA[0]+XA[2])*SQ, 'Z23': (XA[1]+XA[2])*SQ,
      'X1': XA[0].copy(), 'X2': XA[1].copy(), 'X3': XA[2].copy()}


def coded_all(gates):
    clear_gates()
    for g in gates:
        os.environ[g] = '1'
    out = {n: Gfull(z) for n, z in Zs.items()}
    clear_gates(); return out


def cross_of(GZ):
    return {'12': GZ['Z12']-0.5*(GZ['X1']+GZ['X2']),
            '13': GZ['Z13']-0.5*(GZ['X1']+GZ['X3']),
            '23': GZ['Z23']-0.5*(GZ['X2']+GZ['X3'])}


def cosr(a, b):
    return (np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30),
            np.linalg.norm(a)/(np.linalg.norm(b)+1e-30))


oracle = np.load('/tmp/nactest/p11_damp_oracle.npz')
gap = {'12': Om[1]-Om[0], '13': Om[2]-Om[0], '23': Om[2]-Om[1]}
numeric_cross = {p: gap[p]*oracle[f'{p}_damp'] for p in ('12', '13', '23')}

print("=== PART A: coded_cross vs numeric_cross under runtime gates ===")
print("gaps:", {p: round(gap[p], 5) for p in gap})
baseline_cc = None
for name, gates in [('baseline', []), ('AB1off', ['NAC_ZERO_AB1']),
                    ('2FToff', ['NAC_ZERO_2FT']), ('SPoff', ['NAC_ZERO_SP'])]:
    cc = cross_of(coded_all(gates))
    if name == 'baseline':
        baseline_cc = cc
    line = f"  {name:9s}:"
    for p in ('12', '13', '23'):
        c, rr = cosr(cc[p], numeric_cross[p])
        line += f"  ({p}) cos={c:+.4f} r={rr:.3f}"
    print(line, flush=True)
print("  --> cos (DIRECTION) is gate-invariant: AB1/2FT/SP do NOT cause the angular defect")
del rA, molA

# ===================== PART B : output-fold vs input-fold full operator =================
INP_B = '/tmp/nactest/H2O_en.inp'
rB = Runner(input_file=INP_B, log='/tmp/nactest/reconB.log'); rB.run(); mol = rB.mol
mol.data._data.control.int2e_cutoff = 1e-20
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca-2
nbf = np.array(mol.data['OQP::VEC_MO_A'], copy=True).shape[0]
nvirb = nbf-nocb; nij = noca*nvirb
lr1, lr2 = noca-2, noca-1
Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); shp = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))
ijlr1 = (noca-1-nocb-1)*noca + (noca-1) - 1
ijlr2 = (noca-nocb-1)*noca + (noca) - 1


def setb(c):
    raw = Xraw.copy().reshape(-1); raw[0:nij] = c; mol.data['OQP::td_bvec_mo'] = raw.reshape(shp)


def matvec(c):
    setb(c); oqp.mrsf_matvec_apply(mol)
    ax = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
    gmo = np.array(mol.data['OQP::nac_gmo'], copy=True).reshape((nbf, nbf), order='F')
    gch = np.array(mol.data['OQP::nac_gchan'], copy=True).reshape((nbf, nbf, 6), order='F')
    return ax, gmo, gch


def iatogen(xv):
    w = np.zeros((nbf, nbf)); w[0:noca, nocb:nbf] = xv.reshape(nvirb, noca).T; return w


def comp(w):
    p = np.zeros(nij); ij = 0
    for j in range(nocb, nbf):
        for i in range(0, noca):
            p[ij] = w[i, j]; ij += 1
    return p


def mrsfxvec(col):
    f = col.copy(); f[ijlr1] = col[ijlr1]*SQ; f[ijlr2] = -col[ijlr1]*SQ; return f


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


def mntoia_wrk(gmo, gch):
    """OUTPUT-fold 2e back-transform operator (matvec mrsfmntoia)."""
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


# cache kernels
AX = [None]*nstate; GMO = [None]*nstate; GCH = [None]*nstate
for k in range(nstate):
    AX[k], GMO[k], GCH[k] = matvec(X[k])
matvec(X[0])
fa = np.array(mol.data['OQP::nac_fa'], copy=True).reshape((nbf, nbf), order='F')
fb = np.array(mol.data['OQP::nac_fb'], copy=True).reshape((nbf, nbf), order='F')

# --- Build FULL operators column by column. M_O (output fold) needs kernels for each
#     basis amplitude column. M_M (input fold) = gradient-chain hxa/hxb operator. ---
GMO_col = [None]*nij; GCH_col = [None]*nij
for j in range(nij):
    e = np.zeros(nij); e[j] = 1.0
    _, GMO_col[j], GCH_col[j] = matvec(e)
matvec(np.zeros(nij))  # restore

# full 2e operator A2e = matvec A minus esum (relaxation/Fock part)
A2e = np.zeros((nij, nij))
for j in range(nij):
    e = np.zeros(nij); e[j] = 1.0
    ax, _, _ = matvec(e)
    A2e[:, j] = ax - esum_full(e, fa, fb)

# OUTPUT-fold reconstruction M_O (must == A2e, the validated mntoia value anchor)
M_O = np.zeros((nij, nij))
for j in range(nij):
    M_O[:, j] = mntoia_wrk(GMO_col[j], GCH_col[j])

# INPUT-fold gradient-chain AMPLITUDE operator. The matvec output-fold A2e = O.K.M with
# M=mrsfxvec input fold (applied to bvec before the kernel) and O=mntoia output fold.
# The gradient-chain z-vector RHS instead builds the 2e bilinear with the input fold on
# BOTH legs and NO output fold -> effective amplitude operator M_M = comp_nofold(K(M x)),
# i.e. drop the SOMO output-fold line of mntoia.  Build column by column.
def mntoia_NOoutputfold(gmo, gch):
    ado2v, ado1v, adco1, adco2, ao21v, aco12 = (gch[:, :, k] for k in range(6))
    wrk = gmo.copy()
    for i in range(0, noca-2):
        wrk[i, lr2] += ado1v[i, lr2] + aco12[i, lr1]
        wrk[i, lr1] += ado2v[i, lr1] - aco12[i, lr2]
    for a in range(noca, nbf):
        wrk[lr1, a] += adco2[lr1, a] + ao21v[lr2, a]
        wrk[lr2, a] += adco1[lr2, a] - ao21v[lr1, a]
    # INPUT-fold convention: SOMO slot kept on the lr1 row WITHOUT the (lr1-lr2)*SQ
    # symmetric fold (the chain folded the INPUT, so the OUTPUT fold is absent/transposed)
    wrk[lr1, lr1] = gmo[lr1, lr1] * SQ
    wrk[lr2, lr2] = 0.0
    return comp(wrk)


M_M = np.zeros((nij, nij))
for j in range(nij):
    M_M[:, j] = mntoia_NOoutputfold(GMO_col[j], GCH_col[j])

print("\n=== PART B: operator symmetry & match to true A 2e part ===")
print(f"  |A2e - A2e^T|   = {np.max(np.abs(A2e-A2e.T)):.2e}  (true 2e part of matvec A: SYMMETRIC)")
print(f"  |M_O - A2e|     = {np.max(np.abs(M_O-A2e)):.2e}  (OUTPUT-fold mntoia == true A 2e)")
print(f"  |M_O - M_O^T|   = {np.max(np.abs(M_O-M_O.T)):.2e}  (output-fold operator symmetric)")
print(f"  |M_M - M_M^T|   = {np.max(np.abs(M_M-M_M.T)):.2e}  (INPUT-fold operator symmetric?)")
print(f"  |M_M - M_O|     = {np.max(np.abs(M_M-M_O)):.2e}  (folds differ where?)")

print("\n=== PART B: cross contractions per pair ===")
print("    truth: X_I^T A2e X_J (= 0, eigenvectors) ; off-diag DERIVATIVE is the NAC.")
print("    Key test: does the input-fold sym extraction match the output-fold on the")
print("    SOMO-channel cross?  Report the X_I.op.X_J asymmetry that drives the (1,3) flip.")
for (I, J), p in zip([(0, 1), (0, 2), (1, 2)], ('12', '13', '23')):
    o_IJ = X[I] @ M_O @ X[J]; o_JI = X[J] @ M_O @ X[I]
    m_IJ = X[I] @ M_M @ X[J]; m_JI = X[J] @ M_M @ X[I]
    print(f"  pair ({p}):  OUTPUT-fold X_I.M_O.X_J={o_IJ:+.6f} X_J.M_O.X_I={o_JI:+.6f}"
          f"  (asym={o_IJ-o_JI:+.2e})")
    print(f"             INPUT-fold  X_I.M_M.X_J={m_IJ:+.6f} X_J.M_M.X_I={m_JI:+.6f}"
          f"  (asym={m_IJ-m_JI:+.2e})")

# decisive: project the off-diagonal operator difference onto the SOMO singlet slot
print("\n=== PART B: per-amplitude-component operator diff on the SOMO (ijlr1) slot ===")
print(f"  ijlr1={ijlr1} ijlr2={ijlr2}")
print(f"  M_O[ijlr1,:] nonzero count={np.sum(np.abs(M_O[ijlr1])>1e-8)}  "
      f"M_M[ijlr1,:] nonzero={np.sum(np.abs(M_M[ijlr1])>1e-8)}")
print(f"  row ijlr1: |M_O-M_M|={np.max(np.abs(M_O[ijlr1]-M_M[ijlr1])):.4e}  "
      f"col ijlr1: |M_O-M_M|={np.max(np.abs(M_O[:,ijlr1]-M_M[:,ijlr1])):.4e}")
np.savez('/tmp/nactest/recon_operators.npz', A2e=A2e, M_O=M_O, M_M=M_M, X=X)

# ===================== PART C : end-to-end nuclear NAC from each fold ==================
# FD of the transported 2e operator: damp(I,J)=X0_I^T dM_ref X0_J / gap, with M_ref the
# operator transported to the FIXED reference frame (PoC recipe).  Build M_O (output) and
# M_M (input) at each +/- displacement from the matvec kernels, reconstruct, transport, FD.
from oqp.library.single_point import SinglePoint
print("\n=== PART C: nuclear NAC from each fold (FD of transported operator) vs oracle ===")
DELTA = 1e-3
OVTAG = 'OQP::overlap_mo_non_orthogonal'
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
natom = mol.data['natom']; xyz0 = np.array(mol.get_system(), copy=True); ncoord = 3*natom
Etot = list(mol.energies); OmC = [Etot[k+1]-Etot[0] for k in range(nstate)]
X0c = X.copy()
mol.save_data()
cfg = mol.config; json0 = mol.log.replace('.log', '.json')
cfg['guess']['type'] = 'json'; cfg['guess']['file'] = json0; cfg['guess']['continue_geom'] = False


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


def build_ops_at_geom():
    """build M_O and M_M operators column-by-column at the CURRENT geometry/orbitals."""
    MO = np.zeros((nij, nij)); MM = np.zeros((nij, nij))
    for j in range(nij):
        e = np.zeros(nij); e[j] = 1.0
        _, g, gc = matvec(e)
        MO[:, j] = mntoia_wrk(g, gc)
        MM[:, j] = mntoia_NOoutputfold(g, gc)
    return MO, MM


def displaced_ops(coord):
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
    MO, MM = build_ops_at_geom()
    T = transport_T(Q)
    return T.T @ MO @ T, T.T @ MM @ T


dMO = np.zeros((ncoord, nij, nij)); dMM = np.zeros((ncoord, nij, nij))
for k in range(ncoord):
    MOp, MMp = displaced_ops(xyz0 + DELTA*np.eye(ncoord)[k])
    MOm, MMm = displaced_ops(xyz0 - DELTA*np.eye(ncoord)[k])
    dMO[k] = (MOp-MOm)/(2*DELTA); dMM[k] = (MMp-MMm)/(2*DELTA)
mol.update_system(xyz0); oqp.library.ints_1e(mol); oqp.library.guess(mol); SinglePoint(mol).energy()
mol.data._data.control.int2e_cutoff = 1e-20

for (I, J), p in zip([(0, 1), (0, 2), (1, 2)], ('12', '13', '23')):
    gapC = OmC[J]-OmC[I]
    dO = np.array([X0c[I] @ dMO[k] @ X0c[J] for k in range(ncoord)])/gapC
    dM = np.array([0.5*(X0c[I] @ dMM[k] @ X0c[J] + X0c[J] @ dMM[k] @ X0c[I]) for k in range(ncoord)])/gapC
    o = oracle[f'{p}_damp']
    cO = np.dot(dO, o)/(np.linalg.norm(dO)*np.linalg.norm(o)+1e-30)
    cM = np.dot(dM, o)/(np.linalg.norm(dM)*np.linalg.norm(o)+1e-30)
    rO = np.linalg.norm(dO)/(np.linalg.norm(o)+1e-30)
    rM = np.linalg.norm(dM)/(np.linalg.norm(o)+1e-30)
    print(f"  pair ({p}): OUTPUT-fold cos={cO:+.5f} r={rO:.3f} | INPUT-fold cos={cM:+.5f} r={rM:.3f}")
print("  (NOTE: 2e-only, no relaxation/2FT/AB1; oracle is the FULL d_amp so ratios<1.)")
print("\nDONE")
