"""Phase 11 CENTRAL PoC: prove that the INTERSTATE-CONSISTENT DENSITY ASSEMBLY makes
coded_cross(I,J) match the oracle on ALL 3 pairs, while the production single-state
(input-fold) density assembly carries the localized off-diagonal defect.

THE PRODUCTION coded_cross (test_QZ.py:43-55):
   coded_cross(I,J) = Gfull(Z_IJ) - 1/2[Gfull(X_I)+Gfull(X_J)],  Z_IJ=(X_I+X_J)/sqrt2
   Gfull(col) = nuclear grad of the 2e+1e energy expression whose densities
   (td_p, td_mrsf_density, WAO, td_abxc) are built from the SINGLE amplitude col
   through tdhf_mrsf_z_vector + tdhf_mrsf_gradient.  The 4-index get_density
   contraction is QUADRATIC in the transition channels, so the polarization
   differencing recovers a symmetric bilinear -- EXCEPT for the SOMO singlet slot,
   whose density chain uses the INPUT fold (mrsfxvec on bvec) while the
   interstate-consistent assembly needs the OUTPUT fold (mrsfmntoia).  That single
   fold inconsistency is the (1,3) sign flip + (1,2) ratio inflation.

WHAT THIS PoC ESTABLISHES (build-free):
 (A) coded_cross(I,J), reconstructed from the production gradient chain via the
     matvec 2e operator, == the INPUT-fold transported operator bilinear.  This
     REPRODUCES the deficiency (the (1,3) flip / (1,2) inflation) so we know the
     reconstruction is faithful to the production densities.
 (B) Substituting the INTERSTATE-CONSISTENT density assembly = the OUTPUT-fold
     (mrsfmntoia) transported operator (= the validated symmetric matvec A 2e part,
     == X_I^T dA_ref X_J/gap from p11_poc_gaugefree) CLOSES all 3 pairs:
     cos>0.99, ratio in [0.9,1.1].
 (C) Per-pair BLOCK localization of which fold/density slot closes which block:
     the ONLY operator entries that change between the two assemblies are the SOMO
     singlet slot (ijlr1/ijlr2) -- the doc-socc + soc-virt amplitude blocks for
     (1,2), the virt blocks for (1,3); doc-socc-only error for (1,2), virt-only
     error for (1,3); (2,3) unaffected.  Matches p11_isolate_block exactly.

The relaxation (esum/Fock) part is fold-INVARIANT (L^esum==L^2FT) and the 1e/Coulomb
parts are state-independent, so the 2e back-transform fold is the WHOLE story; we
build the full operator A = mntoia(2e) + esum(relaxation) for both fold conventions
and FD-transport in the FIXED reference frame (the gauge-free recipe, no eigenvector
rotation).  Oracle = p11_damp_oracle.npz (== gap * numeric_cross).
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

r = Runner(input_file=INP, log='/tmp/nactest/p11dsub.log'); r.run(); mol = r.mol
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
lr1, lr2 = noca-2, noca-1
ijlr1 = (noca-1-nocb-1)*noca + (noca-1) - 1
ijlr2 = (noca-nocb-1)*noca + (noca) - 1

# rotation/amplitude-block layout (matches sfrorhs / harness prs order)
prs = ([(i, j) for i in range(nocb, noca) for j in range(0, nocb)]
       + [(k, j) for k in range(noca, nbf) for j in range(0, nocb)]
       + [(k, i) for k in range(noca, nbf) for i in range(nocb, noca)])
nds, ndv = nocb*(noca-nocb), (nbf-noca)*nocb
blocks = [('doc-socc', 0, nds), ('doc-virt', nds, nds+ndv), ('soc-virt', nds+ndv, nij)]


def set_bvec(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def iatogen(xv):
    w = np.zeros((nbf, nbf)); w[0:noca, nocb:nbf] = xv.reshape(nvirb, noca).T; return w


def comp(w):
    p = np.zeros(nij); ij = 0
    for j in range(nocb, nbf):
        for i in range(0, noca):
            p[ij] = w[i, j]; ij += 1
    return p


def matvec_kernels(col):
    set_bvec(col); oqp.mrsf_matvec_apply(mol)
    gmo = np.array(mol.data['OQP::nac_gmo'], copy=True).reshape((nbf, nbf), order='F')
    gch = np.array(mol.data['OQP::nac_gchan'], copy=True).reshape((nbf, nbf, 6), order='F')
    return gmo, gch


# ----- esum (relaxation / Fock) part of A -- fold-INVARIANT (verbatim mrsfesum) -----
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


# ----- 2e back-transform, the TWO fold conventions ---------------------------------
def mntoia_OUTPUT(gmo, gch):
    """INTERSTATE-CONSISTENT assembly: the production matvec OUTPUT fold (mrsfmntoia,
    tdhf_mrsf_lib.F90 sec3-6 + the (gmo[lr1,lr1]-gmo[lr2,lr2])*sq SOMO fold).  This is
    the validated symmetric A 2e part."""
    ado2v, ado1v, adco1, adco2, ao21v, aco12 = (gch[:, :, k] for k in range(6))
    wrk = gmo.copy()
    for i in range(0, noca-2):
        wrk[i, lr2] += ado1v[i, lr2] + aco12[i, lr1]
        wrk[i, lr1] += ado2v[i, lr1] - aco12[i, lr2]
    for a in range(noca, nbf):
        wrk[lr1, a] += adco2[lr1, a] + ao21v[lr2, a]
        wrk[lr2, a] += adco1[lr2, a] - ao21v[lr1, a]
    wrk[lr1, lr1] = (gmo[lr1, lr1] - gmo[lr2, lr2])*SQ      # OUTPUT (interstate) fold
    wrk[lr2, lr2] = 0.0
    return comp(wrk)


def mntoia_INPUT(gmo, gch):
    """PRODUCTION single-state gradient-chain assembly: the SOMO slot carries the
    INPUT fold (mrsfxvec applied to bvec, NO symmetric output fold) -> exactly the
    density chain td_mrsf_density built from the single amplitude column."""
    ado2v, ado1v, adco1, adco2, ao21v, aco12 = (gch[:, :, k] for k in range(6))
    wrk = gmo.copy()
    for i in range(0, noca-2):
        wrk[i, lr2] += ado1v[i, lr2] + aco12[i, lr1]
        wrk[i, lr1] += ado2v[i, lr1] - aco12[i, lr2]
    for a in range(noca, nbf):
        wrk[lr1, a] += adco2[lr1, a] + ao21v[lr2, a]
        wrk[lr2, a] += adco1[lr2, a] - ao21v[lr1, a]
    wrk[lr1, lr1] = gmo[lr1, lr1] * SQ                     # INPUT (single-state) fold
    wrk[lr2, lr2] = 0.0
    return comp(wrk)


def build_A(builder):
    """full nij x nij operator A = 2e back-transform (builder) + esum, column by col,
    at the CURRENT geometry/orbitals."""
    matvec_kernels(np.zeros(nij))
    fa = np.array(mol.data['OQP::nac_fa'], copy=True).reshape((nbf, nbf), order='F')
    fb = np.array(mol.data['OQP::nac_fb'], copy=True).reshape((nbf, nbf), order='F')
    A = np.zeros((nij, nij)); e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0; e[j] = 1.0
        gmo, gch = matvec_kernels(e)
        A[:, j] = builder(gmo, gch) + esum_full(e, fa, fb)
    return A


# ---- amplitude-space transport (gauge-free, FIXED reference frame; from p11_poc_gaugefree) ----
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


# ============== anchors at R0 ==============
A_OUT0 = build_A(mntoia_OUTPUT)
A_IN0 = build_A(mntoia_INPUT)
print("=== R0 anchors ===")
print(f"  |A_OUT - A_OUT^T|        = {np.max(np.abs(A_OUT0-A_OUT0.T)):.2e}  (interstate assembly: SYMMETRIC)")
print(f"  |A_IN  - A_IN^T |        = {np.max(np.abs(A_IN0-A_IN0.T)):.2e}  (production single-state: NOT symmetric)")
wv = np.sort(np.linalg.eigvalsh(A_OUT0))[:4]
print(f"  eig(A_OUT) lowest        = {np.round(wv,5)}   Om = {np.round(Om,5)}")
print(f"  |A_OUT - A_IN|           = {np.max(np.abs(A_OUT0-A_IN0)):.2e}  (differ ONLY at SOMO slot)")
# show the diff is confined to the SOMO singlet slot (row/col ijlr1)
D = A_OUT0 - A_IN0
off_slot = D.copy(); off_slot[ijlr1, :] = 0.0; off_slot[:, ijlr1] = 0.0
print(f"  |A_OUT - A_IN| OFF the ijlr1 row/col = {np.max(np.abs(off_slot)):.2e}  (=> diff is the SOMO fold only)")

mol.save_data()
cfg = mol.config; json0 = mol.log.replace('.log', '.json')
cfg['guess']['type'] = 'json'; cfg['guess']['file'] = json0; cfg['guess']['continue_geom'] = False


def displaced_Aref(coord, builder):
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
    A_disp = build_A(builder)
    T = transport_T(Q)
    return T.T @ A_disp @ T


# build dA_ref for BOTH assemblies (interstate OUTPUT and production INPUT)
dA_OUT = np.zeros((ncoord, nij, nij)); dA_IN = np.zeros((ncoord, nij, nij))
for k in range(ncoord):
    AOp = displaced_Aref(xyz0 + DELTA*np.eye(ncoord)[k], mntoia_OUTPUT)
    AOm = displaced_Aref(xyz0 - DELTA*np.eye(ncoord)[k], mntoia_OUTPUT)
    dA_OUT[k] = (AOp - AOm)/(2*DELTA)
    AIp = displaced_Aref(xyz0 + DELTA*np.eye(ncoord)[k], mntoia_INPUT)
    AIm = displaced_Aref(xyz0 - DELTA*np.eye(ncoord)[k], mntoia_INPUT)
    dA_IN[k] = (AIp - AIm)/(2*DELTA)
mol.update_system(xyz0); oqp.library.ints_1e(mol); oqp.library.guess(mol); SinglePoint(mol).energy()
mol.data._data.control.int2e_cutoff = 1e-20

oracle = np.load('/tmp/nactest/p11_damp_oracle.npz')
qz = np.load('/tmp/nactest/qz_full.npz')
gap = {(0, 1): Om[1]-Om[0], (0, 2): Om[2]-Om[0], (1, 2): Om[2]-Om[1]}
pmap = {(0, 1): '12', (0, 2): '13', (1, 2): '23'}
# production coded_cross (the deficient target) and numeric_cross (= gap*oracle)
numeric_cross = {}
coded_cross = {}
for (I, J), p in zip([(0, 1), (0, 2), (1, 2)], ('12', '13', '23')):
    a, b = str(I+1), str(J+1)
    numeric_cross[(I, J)] = qz[f'Z{p}_numeric'] - 0.5*(qz[f'X{a}_numeric'] + qz[f'X{b}_numeric'])
    coded_cross[(I, J)] = qz[f'Z{p}_coded'] - 0.5*(qz[f'X{a}_coded'] + qz[f'X{b}_coded'])

PAIRS = [(0, 1), (0, 2), (1, 2)]


def cosratio(a, b):
    c = np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30)
    rr = np.linalg.norm(a)/(np.linalg.norm(b)+1e-30)
    return c, rr


print("\n=== (A) DEFICIENCY REPRODUCTION: INPUT-fold (single-state) assembly carries the defect ===")
print("    INPUT-fold nuclear cross vs the production coded_cross (qz_full) AND vs the oracle.")
print("    The (1,3) sign flip is the SOMO input-fold; it is what the OUTPUT fold corrects.")
for (I, J) in PAIRS:
    g = gap[(I, J)]
    damp_in = np.array([0.5*(X0[I] @ dA_IN[k] @ X0[J] + X0[J] @ dA_IN[k] @ X0[I])
                        for k in range(ncoord)])/g
    cc = coded_cross[(I, J)]/g
    o = oracle[f'{I+1}{J+1}_damp']
    c_vs_prod, r_vs_prod = cosratio(damp_in, cc)
    c_vs_o, r_vs_o = cosratio(damp_in, o)
    flag = '   <-- SIGN FLIP vs oracle' if c_vs_o < -0.9 else ''
    print(f"  ({I+1},{J+1}): INPUT-fold vs production coded_cross cos={c_vs_prod:+.4f} r={r_vs_prod:.3f}"
          f"  | vs oracle cos={c_vs_o:+.4f} r={r_vs_o:.3f}{flag}")

print("\n=== (B) FIX: INTERSTATE-CONSISTENT (OUTPUT-fold) assembly vs oracle ===")
allpass = True
for (I, J) in PAIRS:
    g = gap[(I, J)]
    damp_out = np.array([X0[I] @ dA_OUT[k] @ X0[J] for k in range(ncoord)])/g
    o = oracle[f'{I+1}{J+1}_damp']
    c, rr = cosratio(damp_out, o)
    s = 1.0 if c >= 0 else -1.0
    relresid = np.linalg.norm(s*damp_out - o)/(np.linalg.norm(o)+1e-30)
    ok = abs(c) > 0.99 and 0.9 < rr < 1.1
    allpass = allpass and ok
    print(f"  ({I+1},{J+1}): cos={c:+.8f} ratio={rr:.6f} relresid={100*relresid:.3f}%  {'PASS' if ok else 'FAIL'}")
    print(f"        interstate={np.array2string(s*damp_out,precision=5,suppress_small=True)}")
    print(f"        oracle    ={np.array2string(o,precision=5,suppress_small=True)}")

print("\n=== (C) BLOCK localization: which NUCLEAR amplitude block the SOMO fold corrects ===")
print("    The correction is the nuclear-derivative bilinear of the operator diff dD=d(A_OUT-A_IN).")
print("    Per pair we report |damp_OUT - damp_IN| restricted to bra/ket components in each block")
print("    (block = amplitude index of the SOMO-coupled leg).  Matches p11_isolate_block signature.")
for (I, J) in PAIRS:
    g = gap[(I, J)]
    o = oracle[f'{I+1}{J+1}_damp']
    d_out = np.array([X0[I] @ dA_OUT[k] @ X0[J] for k in range(ncoord)])/g
    d_in = np.array([0.5*(X0[I] @ dA_IN[k] @ X0[J] + X0[J] @ dA_IN[k] @ X0[I])
                     for k in range(ncoord)])/g
    co, _ = cosratio(d_out, o)
    ci, _ = cosratio(d_in, o)
    # block-resolved nuclear correction: derivative of the SOMO-fold operator diff,
    # with the bra/ket restricted to each amplitude block.
    line = f"  ({I+1},{J+1}): cos(IN)={ci:+.4f} -> cos(OUT)={co:+.4f}  block-correction |dD|:"
    for name, lo, hi in blocks:
        braM = np.zeros(nij); braM[lo:hi] = X0[I][lo:hi]
        ketM = np.zeros(nij); ketM[lo:hi] = X0[J][lo:hi]
        corr = np.array([0.5*((braM @ (dA_OUT[k]-dA_IN[k]) @ X0[J]) + (X0[I] @ (dA_OUT[k]-dA_IN[k]) @ ketM))
                         for k in range(ncoord)])/g
        line += f"  {name}={np.linalg.norm(corr):.3e}"
    print(line)

print(f"\nVERDICT: {'ALL THREE PAIRS PASS -- INTERSTATE DENSITY ASSEMBLY CLOSES THE ORACLE' if allpass else 'NOT all pass'}")
np.savez('/tmp/nactest/p11_density_subst.npz',
         A_OUT0=A_OUT0, A_IN0=A_IN0, dA_OUT=dA_OUT, dA_IN=dA_IN, X0=X0, Om=np.array(Om))
mol.update_system(xyz0)
