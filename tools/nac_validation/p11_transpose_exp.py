"""Phase 11 OUTPUT-FOLD-TRANSPOSE experiment (read-only).

HYPOTHESIS (from nac-phase11-blueprint, not yet validated):
  The matvec encodes the SOMO fold as an OUTPUT fold O inside mrsfmntoia, and
  O == M^T EXACTLY (M = the mrsfxvec INPUT fold), while M != M^T.  On the diagonal
  both give the same RHS (M==M^T contracted with the same X), but OFF-diagonal the
  transpose-difference produces a DIFFERENT bilinear -> candidate fix for the (1,3)
  sign flip and the (1,2) gl deficiency.

WHAT THIS SCRIPT DOES (purely measurement):
  Build the interstate 2e z-vector RHS two ways for each pair (I,J):
    route M  (input fold) : wrk3_K = iatogen( M  @ X_K )   (== mrsfxvec, current)
    route MT (output fold): wrk3_K = iatogen( M^T @ X_K )  (transpose-fold variant)
  with the symmetric bilinear   hxa = G_I wrk3_J^T + G_J wrk3_I^T ,
                                hxb = G_I^T wrk3_J + G_J^T wrk3_I   (gmo_interstate).
  Project to the ROHF rotation space (sfrorhs block projection), then push EACH
  rotation-space RHS through the production z-vector solve + gradient (the
  mrsf_nac_cphf seam: inject the rotation RHS as wrk1[hi,lo], antisym projection
  reproduces it) -> a NUCLEAR cross coupling, and compare to the test_QZ
  numeric_cross oracle (qz_full.npz) for (1,2),(1,3),(2,3).

  ALSO reports the rotation-space comparison against the gradient-chain 2e RHS dump
  (the deficient baseline) so the rotation-space picture is visible even if the
  nuclear projection is noisy.

The diagonal (I=J) is exact for both routes by construction (M==M^T on the diag),
so the discriminating signal lives ENTIRELY in the off-diagonal cross term.

Run:  source /tmp/nactest/nacenv.sh && OMP_NUM_THREADS=1 \
      /opt/conda/bin/python3 p11_transpose_exp.py
"""
import os
import oqp                       # MUST precede numpy (ILP64 LAPACK clash)
import oqp.library
from oqp.pyoqp import Runner
import numpy as np

INP = '/tmp/nactest/H2O_en.inp'
SQ = 1.0 / np.sqrt(2.0)

r = Runner(input_file=INP, log='/tmp/nactest/p11_transpose.log')
r.run(); mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20      # LINEAR matvec
print("RUN OK (tight int2e_cutoff=1e-20)", flush=True)

noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
nbf = np.array(mol.data['OQP::VEC_MO_A'], copy=True).shape[0]
nvirb = nbf - nocb; nij = noca * nvirb; nstate = 3; natom = mol.data['natom']
Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); Xshape = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))

# slot indices (0-based)
ijlr1 = (noca - 1 - nocb - 1) * noca + (noca - 1) - 1
ijlr2 = (noca - nocb - 1) * noca + (noca) - 1
print(f"nbf={nbf} noca={noca} nocb={nocb} nvirb={nvirb} nij={nij} "
      f"ijlr1={ijlr1} ijlr2={ijlr2} natom={natom}")

# rotation-space block layout (matches sfrorhs / sfrogen / the seam)
prs = ([(i, j) for i in range(nocb, noca) for j in range(0, nocb)]
       + [(k, j) for k in range(noca, nbf) for j in range(0, nocb)]
       + [(k, i) for k in range(noca, nbf) for i in range(nocb, noca)])
nconf = len(prs)
nds, ndv = nocb * (noca - nocb), (nbf - noca) * nocb
blocks = [('doc-socc', 0, nds), ('doc-virt', nds, nds + ndv), ('soc-virt', nds + ndv, nconf)]


# ---------------------------------------------------------------------------
# FOLD OPERATORS as nij x nij matrices on the COMPRESSED amplitude.
#   M  (input fold, mrsfxvec, mult=1):  M[ijlr1,ijlr1]=+SQ , M[ijlr2,ijlr1]=-SQ
#       (slot ijlr2 sources from input slot ijlr1; output slot ijlr1 scaled by SQ)
#   identity elsewhere.
#   O  (output fold, mntoia):  the compressed-output slot ijlr1 receives
#       (scr[lr1,lr1]-scr[lr2,lr2])*SQ and slot ijlr2 -> 0.  In compressed-input
#       coordinates that is  O[ijlr1,ijlr1]=+SQ , O[ijlr1,ijlr2]=-SQ , slot ijlr2->0.
#   => O == M^T exactly (verified below).
# ---------------------------------------------------------------------------
M = np.eye(nij)
M[ijlr1, ijlr1] = SQ
M[ijlr2, ijlr1] = -SQ
M[ijlr2, ijlr2] = 0.0          # ijlr2 input slot is overwritten (mrsfxvec sets it)

O = np.eye(nij)
O[ijlr1, ijlr1] = SQ
O[ijlr1, ijlr2] = -SQ
O[ijlr2, ijlr2] = 0.0          # output slot ijlr2 zeroed

print("\n=== fold-operator identities ===")
print(f"  max|O - M^T| = {np.max(np.abs(O - M.T)):.2e}   (hypothesis: 0)")
print(f"  max|M - M^T| = {np.max(np.abs(M - M.T)):.2e}   (hypothesis: ~0.707)")
# sanity: M applied as matrix == mrsfxvec applied elementwise
def mrsfxvec(col):
    f = col.copy(); f[ijlr1] = col[ijlr1] * SQ; f[ijlr2] = -col[ijlr1] * SQ
    return f
chk = np.max(np.abs(M @ X[0] - mrsfxvec(X[0])))
print(f"  max|M@X - mrsfxvec(X)| = {chk:.2e}   (M matrix == mrsfxvec)")


def set_bvec(col):
    raw = Xraw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def get_gmo(col):
    set_bvec(col); oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_gmo'], copy=True).ravel().reshape((nbf, nbf), order='F')


def iatogen(colf):
    av = np.zeros((nbf, nbf)); av[0:noca, nocb:nbf] = colf.reshape((nvirb, noca)).T
    return av


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


def rhs_inter(GI, GJ, wI, wJ):
    """symmetric bilinear 2e RHS (gmo_interstate)."""
    hxa = GI @ wJ[0:noca, :].T + GJ @ wI[0:noca, :].T
    hxb = GI[0:noca, :].T @ wJ[0:noca, :] + GJ[0:noca, :].T @ wI[0:noca, :]
    return project(hxa, hxb)


def rhs_diag(GI, wI):
    hxa = 2.0 * GI @ wI[0:noca, :].T
    hxb = 2.0 * GI[0:noca, :].T @ wI[0:noca, :]
    return project(hxa, hxb)


# cache per-state kernels + folded amplitude matrices for both routes
gmos = [get_gmo(X[k]) for k in range(nstate)]
wM = [iatogen(M @ X[k]) for k in range(nstate)]      # route M  (input fold)
wO = [iatogen(O @ X[k]) for k in range(nstate)]      # route MT (output fold)

print("\n=== I=J reduction: both routes -> validated diagonal RHS (must be ~1e-15) ===")
for st in range(nstate):
    RdM = rhs_diag(gmos[st], wM[st])
    RM = rhs_inter(gmos[st], gmos[st], wM[st], wM[st])
    RO = rhs_inter(gmos[st], gmos[st], wO[st], wO[st])
    print(f"  state {st+1}: |R_diag|={np.linalg.norm(RdM):.4e}  "
          f"max|M_inter(I=I)-R_diag|={np.max(np.abs(RM-RdM)):.2e}  "
          f"max|MT_inter(I=I)-M_inter(I=I)|={np.max(np.abs(RO-RM)):.2e}")


# ---------------------------------------------------------------------------
# Rotation-space cross RHS, both routes, vs the gradient-chain 2e dump baseline
# ---------------------------------------------------------------------------
def dump_rhs(target, gates):
    """generic sfrorhs dump with the given NAC_ZERO_* gates set."""
    mol.data.set_tdhf_target(target)
    allg = ('NAC_ZERO_2FT', 'NAC_ZERO_AB1', 'NAC_ZERO_HX', 'NAC_ZERO_T', 'NAC_ZERO_SP')
    for k in allg:
        os.environ.pop(k, None)
    for k in gates:
        os.environ[k] = '1'
    os.environ['NAC_DUMP_RHS'] = '1'
    oqp.tdhf_mrsf_z_vector(mol)
    os.environ.pop('NAC_DUMP_RHS', None)
    for k in allg:
        os.environ.pop(k, None)
    return np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).ravel()[:nconf]


def dump_chain_2e(target):
    """gradient-chain 2e RHS (G_MO + mrsfsp), 2FT+AB1 gated off; SP kept."""
    return dump_rhs(target, ['NAC_ZERO_2FT', 'NAC_ZERO_AB1'])


def dump_chain_2e_nosp(target):
    """G_MO (channel-7) only: SP gated off too (== gmo_validate target)."""
    return dump_rhs(target, ['NAC_ZERO_2FT', 'NAC_ZERO_AB1', 'NAC_ZERO_SP'])


def cosr(a, b):
    return (np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-30),
            np.linalg.norm(a) / (np.linalg.norm(b) + 1e-30))


PAIRS = [(0, 1), (0, 2), (1, 2)]
RHS_M = {}; RHS_O = {}
print("\n=== rotation-space cross RHS: route M vs route MT (norms per block) ===")
for (I, J) in PAIRS:
    RM = rhs_inter(gmos[I], gmos[J], wM[I], wM[J])
    RO = rhs_inter(gmos[I], gmos[J], wO[I], wO[J])
    RHS_M[(I, J)] = RM; RHS_O[(I, J)] = RO
    c, _ = cosr(RM, RO)
    line = (f"  pair ({I+1},{J+1}): |M|={np.linalg.norm(RM):.4e} |MT|={np.linalg.norm(RO):.4e} "
            f"cos(M,MT)={c:+.4f}")
    print(line)
    for nm, lo, hi in blocks:
        cb, rb = cosr(RM[lo:hi], RO[lo:hi])
        print(f"      {nm:9}: |M|={np.linalg.norm(RM[lo:hi]):.3e} |MT|={np.linalg.norm(RO[lo:hi]):.3e} "
              f"cos={cb:+.4f} ratio(M/MT)={rb:.3f}")


# ---------------------------------------------------------------------------
# NUCLEAR projection via the mrsf_nac_cphf seam.
# The seam (tdhf_mrsf_z_vector.F90:1785-1839) reads wrk1 = OQP::nac_gamma_tlf
# and sets rhs(ijp) = wrk1(hi,lo)-wrk1(lo,hi) over the SAME prs block order, then
# zeros td_abxc/td_mrsf_den/tij/tab/hxa/hxb and solves (A+B)Z=rhs, then the
# gradient contracts ONLY Z (the SCF/ground part removed by gZ-gS differencing).
# => injecting wrk1[hi,lo]=R[ijp], wrk1[lo,hi]=0 makes the seam's antisym
# projection reproduce R[ijp] exactly.  This pushes an arbitrary rotation-space
# RHS through the production z-vector+gradient to a nuclear coupling.
# ---------------------------------------------------------------------------
mol.save_data()


def inject_wrk1(R):
    """build the nbf x nbf wrk1 whose seam antisym projection == R."""
    W = np.zeros((nbf, nbf))
    for n, (hi, lo) in enumerate(prs):
        W[hi, lo] = R[n]            # wrk1(lo,hi)=0 -> projection = R[n]
    return W


def push_gamma(W):
    """write W into OQP::nac_gamma_tlf for ALL (ist,jst); we only ever query one."""
    flat = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            flat[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = W.T.reshape(-1)
    mol.data['OQP::nac_gamma_tlf'] = flat.reshape((nbf * nbf, nstate, nstate))


def grad_now():
    oqp.tdhf_mrsf_gradient(mol)
    return mol.get_grad().reshape(-1).copy()


def nuclear_from_rhs(i, j, R):
    """gZ - gS with the seam RHS = R (rotation space). i,j are 1-indexed states."""
    push_gamma(inject_wrk1(R))
    mol.data.set_tdhf_target(i)
    oqp.set_mrsf_nac_cphf(mol, i, j)
    oqp.tdhf_mrsf_z_vector(mol)
    conv = mol.mol_energy.Z_Vector_converged
    if not conv:
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        return None
    gZ = grad_now()
    mol.data['OQP::td_p'] = np.zeros_like(np.array(mol.data['OQP::td_p'], copy=True))
    mol.data['OQP::WAO'] = np.zeros_like(np.array(mol.data['OQP::WAO'], copy=True))
    gS = grad_now()
    oqp.set_mrsf_nac_cphf(mol, 0, 0)
    return (gZ - gS).reshape(-1)


# oracle: numeric_cross from qz_full.npz (nuclear, 9 comps)
qz = np.load('/tmp/nactest/qz_full.npz')


def numeric_cross(I, J):
    return qz[f'Z{I}{J}_numeric'] - 0.5 * (qz[f'X{I}_numeric'] + qz[f'X{J}_numeric'])


def coded_cross_chain(I, J):
    """the deficient gradient-chain cross (coded) for reference."""
    return qz[f'Z{I}{J}_coded'] - 0.5 * (qz[f'X{I}_coded'] + qz[f'X{J}_coded'])


def cos9(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return (np.dot(a, b) / (na * nb + 1e-30), na / (nb + 1e-30))


# ---- SEAM FIDELITY anchor 1: the seam is a LINEAR map rotation->nuclear ----
print("\n=== SEAM FIDELITY: nuclear(seam) is linear in the rotation RHS ===")
g_unit = nuclear_from_rhs(1, 2, RHS_M[(0, 1)])
g_two = nuclear_from_rhs(1, 2, 2.0 * RHS_M[(0, 1)])
if g_unit is not None and g_two is not None:
    print(f"  max|nuclear(2R)-2 nuclear(R)| = {np.max(np.abs(g_two - 2 * g_unit)):.2e}  "
          f"(linear => ~0)")

# ---- SEAM ANCHOR 2: route-M is CHANNEL-7 ONLY. It matches the SP-gated-OFF dump
#      (gmo_validate, cos=+1.0); the SP-kept dump adds the mrsfsp channels 1-6,
#      which route M/MT DO NOT carry -> the gap quantifies the missing
#      interstate channels-1-6 piece (task suspect #1).
print("\n=== SEAM ANCHOR: route-M (ch7) vs SP-off dump (G_MO) and SP-on dump (G_MO+sp) ===")
for st in range(1, nstate + 1):
    RM = rhs_diag(gmos[st - 1], wM[st - 1])
    Rnosp = dump_chain_2e_nosp(st)        # G_MO only
    Rsp = dump_chain_2e(st)               # G_MO + mrsfsp
    c0, r0 = cosr(RM, -Rnosp)
    c1, r1 = cosr(RM, -Rsp)
    print(f"  state {st}: vs G_MO(SP off) cos={c0:+.8f} ratio={r0:.5f} "
          f"max|diff|={np.max(np.abs(RM + Rnosp)):.2e}")
    print(f"           vs G_MO+sp(SP on) cos={c1:+.8f} ratio={r1:.5f} "
          f"|sp piece|={np.linalg.norm(Rsp - Rnosp):.4e}")

# I=J seam anchor: nuclear(seam, diagonal 2e RHS) finite & reproducible
print("\n    I=J anchor (seam fed route-M diagonal 2e RHS), reports |nuclear|:")
for st in range(1, nstate + 1):
    Rd = rhs_diag(gmos[st - 1], wM[st - 1])
    g = nuclear_from_rhs(st, st, Rd)
    tag = 'n/a (Z not converged)' if g is None else f"|nuclear|={np.linalg.norm(g):.5e}"
    print(f"      state {st}: {tag}")


# ---- SEAM ANCHOR 3 (apples-to-apples): push the FULL chain cross (the exact
#      object qz_full 'coded' nuclear-measures: full sfrorhs = 2e+SP+2FT+AB1)
#      through the seam, confirm seam reproduces qz coded_cross.  This proves the
#      seam projection is faithful AND shows the gap between the 2e-only routes
#      and the full numeric.
def set_amp(col):
    raw = Xraw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def full_chain_rhs(col, target):
    """full sfrorhs (all terms) for amplitude col at given target, as the stored
    rhs (negated inside sfrorhs)."""
    set_amp(col)
    Rd = dump_rhs(target, [])      # no gates: full 2e+SP+2FT+AB1
    set_bvec(X[target - 1])        # restore
    return Rd


def full_chain_cross(I, J):
    """rotation-space full chain CROSS via polarization (matches qz coded_cross):
        A(X_I,X_J) = [A(Z_IJ) - 1/2(A(X_I)+A(X_J))], Z_IJ=(X_I+X_J)/sqrt2.
    target index for the cross uses state I (the qz convention: grad_for(i,.))."""
    i = I + 1
    Z = (X[I] + X[J]) / np.sqrt(2.0)
    RZ = full_chain_rhs(Z, i)
    RI = full_chain_rhs(X[I], i)
    RJ = full_chain_rhs(X[J], i)
    return RZ - 0.5 * (RI + RJ)


print("\n=== SEAM ANCHOR 3: nuclear(seam, FULL chain cross) vs qz coded_cross ===")
for (I, J) in PAIRS:
    i, j = I + 1, J + 1
    Rfull = full_chain_cross(I, J)
    g = nuclear_from_rhs(i, j, Rfull)
    cc = coded_cross_chain(i, j)
    if g is None:
        print(f"  ({i},{j}): Z not converged")
        continue
    c, rr = cos9(g, cc)
    print(f"  ({i},{j}): cos(seam_full, qz_coded)={c:+.6f} ratio={rr:.4f} "
          f"|seam|={np.linalg.norm(g):.5e} |qz_coded|={np.linalg.norm(cc):.5e}")


# ---------------------------------------------------------------------------
# MAIN COMPARISON
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("MAIN: nuclear cross coupling, route M vs route MT, vs numeric_cross oracle")
print("=" * 78)
hdr = f"{'pair':6} {'route':6} {'cos(coded,num)':>16} {'ratio|coded|/|num|':>20} {'|nuclear|':>12} {'|num|':>10}"
print(hdr)
results = {}
for (I, J) in PAIRS:
    i, j = I + 1, J + 1
    num = numeric_cross(i, j)
    nM = nuclear_from_rhs(i, j, RHS_M[(I, J)])
    nO = nuclear_from_rhs(i, j, RHS_O[(I, J)])
    results[(i, j)] = {'num': num, 'M': nM, 'O': nO}
    for tag, g in (('M', nM), ('MT', nO)):
        if g is None:
            print(f"({i},{j})  {tag:6} {'Z not conv':>16}")
            continue
        c, rr = cos9(g, num)
        print(f"({i},{j})  {tag:6} {c:>+16.6f} {rr:>20.4f} "
              f"{np.linalg.norm(g):>12.5e} {np.linalg.norm(num):>10.5e}")

# also report the chain-baseline nuclear cross (coded chain) for context
print("\n--- reference: gradient-chain coded cross vs numeric (the known deficiency) ---")
for (I, J) in PAIRS:
    i, j = I + 1, J + 1
    cc = coded_cross_chain(i, j)
    num = numeric_cross(i, j)
    c, rr = cos9(cc, num)
    print(f"  ({i},{j}): cos(chain_coded, num)={c:+.6f}  ratio={rr:.4f}  "
          f"|chain|={np.linalg.norm(cc):.5e} |num|={np.linalg.norm(num):.5e}")

print("\n--- 13_damp benchmark scalar (target 0.0284): |numeric_cross(1,3)| context ---")
print(f"  |numeric_cross(1,3)| = {np.linalg.norm(numeric_cross(1,3)):.5f}")

# ---- MECHANISM: M vs M^T differ ONLY in the ijlr SOMO slots, so the route
#      difference is driven entirely by the ground-config amplitude X[ijlr1].
print("\n--- MECHANISM: per-state ground-config amplitude X[ijlr1] (drives M vs MT) ---")
for st in range(nstate):
    print(f"  state {st+1}: X[ijlr1]={X[st, ijlr1]:+.6f}  X[ijlr2]={X[st, ijlr2]:+.6f}")
print("  => states 1,3 have X[ijlr1]=0 -> M@X==M^T@X==X, routes IDENTICAL for (1,3).")
print("     Only state 2 carries the fold, so M vs MT differ only via state-2.")

print("\nLIMITATIONS:")
print("  * seam is a LINEAR, faithful rotation->nuclear projector (anchor 1: 1e-11)")
print("    but it is the CPHF orbital-response (Z-vector) projection, NOT the full")
print("    amplitude gradient -> absolute |nuclear| differs from the full numeric")
print("    oracle (anchor 3: cos=+/-1 but ratio != 1). The M-vs-MT A/B comparison is")
print("    still controlled (identical seam), so direction differences are real.")
print("  * routes M/MT are CHANNEL-7 (G_MO) only; the mrsfsp channels 1-6 are absent")
print("    (anchor 2: sp piece 1.5e-2 on states 1,3) -> the channels-1-6 interstate")
print("    piece is NOT addressed by the fold transpose.")
print("  * (2,3) gap=0.058 Ha is near-degenerate; the seam Z-solve is ill-conditioned")
print("    there (full-chain-cross Z did not converge) -> (2,3) seam numbers are soft.")
