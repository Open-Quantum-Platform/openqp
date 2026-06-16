"""Phase 11: frozen-Fock matvec evaluator + orbital-rotation FD of X_I^T(d A)X_J.

Step A (correctness): the standalone matvec mrsf_matvec_apply must reproduce
A.X_K so that X_K . (A X_K) = Omega_K (the excitation energy) for every state.

Step B (frozen-Fock orbital gradient): rotate the ROHF orbitals in a single
occ/virt plane, hold the AO Fock fixed (mrsf_matvec_apply rebuilds C^T F_AO C),
and finite-difference z^T A z to get R^matvec_pq = d(z^T A z)/d theta_pq for
z = X2 (diagonal). The diagonal was expected to agree with sfrorhs by HF.

CONCLUSION (2026-06-16, see PHASE11_fd_findings.md): the diagonal does NOT match
for state 2. ROOT CAUSE = mrsfesum applies the open-shell 1/sqrt2 SOMO fold on the
OUTPUT side, while the production sfrorhs path folds on the INPUT (mrsfxvec). The
fold's orbital-rotation derivative, driven by xlr~0.998, is a large artifact that
swamps every alpha block (zeroing SOMO amplitudes collapses |R| by 78-545x per
block). FD of X^T A_TDA X is therefore the WRONG operator for the rotation gradient
on the MRSF SOMO sector. Conventions that ARE correct: AXIS='col', the (i+,a-)
amplitude layout, the Givens generator sign vs sfrogen. sfrorhs itself is NOT the
bare FD of X^T A X: it = -[ hxa/hxb 2e-projection + 2FT relaxation + ab1 B-matrix ].
The input-side refold of X does NOT fix it (double-folds with mrsfesum's output
fold); a faithful Python-only comparison needs the matvec to export the full-MO
G_MO (a Fortran add, out of scope here).
"""
import oqp                       # MUST precede numpy: numpy's ILP64 LAPACK
from oqp.pyoqp import Runner      # otherwise intercepts oqp's DSYEVD (Huckel fails)
import numpy as np

INP = '/tmp/nactest/H2O_en.inp'
TH = 1e-3

r = Runner(input_file=INP, log='/tmp/nactest/mvrhs.log')
r.run()
mol = r.mol
print("RUN OK", flush=True)

nstate = 3
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2                      # MRSF: 2 SOMOs
C0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = C0.shape[0]
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
nvirb = nbf - nocb
nij = noca * nvirb

Xraw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = Xraw.shape
X = Xraw.reshape(-1).reshape((nstate, nij))
Om = np.array(mol.data['OQP::td_energies'], copy=True).ravel()   # excitation energies

print(f"nbf={nbf} noca={noca} nocb={nocb} nvirb={nvirb} nij={nij}")
print("Omega (Hartree):", np.array2string(Om, precision=6))


def set_bvec(col):
    raw = Xraw.copy().reshape(-1)
    raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)


def Ax(col):
    """A . col  via the frozen-Fock standalone matvec (current orbitals)."""
    set_bvec(col)
    oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()


# ---- Step A: eigenvalue reproduction (validates the standalone matvec) ----
print("\n=== Step A: X_K . (A X_K) vs Omega_K  (matvec correctness) ===")
okA = True
for k in range(nstate):
    xk = X[k]
    axk = Ax(xk)
    nn = np.dot(xk, xk)
    expk = np.dot(xk, axk) / nn
    resid = np.linalg.norm(axk - Om[k]*xk) / np.sqrt(nn)   # eigenvector residual
    err = abs(expk - Om[k])
    print(f"  state {k+1}: <X|AX>/<X|X> = {expk:.6f}  Omega = {Om[k]:.6f}  "
          f"|diff| = {err:.2e}   ||AX-OmX||/||X|| = {resid:.2e}")
    if resid > 1e-3:
        okA = False
print("Step A:", "PASS -- standalone matvec reproduces the spectrum" if okA
      else "FAIL -- matvec mismatch, debug before FD")

# ---- Step B: frozen-Fock orbital-rotation FD -> R^matvec, vs sfrorhs ----
# SPIN-RESOLVED ROHF parameterization, from sfrogen (tdhf_sf_lib.F90:487):
#   doc-socc -> beta density only      (alpha doc,socc both occupied = redundant)
#   doc-virt -> alpha AND beta
#   socc-virt-> alpha density only      (beta socc,virt both empty = redundant)
# pair = (a=high orbital, b=low orbital, spin set to rotate).  Order matches sfrogen.
prs = ([(i, j, 'b') for i in range(nocb, noca) for j in range(0, nocb)]      # doc-socc
       + [(k, j, 'ab') for k in range(noca, nbf) for j in range(0, nocb)]    # doc-virt
       + [(k, i, 'a') for k in range(noca, nbf) for i in range(nocb, noca)]) # socc-virt
nconf = len(prs)
print(f"\n=== Step B: spin-resolved frozen-Fock FD R^matvec  (nconf={nconf}) ===")


def givens(C, a, b, th, axis):
    Cr = C.copy(); c, s = np.cos(th), np.sin(th)
    if axis == 'row':                 # MOs are rows (Python is Fortran^T)
        Cr[a, :] = c*C[a, :] + s*C[b, :]; Cr[b, :] = -s*C[a, :] + c*C[b, :]
    else:                             # MOs are columns
        Cr[:, a] = c*C[:, a] + s*C[:, b]; Cr[:, b] = -s*C[:, a] + c*C[:, b]
    return Cr


def Rmatvec(col, axis):
    R = np.zeros(nconf)
    for n, (a, b, spin) in enumerate(prs):
        def f(th):
            Ca = givens(C0, a, b, th, axis) if spin in ('a', 'ab') else C0
            Cb = givens(C0b, a, b, th, axis) if spin in ('b', 'ab') else C0b
            mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(Ca)
            mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(Cb)
            return float(np.dot(col, Ax(col)))
        R[n] = (f(TH) - f(-TH)) / (2*TH)
    mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(C0)   # restore
    mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(C0b)
    return R


import os
TARGET = 2
nds, ndv = nocb*(noca-nocb), (nbf-noca)*nocb
blocks = [('doc-socc', 0, nds), ('doc-virt', nds, nds+ndv), ('soc-virt', nds+ndv, nconf)]

# NOTE on phase: the Davidson eigenvector sign is arbitrary BETWEEN runs. Within a
# single run, R^matvec (uses this Python X copy) and R^sfrorhs (uses the Fortran
# td_bvec_mo) see the SAME eigenvector, so their relative cos sign is consistent
# within a run; only the OVERALL sign (cos vs -cos) varies run-to-run because
# sfrorhs is X-linear (odd) while the FD of X^T A X is X-even. The magnitude RATIOS
# and the block-collapse factors (the load-bearing results) are phase-independent.
# Do NOT renormalize X here (it would desync from the Fortran td_bvec_mo sfrorhs
# reads). The report below prints |cos| alongside the signed cos for stability.

# FD of X^T A_TDA X w.r.t. the ROHF rotation (only the 'col' axis is correct).
AXIS = 'col'
Rmv = Rmatvec(X[TARGET-1], AXIS)


def dump_sfrorhs(gate=None):
    """Run the production z-vector RHS builder once, optionally with NAC gates
    set, and return the rotation-space rhs (length nconf)."""
    mol.data.set_tdhf_target(TARGET)
    for k in ('NAC_ZERO_2FT', 'NAC_ZERO_AB1', 'NAC_ZERO_HX', 'NAC_ZERO_T'):
        os.environ.pop(k, None)
    if gate:
        for k in gate:
            os.environ[k] = '1'
    os.environ['NAC_DUMP_RHS'] = '1'
    oqp.tdhf_mrsf_z_vector(mol)
    os.environ.pop('NAC_DUMP_RHS')
    for k in ('NAC_ZERO_2FT', 'NAC_ZERO_AB1', 'NAC_ZERO_HX', 'NAC_ZERO_T'):
        os.environ.pop(k, None)
    return np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).ravel()[:nconf]


# Full production rhs and several GATED variants to isolate what the TDA matvec
# FD can and cannot reproduce.
Rsfr_full = dump_sfrorhs(None)
Rsfr_no2ft = dump_sfrorhs(['NAC_ZERO_2FT'])
Rsfr_noab1 = dump_sfrorhs(['NAC_ZERO_AB1'])
Rsfr_aonly = dump_sfrorhs(['NAC_ZERO_2FT', 'NAC_ZERO_AB1'])
print(f"\n|R^sfrorhs full|   = {np.linalg.norm(Rsfr_full):.4e}")
print(f"|R^sfrorhs no2FT|  = {np.linalg.norm(Rsfr_no2ft):.4e}")
print(f"|R^sfrorhs noAB1|  = {np.linalg.norm(Rsfr_noab1):.4e}")
print(f"|R^sfrorhs A-only| = {np.linalg.norm(Rsfr_aonly):.4e}", flush=True)


def cosr(a, b):
    return (np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30),
            np.linalg.norm(a)/(np.linalg.norm(b)+1e-30))


def report(label, Rsfr):
    c, rr = cosr(Rmv, Rsfr)
    print(f"\n[{label}]  overall |cos|={abs(c):.4f} ratio(mv/sfr)={rr:.3f}")
    for name, lo, hi in blocks:
        bc, br = cosr(Rmv[lo:hi], Rsfr[lo:hi])
        print(f"    {name:9} |cos|={abs(bc):.4f} ratio={br:.3f}")


for label, R in (('sfr-FULL', Rsfr_full), ('sfr-no2FT', Rsfr_no2ft),
                 ('sfr-noAB1', Rsfr_noab1), ('sfr-A-only', Rsfr_aonly)):
    report(label, R)

# ---- Sub-task 1: ELEMENT-WISE doc-socc, R^matvec vs each sfrorhs variant ----
print("\n=== Sub-task 1: ELEMENT-WISE doc-socc block ===")
print("pair (socc i, doc j) | R^matvec    R^A-only   ratio(sfr/mv) | R^full")
ds = []
for n in range(nds):
    a, b, sp = prs[n]
    mv = Rmv[n]
    aonly = Rsfr_aonly[n]
    full = Rsfr_full[n]
    rat = aonly/mv if abs(mv) > 1e-12 else float('nan')
    ds.append(rat)
    print(f"  ({a},{b},{sp})  {mv:+.6e}  {aonly:+.6e}  {rat:+.4f}  | {full:+.6e}")
ds = np.array(ds)
print(f"  doc-socc ratio(A-only/mv): mean={np.nanmean(ds):+.4f} "
      f"std={np.nanstd(ds):.4f}  (clean constant => magnitude factor found)")


def elemwise(name, lo, hi, Rsfr, label):
    print(f"\n--- ELEMENT-WISE {name}  vs {label} ---")
    print("  pair (a,b,spin) |  R^matvec      R^sfr        ratio(sfr/mv)")
    rats = []
    for n in range(lo, hi):
        a, b, sp = prs[n]
        mv, sf = Rmv[n], Rsfr[n]
        rat = sf/mv if abs(mv) > 1e-10 else float('nan')
        if abs(mv) > 1e-10 or abs(sf) > 1e-10:
            rats.append(rat)
            print(f"  ({a},{b},{sp})  {mv:+.6e}  {sf:+.6e}  {rat:+.4f}")
    rats = np.array(rats)
    print(f"  {name} ratio(sfr/mv): mean={np.nanmean(rats):+.4f} std={np.nanstd(rats):.4f}")


elemwise('doc-virt', nds, min(nds+ndv, nds+10), Rsfr_aonly, 'A-only')
elemwise('soc-virt', nds+ndv, min(nconf, nds+ndv+10), Rsfr_aonly, 'A-only')

# CENTRAL TEST: FD of X^T A X should = 2FT + (hxa/hxb projection) = sfrorhs WITHOUT
# ab1 (B-matrix), since the TDA matvec has no B but DOES carry orbital-energy fa/fb
# (=> the 2FT difference-density terms appear under FD).  Compare against -Rsfr_noab1.
print("\n=== CENTRAL TEST: Rmv vs -Rsfr_noab1 (2FT on, ab1 off) ===")
for sgn, slab in ((+1, '+'), (-1, '-')):
    target = sgn*Rsfr_noab1
    c, rr = cosr(Rmv, target)
    print(f"  sign {slab}: overall cos={c:+.4f} ratio={rr:.3f}")
    for name, lo, hi in blocks:
        bc, br = cosr(Rmv[lo:hi], target[lo:hi])
        print(f"      {name:9} cos={bc:+.4f} ratio={br:.3f}")

print("\n--- doc-virt element-wise vs -Rsfr_noab1 ---")
for n in range(nds, min(nds+ndv, nds+12)):
    a, b, sp = prs[n]
    mv, sf = Rmv[n], -Rsfr_noab1[n]
    rat = sf/mv if abs(mv) > 1e-10 else float('nan')
    print(f"  ({a},{b},{sp})  mv={mv:+.6e}  -sfr={sf:+.6e}  ratio={rat:+.4f}")

# ====================================================================
#  ROOT-CAUSE: the SOMO 'ground-pair' amplitude xlr feeds an OUTPUT-side
#  1/sqrt2 fold inside mrsfesum (tdhf_mrsf_lib.F90:1496-1544). Its
#  orbital-rotation derivative is a large spurious background that
#  swamps every alpha block. The production sfrorhs applies the SAME
#  fold on the INPUT side (mrsfxvec, :1653-1667) so its derivative is
#  structurally different. We DEMONSTRATE this: zero the SOMO-row
#  amplitudes of x, re-run the FD, and watch the alpha blocks collapse.
# ====================================================================
print("\n=== ROOT-CAUSE DEMO: zero SOMO-row amplitudes -> alpha background dies ===")
Msq = X[TARGET-1].reshape((nvirb, noca)).T   # M[i, a], a->MO a+nocb
M_no = Msq.copy()
M_no[noca-2:noca, :] = 0.0       # i in {SOMO1,SOMO2} occ rows
M_no[:, 0:(noca-nocb)] = 0.0     # a in {MO noccb..noca-1} = SOMO beta-virt cols
col_no = M_no.T.reshape(-1).copy()
print(f"  |x|={np.linalg.norm(X[TARGET-1]):.4f}  |x_noSOMO|={np.linalg.norm(col_no):.4f}")
Rmv_no = Rmatvec(col_no, AXIS)
print(f"  |R^matvec(full x)|        = {np.linalg.norm(Rmv):.4e}")
print(f"  |R^matvec(x_noSOMO)|      = {np.linalg.norm(Rmv_no):.4e}")
for name, lo, hi in blocks:
    print(f"    {name:9}  |full|={np.linalg.norm(Rmv[lo:hi]):.4e}  "
          f"|noSOMO|={np.linalg.norm(Rmv_no[lo:hi]):.4e}  "
          f"(collapse factor {np.linalg.norm(Rmv[lo:hi])/(np.linalg.norm(Rmv_no[lo:hi])+1e-30):.1f}x)")
print("\n  => the alpha-block 'background' is the d/dtheta of mrsfesum's OUTPUT-side")
print("     SOMO fold (driven by xlr~0.998), NOT a clean dA/dtheta. The matvec FD")
print("     of X^T A X therefore CANNOT equal sfrorhs element-wise for any state")
print("     with significant SOMO ground-pair amplitude. doc-socc survives only")
print("     because there the fold derivative IS the dominant correct physics.")

