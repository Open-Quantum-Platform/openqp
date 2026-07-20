# Phase 11 — Frozen-Fock FD vs production sfrorhs RHS (findings)

Goal: make `R^matvec` (FD of `X^T A X` w.r.t. ROHF orbital rotation, using the
standalone frozen-Fock TDA matvec) match the production z-vector RHS `R^sfrorhs`
on the DIAGONAL (state 2). By Hellmann-Feynman they were expected to agree for the
converged eigenvector.

**HEADLINE RESULT (rigorous):** The diagonal match is NOT achievable with the
FD-of-`X^T A X` approach for state 2 (or any MRSF state whose dominant character is
the open-shell SOMO ground-pair). The obstacle is precisely located: the matvec's
`mrsfesum` applies the open-shell `1/sqrt2` SOMO fold on the OUTPUT side, whose
orbital-rotation derivative — driven by the huge ground-pair amplitude
`xlr ~= 0.998` — is a large spurious background that swamps every alpha rotation
block. The production `sfrorhs`/z-vector path applies the SAME fold on the INPUT
side (via `mrsfxvec`), so its derivative is structurally different. The two are
NOT the same orbital functional away from theta=0, even though both equal Omega at
theta=0. This is a hard structural fact, not a tunable convention bug.

--------------------------------------------------------------------------------
## Run env (verbatim)
```
cd /tmp/nactest
export LD_PRELOAD=/opt/soft/install/GCCcore/12.3.0/lib64/libstdc++.so.6
export LD_LIBRARY_PATH=/opt/soft/install/GCCcore/12.3.0/lib64:/bighome/alireza/.local/lib
export OPENQP_ROOT=/bighome/alireza/openqp-nac
export PYTHONPATH=/bighome/alireza/openqp-nac/pyoqp
export OMP_NUM_THREADS=4
unset NAC_ZERO_2FT NAC_DUMP_RHS
timeout 400 /opt/conda/bin/python3 -u .../matvec_rhs.py > out.txt 2>&1
```
System: H2O, 6-31g*, bhhlyp, MRSF triplet ROHF ref, singlet targets (mrst=1).
nbf=19, noca=6, nocb=4, nvirb=15, nij=90; nconf(rotation)=86. TARGET=2.
Omega(Hartree) = [-0.283490, 0.045304, 0.103566].
State-2 amplitude: M[i=4(lr1), a=MO4] = xlr = -0.99847 (pure SOMO1->SOMO1 ground
pair). All three states put ~99.8% of |X| in the SOMO rows.

--------------------------------------------------------------------------------
## CONVENTIONS RESOLVED (these ARE correct / confirmed)

1. AXIS: 'col' (MOs are Python COLUMNS of OQP::VEC_MO_A/B). 'row' axis is wrong
   (overall cos -0.32 vs +0.37). Confirmed; harness now uses AXIS='col' only.

2. Amplitude index layout (iatogen/mrsfesum): `av(1:nocca, noccb+1:nbf)`,
   column-major with i (alpha-occ, 1..nocca incl. 2 SOMOs) FASTEST, j (beta-virt,
   noccb+1..nbf incl. 2 SOMOs) SLOWEST. In Python: `M = Ax.reshape((nvirb,nocca)).T`
   gives M[i,a] with MO(a)=a+noccb. Verified (Step A: X.AX=Omega to 1e-16 states 2,3).

3. Givens rotation generator vs sfrogen: for a pair (a=high, b=low),
   `givens(C,a,b,th)` does C[:,a]=c*C[:,a]+s*C[:,b]; C[:,b]=-s*C[:,a]+c*C[:,b].
   This is the generator kappa(low=b, high=a)=+th, matching sfrogen's
   ava(low,high)=pv / avb(low,high)=pv layout. SIGN AND INDEX ARE CONSISTENT.
   (Swapping a<->b or th->-th only flips the global sign and degrades cos, so the
   current convention is the right one.)

4. Spin-resolved blocks (sfrogen tdhf_sf_lib.F90:504-528) — CONFIRMED:
   doc-socc rotates BETA only (avb(j,i)); doc-virt rotates BOTH (ava=avb(j,k));
   socc-virt rotates ALPHA only (ava(i,k)). prs[] order in the harness matches.

--------------------------------------------------------------------------------
## SUB-TASK 4 RESOLVED FIRST (it reframes everything): sfrorhs is NOT the bare
## FD of X^T A_TDA X.  (source/tdhf_sf_lib.F90:345-434, z_vector.F90:1640-1775)

`sfrorhs` returns  rhs = -( hpt + xhx_blocks )  where, per block (i=socc, j=doc,
k=virt):
  doc-socc : rhs = -[ hptb(j,i) + hxa(i,j) - hxa(j,i) - hxb(j,i) + 2FT ]
  doc-virt : rhs = -[ hpta(j,k) + hptb(j,k) + hxa(k,j) - hxb(j,k) + 2FT ]
  socc-virt: rhs = -[ hpta(i,k) + hxa(k,i) + hxb(k,i) - hxb(i,k) ]
Pieces:
  - hxa/hxb (= hxa,hxb at z_vector.F90:1707-1728) are the antisymmetric projection
    of the FULL-MO A.x kernel G_MO = mo_a^T * agdlr^AO * mo_a, materialized as
    hxa = 2*G_MO*X^T (alpha) and hxb = 2*G_MO^T*X (beta). i.e. the 2-ELECTRON
    part of A.x, projected to the rotation space as [G(p,q)-G(q,p)] with the spin
    channel fixed by sfrogen.
  - 2FT = the 2*Fa*Tij + 2*Fb*Tab ORBITAL-RELAXATION terms (sfrorhs:384-403),
    built from the difference density T (Tij=-X X^T occ-occ, Tab=+X^T X vir-vir).
    These carry the orbital-ENERGY dependence (the analog of the matvec's fa/fb
    diagonal). Gated by NAC_ZERO_2FT.
  - hpta/hptb (= ab1_mo_a/b, z_vector.F90:1668-1671) are the (A+B) / B-matrix
    response density via plain mntoia. ABSENT from a pure TDA A.x matvec.
    Gated by NAC_ZERO_AB1.

CONSEQUENCE: R^matvec (FD of X^T A_TDA X) and R^sfrorhs differ by KNOWN pieces.
A clean comparison would be: FD-of-(X^T A_TDA X)  vs  -[ hxa/hxb projection + 2FT ]
= sfrorhs with AB1 gated off. (AB1 is genuinely absent from the TDA matvec.)
The √2 SOMO fold is injected on the INPUT (X via mrsfxvec, :1653-1667) in the
production path, but on the OUTPUT (wrk1 fold) in the matvec's mrsfesum
(:1496-1544) — this placement difference is the root cause (below).

--------------------------------------------------------------------------------
## SUB-TASK 1: ELEMENT-WISE doc-socc (where cos was already +0.977)

prs doc-socc order: (4,0),(4,1),(4,2),(4,3),(5,0),(5,1),(5,2),(5,3)
[i=4=lr1(SOMO1), i=5=lr2(SOMO2); j=0..3 doc]. Compared vs sfrorhs A-only:
  (4,0,b)  mv=-3.96e-04  sfr=~0          <- ZERO-PATTERN MISMATCH (lr1 rows)
  (4,1,b)  mv=+3.41e-04  sfr=~0          <- mv nonzero, sfr zero
  (5,0,b)  mv=-7.11e-03  sfr=-9.56e-03   ratio 1.34
  (5,1,b)  mv=+1.38e-02  sfr=+1.04e-02   ratio 0.75
  (5,3,b)  mv=~0         sfr=-3.50e-03   <- sfr nonzero, mv zero
NO clean constant ratio (mean 0.52, std 0.56). The match is NOT a magnitude
factor: it is a zero-pattern + per-element scramble localized on the SOMO rows
(lr1 vs lr2). This is the input-side (mrsfxvec) vs output-side (mrsfesum fold)
placement of the √2 SOMO normalization — exactly where the design analysis
predicted the discrepancy would live.

The doc-socc cos=+0.977 in the FULL comparison is dominated by the two large
lr2-row elements (5,0),(5,1); it is NOT an element-wise match.

--------------------------------------------------------------------------------
## SUB-TASK 3: ALPHA blocks (doc-virt both-spin, socc-virt alpha-only)

doc-virt matvec is dominated by a spurious ~+1.30e-2 background:
  (6,2,ab)=(6,3,ab)=(7,2,ab)=(7,3,ab) ~= +1.303e-2   (sfrorhs A-only ~ 1e-3..0)
Decomposition (single-pair probe): (6,2,ab)= (6,2,a)=+4.92e-2 + (6,2,b)=-3.62e-2.
The alpha piece +4.92e-2 is the artifact. socc-virt is worse (ratios -250..+2.6).

--------------------------------------------------------------------------------
## ROOT CAUSE (DECISIVE, sub-task 3 + 4): output-side SOMO fold derivative

Probe: zero the SOMO-row amplitudes of X (set M[i in {lr1,lr2}, :] = 0 and
M[:, a in {SOMO beta-virt cols}] = 0), then re-run the FD. Result:
  |R^matvec(full x)|   = 1.234e-01
  |R^matvec(x_noSOMO)| = 6.26e-04         (197x smaller)
  per block collapse:  doc-socc 545x,  doc-virt 78x,  socc-virt 277x.
So ~99% of the FD signal in EVERY block is the orbital-rotation derivative of the
SOMO fold, driven by xlr~0.998 (state 2 is the pure open-shell-singlet ground
pair). This is mrsfesum's OUTPUT-side fold (:1496-1513): terms like
fab(j,lr)*xlr*√2, fij(i,lr)*xlr*√2, and the wrk1(lr1,lr1) energy term
xlr*(fab+fab-fij-fij)*0.5 — all proportional to xlr and all functions of the
rotating fa/fb. Their d/dtheta is a real number but it is NOT a clean dA/dtheta;
it is an artifact of holding xlr fixed while the fold's Fock factors rotate.

Because the production path folds on the INPUT (mrsfxvec rescales the amplitude
itself by ±√2 BEFORE building hxa/hxb), its derivative w.r.t. theta does NOT
generate this xlr*d(fa/fb)/dtheta background in the same form. The two evaluators
agree at theta=0 (X.AX=Omega to 1e-16) but their orbital-gradients diverge on the
SOMO sector. doc-socc partially survives because there the fold derivative happens
to BE the dominant correct physics; the alpha blocks (doc-virt, socc-virt) are
pure artifact.

--------------------------------------------------------------------------------
## DIAGONAL MATCH: NOT ACHIEVED. Final per-block cos/ratio (state 2, AXIS=col):

vs sfrorhs FULL:     overall cos -0.022;  doc-socc +0.977/0.303,
                     doc-virt -0.244/1.28, soc-virt -0.144/2.71
vs sfrorhs A-only:   overall cos -0.042;  doc-socc +0.935/1.071,
                     doc-virt +0.177/15.5, soc-virt -0.054/0.534

Only doc-socc aligns in direction, and even there it is a 2-element coincidence,
not an element-wise match (ratio scatter mean 0.52 std 0.56, zero-pattern wrong on
lr1 rows). The alpha blocks do not align under any sign/factor.

--------------------------------------------------------------------------------
## PRECISE REMAINING OBSTACLE & RECOMMENDED NEXT STEP

OBSTACLE: FD of `X^T A_TDA X` using the standalone matvec is the WRONG operator for
the rotation gradient on the MRSF SOMO sector, because mrsfesum's √2 SOMO fold is
on the output. To get a faithful `dOmega/dtheta` you must differentiate the SAME
functional the production code does, i.e. fold on the INPUT.

FIX ATTEMPT THAT FAILS (recorded so it is not retried): feeding the matvec the
mrsfxvec-folded amplitude (xlr -> ±xlr/√2 in the two SOMO diagonal slots) does NOT
fix it — it DOUBLE-folds, because mrsfesum STILL applies its own output fold on top.
Measured: doc-virt cos goes from +0.18 to -0.42 (worse), soc-virt ~0. So the input
refold alone is wrong.

CORRECT FIX (requires exporting the full-MO A.x kernel — a Fortran add, hence OUT
OF SCOPE for this Python-only task):
  The right comparison object is the antisymmetric projection [G(p,q)-G(q,p)] of
  the matvec's OWN full-MO A.x matrix G_MO = mo_a^T * agdlr^AO * mo_a, contracted
  against the mrsfxvec-folded X — i.e. exactly hxa = 2*G_MO*X_folded^T and
  hxb = 2*G_MO^T*X_folded. This reproduces sfrorhs's hxa/hxb (A-only) by
  construction, with the √2 fold on the INPUT only (no FD, no output fold). The
  standalone matvec presently exports only the (occ_a,virt_b) amplitude block
  amo(:,1)=A.x, NOT the full nbf x nbf G_MO/wrk2. To take this route, add ONE
  tagarray export of wrk2=G_MO inside mrsf_matvec_apply (3 lines), then in Python:
  form hxa/hxb from G_MO and X_folded, build the 3 antisymmetric block projections
  per the sub-task-4 table, add the 2FT term (also computable from fa/fb + Tij/Tab,
  both already exportable), and compare to sfrorhs noAB1. AB1 (B-matrix) is the only
  genuinely-absent piece (TDA matvec has no B), so the closest achievable target is
  sfrorhs with NAC_ZERO_AB1.

DO NOT spend more effort on sign/factor tweaks of the current FD — the root cause is
structural (output- vs input-side SOMO fold), proven by the 78-545x amplitude
collapse and the failed double-fold experiment.

State 1 / state 3 are NOT cleaner test cases: they also carry ~99.8% SOMO-row
amplitude (SOMO1->SOMO2 and SOMO->virt), so the same fold contamination applies;
state 2 is merely the extreme (pure ground pair).

================================================================================
## 2026-06-17 — ORACLE CLOSED. The analytic amplitude term = the matvec-A bilinear.
================================================================================

RESULT (p11_poc_gaugefree.py): the semi-numerical amplitude term of the working
analytical NAC (benchmark_full_nac.py d_amp = X_I^T dX_J, the ONE non-analytic
piece) is reproduced by the matvec operator form
    d_amp(I,J) = X0_I^T (dA/dR) X0_J / (Omega_J - Omega_I)
for ALL THREE pairs to the FD floor:
    (1,2) cos=+1.00000000 ratio=1.000000 resid 1.9e-7 (0.000%)
    (1,3) cos=+1.00000000 ratio=1.000001 resid 3.7e-8 (0.000%)   [headline |d_amp|=0.0284]
    (2,3) cos=+1.00000000 ratio=1.000023 resid 1.5e-5 (0.003%)
I=J diagonal self-test (gmo_validate.py) STILL cos=+1.0 / 1e-15. No Fortran edited.

KEY CORRECTIONS to the prior synthesis (which was BLOCKED on an antisymmetric-sp lead):
 1. The test_QZ "deficiency" (D=coded-numeric, cos -0.9997 etc.) was a RED HERRING:
    test_QZ's hand-rolled numeric reference (sum_K Omega_K <Z|X_K>^2 with SVD MO
    alignment) is NOT the validated oracle. The REAL oracle is the benchmark's
    semi-numerical d_amp = X0_I^T dX_J (TLF FD of displaced eigenvector amps), which
    drives the cos=1.0 total NAC. test_QZ's symmetric polarization Gfull(Z) is a
    DIFFERENT object (gradient of the quadratic form, not X_I^T dX_J).
 2. The matvec A in the full 90-dim amplitude space is SYMMETRIC to |A-A^T|=2.8e-13
    (eigvalsh gives exactly Om). So the "input-vs-output SOMO fold makes A
    non-symmetric / antisymmetric-sp is the missing physics" hypothesis is WRONG:
    A is symmetric, the PT identity X_I^T dX_J = X_I^T dA X_J/(Omega_J-Omega_I) is
    exact, and the genuinely-new object is NOT antisymmetric.
 3. The whole rotation-space / seam machinery was the wrong arena. The amplitude
    term is a NUCLEAR derivative of A; the clean evaluation is gauge-free.

THE TWO THINGS THAT ACTUALLY CLOSED IT (both pure gauge fixes in the FD harness,
NOT physics — the operator X_I^T dA X_J/gap was correct all along):
 (a) GAUGE-FREE A-transport: build the displaced matvec A column-by-column (A is
     linear at int2e_cutoff=1e-20), transport it into the FIXED reference MO frame
     A_ref = T^T A_disp T (T = amplitude-space rep of the per-block orbital
     transport), FD A_ref, contract with FIXED X0. This removes the quadratic
     gauge-contamination of the naive X_I(R)^T A(R) X_J(R) bilinear FD (which only
     ever closed 2 of 3 pairs in any per-state-sign gauge).
 (b) The orbital transport block Q must be SIGN-CONTINUOUS: use LOEWDIN symmetric
     orthogonalization Q = sub (sub^T sub)^{-1/2} of the ref x displaced MO-overlap
     sub-block, NOT the SVD Procrustes Vt^T W^T (whose per-direction sign flip
     corrupted (1,2)/(1,3)). Loewdin alone took (1,3),(2,3) to cos=1.0 exact and
     (1,2) to cos=1.0 ratio=0.7090 (=1/sqrt2).
 (c) The remaining (1,2) 1/sqrt2: the spin-adapted GROUND-CONFIG slot ijlr1 does
     NOT transport like a plain grid entry. transport_T must UNFOLD ijlr1 into the
     two SOMO-SOMO determinant slots (+/-1/sqrt2), transport in the DETERMINANT
     grid, then REFOLD (sqrt2 back on ijlr1, ijlr2 slot empty) — exactly the
     benchmark's displaced_amps fold handling. This restored the sqrt2 ONLY for
     the SOMO-sector pair (1,2) (state1=SOMO2->SOMO1 idx5, state2=ground ijlr1 idx4,
     both in the folded SOMO block); (1,3)/(2,3) involve state3 (SOMO->virt, idx17)
     outside the fold and were unaffected. Final: ALL THREE cos=1.0 ratio=1.0.

WHY (1,2) needed sqrt2 but (2,3) did not (both contain the ground state 2): the
fold couples states whose dominant configs BOTH live in the SOMO-SOMO block.
state1(idx5 SOMO2->SOMO1) and state2(idx4 ground SOMO1->SOMO1) are both in-block ->
sqrt2. state3(idx17 SOMO->virt) is out-of-block -> no fold factor.

FORTRAN PORT (mechanical, the milestone is the operator validation above):
 - The analytic dA/dR is the matvec A's nuclear derivative = the standard CPHF
   orbital-response (U^x) of the frozen-Fock matvec, which the project's validated
   analytic-gradient machinery already computes. Replace the displaced-SCF FD by
   the U^x response: dA = d(orbital)A via the existing z-vector/CPHF, contracted
   with the exported kernels (nac_gmo/nac_gchan/nac_fa/nac_fb). The SOMO fold
   (unfold->transport->refold) maps to applying mrsfxvec on the ground-config slot
   in the interstate contraction — i.e. the bilinear must fold the SOMO-SOMO sector
   exactly as mrsfxvec does on input. No new export needed; the kernels suffice.
 - Validation gate for the port: reproduce p11_poc_gaugefree.py's three cos=1.0
   ratio=1.0 against p11_damp_oracle.npz, keep gmo_validate I=J at 1e-15, keep the
   gradient anchor O dE/dZ=-0.182993575.

Scripts: p11_poc_gaugefree.py (THE closing PoC), p11_poc.py (per-state-sign gauge,
shows 2/3 + the gauge diagnosis), p11_damp_analytic.py (naive bilinear FD, the
gauge-sensitivity demo), p11_grad_decomp.py / p11_ab1_nuclear.py (nuclear RHS-piece
decompositions that ruled out single-term fixes). Oracle: p11_damp_oracle.npz
(= benchmark_bhh.npz d_amp). A matrix: p11_Amat.npy (symmetric, eigs==Om).
