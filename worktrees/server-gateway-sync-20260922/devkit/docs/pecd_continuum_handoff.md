# PECD development handoff (openqp-spec)

Status: engine complete and validated (M1-M5), manuscript v0-3 pushed to
Overleaf (git.overleaf.com/6a465178be75058390d95b4f, commit 2374599).
This file lists the open follow-ups with enough context to resume cold.

## Where things live
- Engine: `src/` in THIS repo (fetched into the OpenQP build when `-DENABLE_SPEC=ON`; `[pecd] engine=1`)
- Standalone validation tests: `tests/` in THIS repo
  (M2_SPEC.md / M3_SPEC.md = pinned conventions; test_m2/test_m3g*/test_m4g* = the test suite)
- Methyloxirane results + literature comparison: `tests/moxirane_validation/`
  (RESULTS_DOSSIER.md = numbers source-of-truth; LITERATURE_GATES.md = verified citations)
- Python API: `Molecule.get_pecd_results()`; example `examples/pecd_pythonic.py`
- Docs: openqp-docs branch `pecd-docs` (keywords/pecd.md + workflows/pecd.md)
- Run inputs: session dir `run_moxirane/bs_{R,S}_l*.inp`, `lumo_{R,S}.inp`

## Follow-ups

### 1. Velocity gauge (pre-submission)
Length gauge only today. Add gradient dipole matrix elements in
`src/pecd_bs_m4.f90` (I_raw variant with p operator), verify
length/velocity agreement on H2+ (exact system: L/V equality is a real test
there; target <0.5%), then report the L/V spread per energy for production
molecules as a model-quality indicator (NOT a convention arbiter).

### 2. Excited-state spectra with engine=1 (pre-submission)
- Valence S1 channels (pole 0.76-0.97) should run as-is: `[tdhf] target=2`,
  pick iroot from the printed EKT table.
- BLOCKED: LUMO-electron channel (S1 iroot=17, pole 0.12, eBE 0.41 eV)
  gave Dyson norm 44 (unphysical). Needs (a) normalization handling of
  weak-pole EKT Dyson vectors in `build_dyson_orbital`/projection, and
  (b) diffuse basis (aug-cc-pVDZ) for the quasi-Rydberg orbital.
  This is the TR-PECD probe channel - worth fixing properly.

### 3. LB94-class potential upgrade (revision stage)
Replace/augment Xalpha in `src/pecd_bs_potential.f90` (pointwise vxc
on the angular grid is the only touch point). Mostly shifts b1 node
positions (phases).

### 4. Absolute helicity<->LCP sign pin (revision stage)
Use the Tia et al. JPCL 8, 2780 (2017) O1s fixed-in-space measurement as
the anchor (manuscript currently hedges this correctly, Sec. results item iv).

### 5. Manuscript logistics
Co-authors/affiliations, final proofread, figure polish. v0-3 = current;
v0-2/v0-1 kept for history. Style rules: no "honest", no CS jargon
(gate/fingerprint/pipeline/falsifiability), validation in SI except the
single H2+ table, main text = physics of the system.

### 6. openqp-docs PR
Branch `pecd-docs` on Open-Quantum-Platform/openqp-docs is ready; open the
PR when this feature branch merges upstream.

### 7. BLAS platform policy PR (separate thread)
Deterministic platform->BLAS map (arm mac Accelerate / intel MKL / arm linux
OpenBLAS), wheel/CI alignment; see memory/branch from the "Harden OpenBLAS"
session.

### 8. liboqp_spec as a literal static archive (cosmetic)
`spec/` + ENABLE_SPEC exist; converting to a standalone .a is an
object-library split in source/CMakeLists.txt (compile-property migration).

### 9. G9 robustness sweep on methyloxirane (cheap insurance)
h -> h/2, rmax +20, matching-window shift, angular-grid doubling at
production settings; lmax scan + box extension already done.

### 10. TR-PECD (follow-up paper)
b1(E, tau) maps over NAMD snapshots; needs item 2 first. Python-prototype
driver logic exists in the session tree (`mrsf_pecd/trpecd.py`).

## Performance of the multi-centre / B-spline path (2026-09-20/21)

Branch `agent/pecd-mc-analytic-coupling-20260920`.  N2 3sigma_g, lmax 10, MC=1,
Lpot=20, four photon energies, six threads, Ultra:

| phase | before | after |
|---|---|---|
| multi-centre setup (mc_setup + mc_add_short_range + mc_build) | ~140 s | 3.9 s |
| vpot_build | 60.8 s | 9.2 s |
| project_dyson | 26.6 s | 5.6 s |
| cc_solve | ~70 s | ~70 s (LAPACK-bound) |
| **total** | **290 s** | **90 s** |

All four partial cross sections are BIT-IDENTICAL across every step
(5.128879 / 10.443763 / 3.411166 / 2.241327 Mb).  Thymine root 33 also
reproduces the old numbers to 2 decimals, with lmax 12 falling 5226 s -> 1168 s.

**The single recurring defect**: every hot spot was a large contraction written
either as a rank-1 update in the innermost loop, or through the Fortran MATMUL
intrinsic, which does not reach the vendor BLAS for these shapes (measured
~0.5 GFLOP/s against Accelerate's tens of GFLOP/s).  Five separate sites had it:
the atom/single-centre ring in `mc_build`, the short-range channel coupling in
`mc_add_short_range` (121 x 8256 x 121 per radial node, ~1100 nodes), the three
Schur updates in `mc_schur` (pecd_mc_types.f90), the multipole projection in
`vpot_build`, and the harmonic projection in `project_dyson`.  A new
`oqp_dgemm_i64` wrapper was added to the host `source/mathlib/lapack_wrap.F90`
(branch `agent/pecd-mc-analytic-host-clean-20260920` on GitLab `cheol/openqp`).

If you profile this code again: **instrument first**.  Three successive
structural guesses about the bottleneck were all wrong, and one "fix" (a cache
-locality loop interchange) made it slower.  `system_clock` around each phase
found it immediately.

`mc_build` also gained an analytic-azimuth path for centres on the z axis: with
the centre on z the global azimuth equals the local one and r, cos(theta_g) are
azimuth-independent, so every real harmonic factorises as A(cos theta)*Phi_mu(phi)
and the azimuthal quadrature collapses onto two small matrices, G for the
overlap/kinetic coupling and T for the potential, both formed on the SAME uniform
phi grid so the reduction is exact rather than an approximation.  Off-axis centres
still use the grid path; generalising needs the real-harmonic rotation matrices
(direct recursion, Ivanic-Ruedenberg / Choi).  Verified bit-identical on N2, where
both centres sit at (0,0,+-1.037).

## Correction (2026-09-21): the captured-norm anomaly was a Dyson-amplitude bug

An earlier version of this file read the captured Dyson norm exceeding 1 as an
aliasing signature in the projection quadrature.  That was wrong.  `tdhf_mrsf_ekt`
stores the metric-normalised eigenvector x in `OQP_mrsf_ekt_orbitals_mo`
(`cdys = U_keep x`) and computes the pole strengths separately; `pecd_dyson`
back-transformed that column straight to AO, so the engine ionised x instead of
the Feynman-Dyson amplitude d = P x.  Applying the stored EKT metric density
first fixes it:

| thymine, captured Dyson norm | before | after |
|---|---|---|
| lmax 12 | 1.0143 | 0.9542 |
| lmax 16 | 1.0339 | 0.9720 |
| lmax 20 | 1.0422 | 0.9797 |

The norm now stays below 1 and rises towards it with lmax, which is the correct
behaviour.  N2 partial cross sections move by only 0.25% because those channels
have pole strength ~1, where P x and x nearly coincide.

**The convergence drift is a separate problem and got worse, not better.**
At 10 eV the lmax 16 -> 20 change went from -27.2% to -37.5%.  So the captured-norm
anomaly and the cross-section drift had different causes; do not treat them as one.

## Open: thymine lmax convergence -- BOTH continuum methods (NOT solved)

Neither the B-spline molecular continuum nor ezDyson has a converged thymine
partial cross section.  They fail differently and the ezDyson case is the
cheaper one to close.

### ezDyson: convergence test never completed (blocked, then unblocked, not redone)

Thymine root 30, single-centre ezDyson, Sigma_tot (Mb):

| l_max | 21.2 eV | 100 eV |
|---|---|---|
| 3  | 0.327 | 0.0121 |
| 10 | 0.994 | 0.2307 |
| 14 | died  | died |
| 18 | died  | died |

A factor of three between l=3 and l=10 with no third point, so nothing can be
said about convergence.  l=14 and l=18 aborted with `Invalid lmax` from
`ezdyson_code/klmgrid.C:80` (`check(!(lmax<theSPH().LMax()),"Invalid lmax")`):
stock ezDyson 2021 carries hard-coded spherical-harmonic and Clebsch-Gordan
tables that stop at l=10.

That limit was REMOVED in the local build at
`~/ezdyson-build-20260916/ezDyson_2021` (general-l formulas replacing the
tables in `sph.C` and `clebsh_gordan_coeff.C`; agreement with the original
1.8e-13 within l<=10, and all 8 shipped sample inputs reproduced to six
digits).  **The l>10 thymine scan was never rerun with the patched binary.**
That is the first thing to do here -- it is minutes of compute, and it decides
whether ezDyson's thymine number means anything.

Separately, ezDyson's effective charge stays a free parameter (Gozem 2015).
At Z=1 it misses the measured H2O/N2 channels by factors of 2 to 8, so it is a
comparison baseline in the manuscript, not a parameter-free prediction.

### B-spline molecular continuum: converges for small molecules, not for thymine

Thymine root 33, multi-centre, partial cross sections (Mb):

| lmax | 6 eV | 10 eV | 14 eV | 20 eV | 30 eV |
|---|---|---|---|---|---|
| 12 | 21.46 | 23.05 | 12.05 | 6.85 | 3.63 |
| 16 | 18.18 | 16.60 | 9.29 | 5.80 | 3.10 |
| 20 | 17.99 | 12.08 | 7.78 | 5.39 | 2.88 |

Still falling ~27% from lmax 16 to 20 at 10 eV, with no sign of flattening.
What is ruled out:
- NOT angular convergence of the continuum: at 10 eV the per-l cross section
  (PECD-SIGL) is below 1e-3 Mb beyond l=10, yet the l<=6 contributions themselves
  halve as lmax grows from 12 to 20.  Adding channels that carry no intensity
  changes the converged ones.
- NOT the Dyson projection order: freezing it (`OQP_PECD_LDYSON=11`) while raising
  the continuum lmax leaves the drift, and flips its sign at high energy
  (+15.5% at 20 eV, +19.7% at 30 eV from lmax 16 to 20).
Still suspect: conditioning of the anchored banded solve (rcond ~1e-23 even after
the symmetric diagonal equilibration), and the Dyson projection quadrature -- the
captured Dyson norm EXCEEDS 1 and grows with lmax (1.014 at l12, 1.034 at l16,
1.042 at l20) where N2 gives 0.9993.  A projection onto orthonormal channels
cannot exceed the true norm, so that is an aliasing signature on the angular grid
(`nth = 2*lmx+32`, `nph = 4*lmx+65` in project_dyson) for a molecule with 12
off-centre nuclei.  Small molecules are unaffected.

NOTE: no number in the JCTC manuscript depends on this.  The thymine figures use
the Gelius model, which has no lmax; the B-spline continuum appears only for H2O
and N2, which are converged (H2O 1b1 lmax 16 and 20 agree to three digits).

## Gotchas (hard-won)
- NEVER `cp` over the installed (memory-mapped) liboqp.dylib: `rm` first,
  then copy (and `codesign -f -s -` if needed), else dyld SIGKILLs every
  new loader (exit 137).
- basis%cc is already axial-normalized: primitive builders must apply only
  the angular ratio gto_norm(lx,ly,lz,1)/gto_norm(L,0,0,1).
- Tight Gaussians (zeta>25) off-centre alias any angular grid: keep the
  analytic point-charge / displaced-Gaussian split (see src/pecd_bspline.F90).
- b1 near a node of b1(E) is maximally phase-sensitive: quote node
  positions, not point values (methyloxirane 4 eV).
- Pythonic Runner(input_dict=...): all values strings; inline geometry
  starts with a blank line; usempi=False for serial scripting.
