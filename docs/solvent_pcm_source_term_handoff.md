# PCM QM source-term fix — handoff to hermes

Pick-up spec for the next narrow coding task on the ddX/PCM line. Everything
needed to reproduce the state and make the fix is here. Read this top to bottom.

## Objective (one task only)

Fix the **QM-SCF source term** handed to ddX so the solve sees the **total**
solute source (nuclei **+** electronic density), not nuclear monopoles only.
Do **not** touch gradients, TDDFT, MRSF, excited-state PCM, nonequilibrium
solvation, or anything else.

## Accepted current state (verified this session, real ddX)

- Branch: **`feat/solvent-backend-spike`** (private repo `karmachoi/openqp-private`,
  remote `origin`). Develop on the designated working branch and fast-forward
  `feat/solvent-backend-spike`. **No PR. No public fork / upstream. No force.**
- The Fortran PCM diagnostic in `source/solvent_pcm.F90` **compiles** (gfortran 13.3).
- OpenQP **links against ddX v0.8.0**.
- **Tier-1** point-charge adapter regression **PASSES** with real ddX
  (`-0.0001789226 Eh`, within `1e-6` of the pinned value; `chg = -0.5*q_cav`
  confirmed by finite difference).
- **PCM-off H₂O vacuum is correct** (`-76.0107465`).
- **Tier-2 H₂O QM-SCF PCM is a known-blocked negative result** (do NOT treat as a pass):
  - FMM-on: ddX ddPCM Jacobi solver **does not converge** (raising `maxiter`
    100→2000 does not help).
  - FMM-off (exploratory, **not committed**): converges but unphysical —
    `e_pcm ≈ -32 Eh`, `½Tr[D·V] ≈ -16 Eh` (so `e_pcm ≈ 2·½Tr[D·V]`),
    `q_cav_sum ≈ +9.87 ≈ total nuclear charge`, `phi_cav` all-positive,
    `psi_source = nuclear_monopoles_only`.
  These are **diagnostic failure data, not reference values** — they are NOT in
  the benchmark table and must stay out of it until physically meaningful.

## Root cause (confirmed against ddX source + example)

ddX's own example (`/tmp/ddx/examples/run_ddx_low_level.c`) builds **both** the
potential and `psi` from **one consistent multipole set**:

```c
ddx_multipole_electrostatics(model, nsph, nmultipoles, solute_multipoles, ...); // phi
ddx_multipole_psi(model, nbasis, nsph, nmultipoles, solute_multipoles, psi, ...); // psi
```

and `multipole_psi` (ddX `src/ddx_multipolar_solutes.f90`) is:

```
psi(lm, j) = 4*pi / ((2l+1) * rsph(j)^l) * multipoles(lm, j)
```

i.e. `multipoles(nbasis, nsph)` are the solute's **per-sphere real-solid-harmonic
multipole moments**, and only the `l=0` (monopole) term is currently populated.

OpenQP's QM path is **internally inconsistent**:
- `phi_cav` is built **exactly** from the full QM density + nuclei
  (`source/solvent_pcm.F90`, via `int1::electrostatic_potential_unweighted` plus
  the nuclear sum), but
- `psi` is built from **nuclear monopoles only**
  (`source/solvent_ddx_adapter.c :: oqp_ddx_pcm_solve`, the
  `solute_multipoles[i] = charges[i]/sqrt(4*pi)` block).

So the forward solve sees the full solute while the adjoint/energy see only the
nuclei → `q_cav_sum ≈ Z_total`, `e_pcm` absurd. **The fix is to make the source
consistent and total.**

## Where to change (production path stays Fortran-side)

Keep the runtime path: `pcm_enabled → add_pcm_reaction_field → e_pcm →
Fortran/C ddX adapter`. Python only generates input / launches / parses.

- `source/solvent_pcm.F90 :: add_pcm_reaction_field` — owns the QM source; build
  the total per-sphere source here (nuclear + electronic) with a clearly defined
  sign convention, and pass it to the adapter.
- `source/solvent_ddx_adapter.c` — add a QM solve entry (or extend
  `oqp_ddx_pcm_solve`) that **accepts a caller-provided source** and drives both
  `ddx_multipole_electrostatics` and `ddx_multipole_psi` (or accepts a
  caller-provided `psi` array directly — `ddx_pcm_setup` takes `psi` verbatim).
  **Keep the nuclear-monopole path only for the Tier-1 point-charge smoke
  (`oqp_ddx_run_point_charge_smoke`).**
- Diagnostic line (`source/solvent_pcm.F90`): change `psi_source` so it reports
  `total_qm_multipole_source` (or equivalent) for the QM path, and keep
  `nuclear_monopoles_only` only for Tier-1. **Do not relabel until the source is
  actually changed — a false label is worse than none.**

## Two approaches — pick one (this is the open design decision)

**A. Consistent monopole net-atomic-charge source (narrow, low-risk, approximate).**
Source = net atomic charges `Z_j − pop(j)` (electronic population from
Mulliken/Löwdin; `population_analysis.F90`, needs `OQP_SM` overlap + AO→atom map).
Drive **both** phi and psi from these monopoles. Uses only the **known** `l=0`
normalization (`charge/sqrt(4π)`), so very low physics risk. `q_cav_sum → ~0`,
`e_pcm` small. **Caveat:** monopole-only drops water's dipole, so `e_pcm` will be
unphysically *small* — it satisfies the behavioral acceptance below but is not
yet physically accurate. Good as a safe first step.

**B. Full per-atom solid-harmonic multipole source (correct, larger, higher-risk).**
Compute the electronic density's **real-solid-harmonic moments per atom up to
`mmax`** (≥ dipole; ideally a few `l`), convert to ddX's normalization, add the
nuclear monopole, and drive both phi and psi via `ddx_multipole_*`. This is the
physically correct source and captures water's dipole. **Risk:** OpenQP's
existing `int1::multipole_integrals` is Cartesian, ≤ octupole, single-center —
needs per-atom + Cartesian→real-solid-harmonic conversion matching ddX exactly.

**Mandatory self-validation for B (and recommended for A):** before trusting the
energy, verify the ddX multipole-derived `phi` reproduces OpenQP's *exact*
`phi_cav` at the cavity points (within multipole truncation). If they disagree,
the normalization/convention is wrong — fix before proceeding. This is the only
internal check available (no literature reference in the build env).

Recommendation: if time-boxed, do **A** first (safe, satisfies acceptance, real
progress), then **B** for physical accuracy. Confirm the choice before coding.

## Acceptance criteria for this fix (NOT literature agreement yet)

1. Tier-1 adapter regression still passes.
2. PCM-off H₂O vacuum unchanged (`-76.0107465`).
3. H₂O PCM-on no longer shows nuclear-only behavior.
4. `q_cav_sum` is **no longer ≈ total nuclear charge** (should be ≈ net molecular
   charge, ~0 for neutral water).
5. `e_pcm` is **no longer tens of Hartree**.
6. Diagnostics clearly report the source convention (`psi_source`).
7. Re-run H₂O with **committed/default FMM-on first**. Do **not** "fix" by
   committing FMM-off; use FMM-off only as a diagnostic comparison if FMM-on
   still fails. (FMM-on non-convergence for tiny molecules is a *separate* open
   item — see below.)
8. **Do not populate the H₂O benchmark table** unless the value is physically
   meaningful AND protocol-matched.

## Deferred to AFTER the source-term fix (do not start these now)

`e_pcm` vs `½Tr[D·V]` bookkeeping (measured factor ~2 + scale error); `q_cav`
sign/scale; `phi_cav` convention; **FMM-on convergence for small molecules**
(FMM-off converges, FMM-on does not — likely fmm_lmax/cavity related);
literature/trusted-reference H₂O comparison.

## Reproducible ddX-enabled build recipe (this worked, from scratch)

```bash
# toolchain (root; apt main archive — ignore the deadsnakes/ondrej PPA 403s,
# or delete /etc/apt/sources.list.d/{deadsnakes,ondrej}*.list first)
apt-get install -y --no-install-recommends gfortran libopenblas-dev liblapack-dev
python -m pip install numpy cffi scipy basis_set_exchange libdlfind   # pyoqp runtime deps (the independent ddPCM reference NOT needed to run)

# ddX v0.8.0
git clone --depth 1 https://github.com/ddsolvation/ddX.git /tmp/ddx
cmake -S /tmp/ddx -B /tmp/ddx-build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DTESTS=OFF -DEXAMPLES=OFF -DWARNING_FLAGS=OFF -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_CXX_COMPILER=g++
cmake --build /tmp/ddx-build -j4
mkdir -p /tmp/ddx-install/include /tmp/ddx-install/lib
cp /tmp/ddx/src/ddx.h /tmp/ddx-install/include/
cp /tmp/ddx-build/src/libddx.so /tmp/ddx-install/lib/

# OpenQP with ddX (internal ERIs: USE_LIBINT=OFF; 32-bit BLAS int)
cmake -B /tmp/oqp-build -G Ninja -DUSE_LIBINT=OFF \
  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_INSTALL_PREFIX=$PWD -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=OFF \
  -DENABLE_DDX=ON -DDDX_ROOT=/tmp/ddx-install
cmake --build /tmp/oqp-build --target oqp -j4     # liboqp.so (~75 MB); confirms solvent_pcm.F90 compiles
cmake --build /tmp/oqp-build --target install -j4 # populates lib/, share/basis_sets/

# Tier-1 adapter smoke standalone (only needs ddX + 2 C files)
gcc -O2 -DOQP_ENABLE_DDX -I source -I /tmp/ddx-install/include \
  tests/ddx_adapter_smoke.c source/solvent_ddx_adapter.c \
  -L /tmp/ddx-install/lib -lddx -lgfortran -lm -fopenmp -o /tmp/adapter_smoke
LD_LIBRARY_PATH=/tmp/ddx-install/lib OMP_NUM_THREADS=1 /tmp/adapter_smoke
```

### Running the H₂O QM PCM path / the tests

`pytest` must run on the interpreter that has the deps (`/usr/local/bin/python`
here; install pytest into it). The OpenQP subprocess uses `sys.executable`.

```bash
OPENQP_ROOT=$PWD PYTHONPATH=$PWD/pyoqp LD_LIBRARY_PATH=/tmp/ddx-install/lib \
  OMP_NUM_THREADS=1 OQP_DDX_ADAPTER_SMOKE=/tmp/adapter_smoke \
  /usr/local/bin/python -m pytest -rsA -q tests/test_pcm_literature_benchmarks.py
```

Expected today: H₂O vacuum + Tier-1 + parser PASS; Tier-2 diagnostics SKIP with
the KNOWN-BLOCKER message; references SKIP pending. After the fix, the Tier-2
diagnostics test should run (FMM-on permitting) and show the corrected
`q_cav_sum`/`e_pcm`/`psi_source`.

The `PCM diag` lines parsed by the harness come from
`source/solvent_pcm.F90 :: add_pcm_reaction_field`. Reference table:
`tests/data/pcm_literature_benchmarks.json`; doc: the empirical section in
`docs/solvent_pcm_literature_benchmarks.md`.

## Guardrails (standing)

- Production PCM physics stays in Fortran/C/ddX; Python only configures/parses.
- No fabricated numbers. No benchmark reference populated unless physically
  meaningful and protocol-matched.
- Keep all open limitations visible in diagnostics/doc.
- No PR; push only to the private branch after the targeted tests pass.
- Do not add gradients/TDDFT/MRSF/excited-state/nonequilibrium/architecture.
