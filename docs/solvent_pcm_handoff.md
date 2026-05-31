# PCM work — session handoff

Entry point for whoever picks up the PCM/ddX solvent work next. Read this
first, then the two detail docs it points to.

## TL;DR

The first closed RHF/DFT PCM single-point **energy** path is wired, builds, and
is committed on the private branch. It is **not physically validated**, and there
are now **two coexisting PCM implementations** that must be reconciled. The next
activity is validation + reconciliation, **not more architecture**.

## Where the work lives

- Repo / remote: **`openqp-private`** (`git@github.com:karmachoi/openqp-private.git`),
  remote name `private`. **Do not push to `origin` (the fork) or `upstream`.**
- Branch: **`feat/solvent-backend-spike`**, tip **`38f2186`**.
- Key commits:
  - `cdddfba` — `--no-ff` merge integrating the local ddX/Fortran PCM path onto
    the private solvent line (parents: private `1f5b173` + local PCM `9d3026e`).
  - `9f0777d` — ddX-gated PCM energy-path wiring.
  - `9d3026e` — validation-gate doc + provisional-convention inventory.
- Local worktree: `openqp-pcm` (kept intact). Safety refs:
  `backup/local-pcm-before-private-merge` (`9d3026e`),
  `backup/private-solvent-backend-spike` (`1f5b173`).

## Build & run (local gotchas)

- Configure with **`-DLINALG_LIB_INT64=OFF`** (Accelerate rejects ILP64) and
  **`-DCMAKE_CXX_COMPILER=g++-15`**; toolchain gcc-15/gfortran/ninja, `LINALG_LIB=auto`.
- `ninja -C build oqp` outputs `build/source/liboqp.dylib`; **copy it to `lib/`**
  for the RTLD Python loader. Populate **`share/basis_sets/`** from `basis_sets/`.
- Run/tests use `openqp/.venv/bin/python` with `OPENQP_ROOT=<worktree>`,
  `PYTHONPATH=<worktree>/pyoqp`, `OMP_NUM_THREADS=1`.
- ddX is **not installed here** (`OQP_ENABLE_DDX` off), so the ddX solve cannot
  run locally; the PCM-on test is skip-guarded.

## Verified status

- Compiles + links. Full Python suite passes; the ddX-gated PCM-on test
  (`tests/test_pcm_energy_path.py`) is **skipped** (expected, no ddX).
- **PCM-off reproduces the vacuum SCF exactly** (H2O RHF/6-31G* = `-76.0107465133`);
  `pcm.enabled=false` is identical with no PCM term.

## Two open workstreams (do these before any gradient/TDDFT/MRSF/excited-state)

1. **ddX-enabled validation** — see `docs/solvent_pcm_validation_gate.md`.
   Build with `-DENABLE_DDX=ON -DDDX_ROOT=<install>`, un-skip the closed-loop
   test, and compare energies vs a trusted ddPCM reference (H2O, NH3). Resolve
   the provisional conventions (q sign/scale, phi_cav, psi, radii, `E_pcm` vs
   Fock-embedded `Tr[D·V_pcm]` bookkeeping).
2. **Dual-path reconciliation** — see `docs/solvent_pcm_reconciliation.md`.
   Two PCM paths coexist (private `epcm`/`add_reference_pcm_reaction_field`/
   `present(pcm_reaction_potential)`/`solvent.py`, vs new `e_pcm`/
   `add_pcm_reaction_field`/`pcm_enabled`/Fortran-C ddX adapter). Independently
   gated, so tests pass — but pick a canonical architecture (Python-supplied
   reaction potential, Fortran-side ddX solve, or layered) and collapse to one
   field + one hook.

## Guardrails (standing)

- No PR, no GitHub issue, no push to `origin`/`upstream`. Private only,
  fast-forward only — no `--force`/`--force-with-lease`/`--mirror`.
- Do not unify the two PCM paths yet (decision pending).
- No surrogate non-ddX solver. PCM is **not production-ready** until one
  validated canonical energy path exists.

## Docs map

- `docs/solvent_pcm_handoff.md` — this file (start here).
- `docs/solvent_pcm_validation_gate.md` — provisional conventions + ddX validation gate.
- `docs/solvent_pcm_reconciliation.md` — dual-path reconciliation task + design decision.
- `docs/solvent_ddx_scf_integration_seam.md` — seam mapping + ddX API findings.
- `docs/solvent_pcm_backend_spike.md` — original architecture spike notes.
