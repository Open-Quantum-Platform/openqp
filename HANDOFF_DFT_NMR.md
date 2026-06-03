# OpenQP NMR shielding (RHF/DFT base) — branch handoff

**Branch:** `feat/dft-nmr` (off `main`).
**Tip:** `d5de87f` — *docs: Phase-0 validation summary + magnetic-symmetry derivation note*.
**Home repo:** **private** — `git@github.com:karmachoi/openqp-private.git` (remote `private`).
  This branch was **moved off the public** `origin` (`karmachoi/openqp`) and deleted there
  (which closed public PR #2). NMR is private from here on. (`upstream` = `Open-Quantum-Platform/openqp`.)
**Date:** 2026-05-31

This is the **foundational RHF/DFT** NMR branch. The MRSF-TDDFT extension continues on
**`feat/mrsf-nmr-gate5a`** (branched off this one; lives in worktree
`../nmr-wt`, also private). See that branch's `HANDOFF_NMR.md` for the MRSF state.

---

## 1. Scope (what's on this branch)

Native (Fortran) NMR nuclear magnetic shielding for **closed-shell RHF/DFT**,
**Common Gauge Origin (CGO)** at the center of mass. Not bridged to any external code — the independent reference is
only an external validation oracle.

```
  σ_ts(N) = σ_dia (density expectation)  +  σ_para (response to B, contracted with PSO)
  σ(ppm) = α²·1e6 × … ,   α = 1/137.035999084
```

**Done & validated:**
- **1e integrals** vs the independent reference: angular momentum (`int1e_cg_irxp`), diamagnetic
  (`dia()` to ~6 figs), PSO. PSO off-center antisymmetry bug fixed (`8e24156`):
  `pso_integrals` returns the exact antisymmetric part `A=(M−Mᵀ)/2`.
- **Diamagnetic** shielding, full tensor.
- **Paramagnetic** — uncoupled (vs the independent reference gold) and **Phase-0 coupled** HF/hybrid
  ground-state magnetic response: the coupled CPHF/CPKS response is the
  exact-exchange response of the imaginary/antisymmetric `P^B` (Coulomb = 0,
  fxc = 0, confirmed numerically), scaled by `c_x = HFscale`; fixed-point solve
  reusing the `int2` A−B (`int2_td_data_t int_amb`) exchange image. HF matches the
  independent coupled CGO reference exactly (O para −230.63); PBE coupled ≡ uncoupled; the
  BHHLYP/PBE0 coupling Δ matches the oracle. Gates 0–4 automated and green.

**Run:** `scf_prop=nmr`.

## 2. Build & run

Identical recipe to the MRSF worktree (see `../nmr-wt/HANDOFF_NMR.md` §2–§3 for the
gotchas). In short, from a built checkout of this branch:

```bash
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-15 \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_INSTALL_PREFIX=$(pwd) \
  -DENABLE_PYTHON=ON -DENABLE_OPENMP=ON -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=ON \
  -DLINALG_LIB_INT64=OFF
ninja -C build && ninja -C build install        # install once: needs share/basis_sets
export OPENQP_ROOT=$(pwd) PYTHONPATH=$OPENQP_ROOT/pyoqp
.venv/bin/python -m oqp.pyoqp <input.inp>        # reuse the hub/nmr-wt venv (Python 3.9)
```
Gotchas (all real): `g++-15` for CXX (clang lacks `-fopenmp`); `-DLINALG_LIB_INT64=OFF`
(Accelerate; `ON` fails to link complex BLAS / corrupts — `rm -rf build` if flipped);
`ninja install` once for the Hückel-guess basis data.

**NOTE:** this branch is **not pre-built** here (it's the archival base). The active,
pre-built worktree is `../nmr-wt` on `feat/mrsf-nmr-gate5a`.

## 3. Tests

```bash
export OPENQP_ROOT=$(pwd) PYTHONPATH=$(pwd)/pyoqp
.venv/bin/python -m unittest tests.test_nmr_shielding tests.test_nmr_coupled
```
- `tests/test_nmr_shielding.py` — dia + uncoupled para vs the independent reference; PSO antisymmetry.
- `tests/test_nmr_coupled.py` — Phase-0 gates 0–4 + `c_x` exchange scaling
  (HF / BHHLYP / PBE0 / PBE).
- `tests/fixtures/nmr/cgo_reference.json` — oracle.

## 4. Key files

- `source/integrals/mod_1e_primitives.F90` — `comp_amom_int1_prim`,
  `comp_nmr_dia_int1_prim`, `comp_pso_int1_prim`.
- `source/integrals/int1.F90` — `angular_momentum_integrals`, `nmr_dia_shielding`,
  `pso_integrals` (FULL antisymmetric PSO).
- `source/modules/nmr_shielding.F90` — module: dia + coupled CPHF/CPKS paramagnetic.
- `include/oqp.h` — declares `nmr_shielding`.
- `pyoqp/oqp/library/runfunc.py` (dispatch `nmr`), `pyoqp/oqp/utils/input_checker.py`
  (`SCF_PROPS` ⊇ `nmr`).
- Docs: `source/modules/NMR_SHIELDING_STATUS.md`,
  `source/modules/NMR_MAGNETIC_SYMMETRY_NOTE.md` (why J=0, fxc=0, exchange-only),
  `source/modules/MRSF_NMR_DESIGN.md` (MRSF design draft, gates 0–5).
- Derivation handoff: `source/modules/NMR_GATE_2B_A_DERIVATION.md` (Gate 2b-A —
  formal amplitude-coupled MRSF paramagnetic derivation; design-only, no code).

## 5. Lineage & next

- **`feat/mrsf-nmr-gate5a`** (descendant, private, worktree `../nmr-wt`) carries the
  MRSF-TDDFT work: gates G1–G4, 5a-1 reduction, 5a-2 diamagnetic (active), 2a
  orbital-only paramagnetic (experimental). Next there is **Gate 2b-A** (formal
  amplitude-coupled derivation, design-only) — now written up in
  `source/modules/NMR_GATE_2B_A_DERIVATION.md`: derives the amplitude response
  `X^B`, the Z-vector response `z^B`, the assembly of `γ^{B}_I`, the spin-flip
  selection rules, the single-reference reduction to the validated Phase-0 coupled
  CPHF/CPKS, and the Gate 2b-B implementation handoff (files, oracles, risks,
  acceptance criteria).
- Possible future on *this* base branch: GIAO (gauge-including atomic orbitals) to
  remove the CGO origin dependence.

## 6. House rules

- **Push to `private` only** — never the public `origin` (NMR was moved off it).
- Commit messages end with `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.
- Keep PSO antisymmetric (`A=(M−Mᵀ)/2`); don't change the α² prefactor or the
  `(r−c)/r³` field factor (the off-center H-para bug).
- Don't touch sibling worktrees/branches (solvent / GPU / EKT).
