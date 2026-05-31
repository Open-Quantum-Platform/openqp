# PCM energy path — provisional conventions & ddX-enabled validation gate

Status as of commit `feat(pcm): wire ddX-gated PCM energy path`.

The first closed RHF/DFT PCM single-point energy loop is **wired and builds**, and
the vacuum path is preserved exactly. It is **NOT yet physically validated**. This
note records what is provisional and the exact gate that must pass before the
implementation may be treated as correct.

## What is implemented

The closed loop, per SCF Fock build:

```
D  ->  phi_cav  ->  ddX q_cav  ->  V_pcm  ->  Fock / E_pcm  ->  D'
```

| Stage | Code |
| --- | --- |
| `pcm_enabled` / `pcm_epsilon` input | `source/types.F90`, `include/oqp.h`, `pyoqp/oqp/molecule/oqpdata.py` |
| Orchestration (the loop) | `source/solvent_pcm.F90` :: `add_pcm_reaction_field` |
| ddX cavity + solve (C) | `source/solvent_ddx_adapter.c` :: `oqp_ddx_pcm_cavity`, `oqp_ddx_pcm_solve` |
| `phi_cav` (electronic) | `int1.F90` :: `electrostatic_potential_unweighted` |
| `V_pcm` AO matrix | `int1.F90` :: `external_charge_potential` |
| Fock hook + `E%e_pcm` | `source/scf_addons.F90` :: `calc_jk_xc` |

Gated on `infos%control%pcm_enabled`. With a non-ddX build the adapter returns
status 2 and the run aborts cleanly at the adapter boundary
(`PCM (ddX) cavity build failed: OpenQP was built without OQP_ENABLE_DDX`).

## Verified so far (ddX unavailable)

- Compiles and links (requires `-DLINALG_LIB_INT64=OFF`; Accelerate rejects ILP64).
- **PCM-off is unchanged**: H2O RHF/6-31G* = `-76.0107465133`; a run with no
  `[pcm]` section and one with `pcm.enabled=false` are identical, no PCM term.
- PCM-on closed-loop test (`tests/test_pcm_energy_path.py`) is **skip-guarded**,
  not failing, when ddX is absent.

## Provisional conventions — DO NOT treat as physically correct

Each is marked in code; none has been checked against a reference.

1. **`q_cav` sign/scale** — `oqp_ddx_pcm_solve` returns ddX `ddx_get_xi`
   (cavity-projected adjoint `state%q`) and it is used *directly* as the
   external-charge vector for `external_charge_potential`. The branch only
   finite-difference-validated `dE/dphi_cav ≈ -0.5·q_cav` for the *energy*, not
   the Fock-operator charge. The needed factor/sign (e.g. `-0.5·q_cav`) is unconfirmed.
2. **`phi_cav` sign** — `phi_total = Σ_k Z_k/|r-R_k| − phi_elec` (`solvent_pcm.F90`).
   The electronic/nuclear sign split is assumed, not verified.
3. **ddX `psi`** — built from **nuclear monopoles only** (`charges/√(4π)`),
   not the electronic density. The energy expression `0.5·⟨xs,psi⟩` is therefore
   incomplete for a QM solute.
4. **Cavity radii** — a small built-in Bondi vdW table keyed by Z
   (`vdw_radius_angstrom`), not an OpenQP-curated/UFF radius set; heavier elements
   fall back to 1.80 Å.
5. **`E_pcm` bookkeeping** — the ddX solvation energy is reported as a separate
   `E%e_pcm` term added to `etot` *after* the vacuum HF energy. Whether this
   double-counts the Fock-embedded `Tr[D·V_pcm]` (or under-counts) is unresolved.
6. **Model/discretization** — ddPCM (`model=2`), `lmax=8`, `n_lebedev=302`,
   `eta=0.1`, FMM on; hard-coded in `build_pcm_model`. Not convergence-tested.
7. **Cost** — the model is rebuilt twice per SCF iteration (cavity, then solve).
   Correctness-only; not optimized.

## ddX-enabled validation gate (the next review)

Do NOT mark PCM physically validated until ALL of the following pass.

1. **Build with ddX.**
   - Obtain/build ddX (ddsolvation/ddX) and configure with
     `-DENABLE_DDX=ON -DDDX_ROOT=<install>` (see `cmake/FindDDX.cmake`) plus the
     existing `-DLINALG_LIB_INT64=OFF`. Confirm `OQP_ENABLE_DDX` is defined and
     `DDX::ddx` links.
2. **Un-skip the closed-loop test.**
   - `tests/test_pcm_energy_path.py::test_pcm_on_closes_loop` must now RUN (not
     skip) and pass: finite nonzero `E%e_pcm`, total moved off the vacuum value,
     SCF converged.
3. **Reference comparison (the real correctness check).**
   - Reproduce a trusted ddPCM/ddCOSMO single-point solvation energy for at least
     H2O and one ion, e.g. via standalone `pyddx` driving the same geometry/radii,
     and/or a published ddX/GAMESS/psi4-PCM number. Agreement to a stated
     tolerance fixes items (1)–(5) above; until then they stay provisional.
   - In particular, settle the `q_cav` factor/sign against the FD relation and the
     reference, and decide the `E_pcm` vs Fock-embedded `Tr[D·V_pcm]` bookkeeping.
4. **psi from the electronic density** — replace the nuclear-monopole `psi` with a
   QM-density-derived `psi` (or confirm the monopole approximation is acceptable
   for the energy at the chosen tolerance).
5. **Radii decision** — adopt a documented radius set (UFF/Bondi) and record it.

Only after (1)–(3) pass should follow-on scope (gradients, TDDFT, MRSF,
excited-state solvation) be considered — all still explicitly out of scope here.

## Related

- `docs/solvent_ddx_scf_integration_seam.md` — seam mapping and ddX API findings.
- `docs/solvent_pcm_backend_spike.md` — original architecture spike notes.
