# PCM ddX literature/reference validation benchmarks

This document defines the **scientific** validation gate for OpenQP's energy-only
ddPCM path. It exists so that "PCM works" can only be claimed once the computed
solvation energy matches a literature or independently-generated reference value
under an **exactly matched** protocol — not merely because `e_pcm` is finite and
nonzero.

## Implementation-path invariant (do not regress)

The production PCM physics is **Fortran-side**. The runtime path is:

```
infos%control%pcm_enabled  ->  add_pcm_reaction_field (source/solvent_pcm.F90)
                           ->  E%e_pcm  ->  Fortran/C ddX adapter (source/solvent_ddx_adapter.c)
```

Python may own input/configuration and may parse output to compare numbers.
**Python must not compute the production reaction field, the Fock contribution,
or the solvation energy.** Reference values in this document and in
`tests/data/pcm_literature_benchmarks.json` are comparison targets only; they do
not feed the SCF.

## Current validation status

**NOT physically validated.** No benchmark below has yet been matched against a
reference value, because ddX is not available in the current environment and
because several conventions in the present adapter are still provisional (see
`docs/solvent_pcm_validation_gate.md`):

- `q_cav` sign/scale into `external_charge_potential` (only finite-difference
  evidence for `dE/dphi_cav ≈ -0.5·q_cav`, not the Fock-operator charge);
- `phi_cav` sign split (`phi_total = Σ_k Z_k/|r-R_k| − phi_elec`);
- ddX `psi` is built from **nuclear monopoles only**, not the electronic density,
  so the energy expression `0.5·⟨xs,psi⟩` is incomplete for a QM solute;
- cavity radii come from a small built-in Bondi table, not a curated UFF/Bondi set;
- `E_pcm` (= ddX `esolv`) vs the host-side `½·Tr[D·V_pcm]` bookkeeping is unresolved.

These provisional items are exactly what the benchmarks must expose. A benchmark
that fails when ddX is enabled is a **correct gate result**, not a test bug.

## Reference lineage (ddPCM / ddX)

The ddPCM domain-decomposition lineage, as listed on the ddX references page
(ddX provides ddCOSMO, ddPCM, and ddLPB, and computes the electrostatic
solvation energy plus force and Fock/Kohn-Sham contributions for host codes):

1. **Stamm, Cancès, Lipparini, Maday (2016)** — ddPCM discretization within the
   domain-decomposition paradigm.
2. **Gatto, Lipparini, Stamm (2017)** — analytic forces for ddPCM.
3. **Nottoli, Stamm, Scalmani, Lipparini (2019)** — QM calculations in solution
   with ddPCM (the host-code QM coupling this path mirrors).
4. **Mikhalev, Nottoli, Stamm (2022)** — linearly scaling ddPCM energy and forces.

Plus the **ddX project documentation and reference data** for the C/Fortran API
used behind `source/solvent_ddx_adapter.c`.

> Exact bibliographic details (journal/volume/DOI) must be transcribed from the
> ddX references page before a literature number is entered below; they are
> intentionally not hand-typed here to avoid an unverified citation.

## Quantities compared (each validated separately)

| # | Quantity | Source | How obtained |
| --- | --- | --- | --- |
| Q0 | Vacuum total energy is unchanged when PCM is off | Fortran | run with no `[pcm]` / `enabled=false`; equals the established vacuum value |
| Q1 | ddX electrostatic solvation energy `e_pcm` | Fortran (ddX `esolv`) | parsed `PCM solvent energy (prov)` line |
| Q2 | Independent diagnostic `½·Tr[D·V_pcm]` | Fortran (future diagnostic) | not yet exposed; tracked as a gate item, not computed in Python |
| Q3 | Final total energy with PCM | Fortran | parsed `TOTAL energy` line |
| Q4 | SCF convergence behaviour | Fortran | run returns 0, converged |
| Q5 | `q_cav` sign/scale | ddX adapter | cross-check vs reference / FD relation |
| Q6 | `phi_cav` sign convention | Fortran | cross-check vs reference |
| Q7 | electronic contribution to ddX `psi` | adapter | confirm `psi` includes density, not just nuclear monopoles |
| Q8 | cavity/radii/model choices | adapter | recorded and matched to the reference protocol |

The first machine-checked quantity is **Q1** (`e_pcm`), with **Q0** as the
mandatory no-regression guard and **Q3** as a consistency cross-check. Q2 and
Q5–Q8 are convention items resolved as references are matched.

## How a reference value is made admissible

A row may carry a numeric `reference_value` **only** when all of the following
are pinned to the same values OpenQP uses for that run:

- geometry (and units) and its source,
- charge / multiplicity,
- method and basis,
- solvent dielectric `epsilon`,
- cavity/radii convention,
- ddPCM discretization (`model`, `lmax`, `n_lebedev`, `eta`, FMM),
- the energy quantity and its sign/scale convention.

Two admissible reference kinds:

1. **Independent ddX cross-check** — drive the *same* geometry/radii/discretization
   through an independent ddX front end (e.g. `pyddx`) and record its
   `esolv`. Tight tolerance (≈ `1e-5` Eh) because it is the same physics, same
   discretization, different driver. This is the recommended first reference and
   directly isolates OpenQP↔ddX coupling bugs (Q5–Q7).
2. **Published ddPCM QM value** — a number from the lineage above, transcribed
   only when the paper's protocol is reproduced exactly. Looser tolerance
   (≈ `1e-3` Eh / few × 0.1 kcal/mol) reflecting protocol/geometry uncertainty.

Do **not** use bare experimental hydration free energies as pass/fail targets;
they include cavitation/dispersion/repulsion and thermodynamic-state terms this
electrostatic-only path does not model.

## Benchmark target table

Discretization columns reflect the **current adapter defaults**
(`source/solvent_ddx_adapter.c :: build_pcm_model`): `model=ddpcm`, `lmax=8`,
`n_lebedev=302`, `eta=0.1`, FMM on; radii from the built-in Bondi-style table
(provisional). `epsilon=78.3553` is water.

| ID | Molecule | Geometry source | chg/mult | Method/basis | ε | Radii | ddPCM params | Reference quantity | Reference source | Tolerance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| h2o_rhf_631gs_water | H₂O | repo `examples/HF/H2O` geometry (Å) | 0 / 1 | RHF / 6-31G* | 78.3553 | Bondi (built-in, provisional) | ddpcm, lmax=8, nleb=302, eta=0.1, FMM | `e_pcm` (ddX esolv) | pending — pyddx cross-check, same protocol | 1e-5 Eh |
| nh3_rhf_631gs_water | NH₃ | standard C₃ᵥ equilibrium (Å), recorded in JSON | 0 / 1 | RHF / 6-31G* | 78.3553 | Bondi (built-in, provisional) | ddpcm, lmax=8, nleb=302, eta=0.1, FMM | `e_pcm` (ddX esolv) | pending — pyddx cross-check, same protocol | 1e-5 Eh |
| hf_rhf_631gs_water | HF | r(HF)=0.917 Å (recorded in JSON) | 0 / 1 | RHF / 6-31G* | 78.3553 | Bondi (built-in, provisional) | `e_pcm` (ddX esolv) | pending — low-polarity check, pyddx cross-check | 1e-5 Eh |

`reference_value: null` / `status: "pending_reference"` rows are **skipped** by
the validation test even when ddX is available, with a message pointing here.
A row becomes a real pass/fail gate only after its reference is populated and the
protocol is confirmed matched.

## Populating a reference (procedure)

For the recommended ddX cross-check, with a ddX-enabled environment and `pyddx`:

1. Build the cavity for the exact geometry and the **same radii** the OpenQP
   adapter assigns (Bondi-by-Z today).
2. Build `phi_cav` from the **same** total solute potential convention OpenQP
   uses, and `psi` consistently (note: OpenQP's `psi` is currently nuclear-only —
   when Q7 is fixed, the reference must move in lockstep).
3. Run `model=ddpcm`, `lmax=8`, `n_lebedev=302`, `eta=0.1`, FMM on.
4. Record `esolv` into `tests/data/pcm_literature_benchmarks.json` with
   `status: "verified"` and the generation provenance.

Until then the gate stays open and PCM remains **not physically validated**.
