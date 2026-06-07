# 001: ddX API probe

## Question

Given OpenQP needs an optional continuum-solvent backend, when we inspect and run ddX/pyddx, then can it support the first OpenQP target: RHF/ROHF energy-only `reference_scf` PCM coupling, with a path toward Fock/KS contributions and later forces?

## Why this matters

The previous scaffold reserved `[pcm]` input keywords but deliberately blocked runtime use. The next decision is whether ddX is a better first backend than PCMSolver before adding build-system and SCF-loop coupling.

## Research findings

Upstream ddX documents the relevant host-code workflow:

- Models: ddCOSMO, ddPCM, ddLPB.
- Interfaces: Fortran core, C interface, Python interface (`pyddx`).
- Host-code target: QM, QM/MM, polarizable MM, MM, or standalone.
- Advertised outputs: electrostatic solvation energy, forces, and contribution to the Fock/DFT Kohn-Sham matrix.
- C API includes model allocation, setup/solve/adjoin solve, energy, solvation force terms, and a `ddx_ddrun` convenience routine.

Important Python API shape from upstream example:

```python
model = pyddx.Model("pcm", centres_bohr, radii_bohr, solvent_epsilon=78.3553)
solute_field = model.multipole_electrostatics(solute_multipoles)
solute_psi = model.multipole_psi(solute_multipoles)
state = pyddx.State(model, solute_psi, solute_field["phi"])
state.solve()
state.solve_adjoint()
energy = 0.5 * np.sum(state.x * solute_psi)
force = state.solvation_force_terms(solute_field)
force += state.multipole_force_terms(solute_multipoles)
```

The C API exposes the lower-level operations OpenQP would likely call from Fortran/C bindings:

- `ddx_allocate_model(...)`
- `ddx_allocate_electrostatics(...)`
- `ddx_setup(...)`
- `ddx_solve(...)`
- `ddx_solve_adjoint(...)`
- `ddx_energy(...)`
- `ddx_solvation_force_terms(...)`
- `ddx_ddrun(...)`

## Local run

I tested `pip install pyddx` in isolated temporary virtual environments on this macOS host.

Result: the wheel installs, but importing fails on both Python 3.9 and Python 3.11 with:

```text
ImportError: dlopen(.../pyddx.cpython-311-darwin.so, 0x0002): symbol not found in flat namespace '_PyObject_ClearManagedDict'
```

I also built ddX from source with CMake/Ninja and ran its C `examples/run_ddx` target successfully:

```text
status=0
Version: 0.8.0
nsph   =   12
nbasis =   81
ncav   = 1290
```

So the core library builds and runs on this macOS host; the observed blocker is specific to the PyPI Python wheel/import path, not the ddX library itself.

For OpenQP coupling, the C header also exposes these useful state getters after solve/adjoin solve:

- `ddx_get_x(state, nbasis, nsph, x)` for the forward solution.
- `ddx_get_s(state, nbasis, nsph, s)` for the adjoint solution.
- `ddx_get_xi(state, ddx, ncav, xi)` for the adjoint solution projected onto cavity points.
- `ddx_get_cavity(ddx, ncav, cavity)` for cavity data.

This is enough to continue a Fortran/C-level integration spike without requiring `pyddx`.

The probe script in this spike records the expected point-charge example and fails gracefully if `pyddx` cannot import:

```bash
python3.11 spikes/001-ddx-api-probe/probe_pyddx_point_charges.py
```

## Verdict: PARTIAL

### What worked

- ddX is still the strongest backend candidate on API/maintenance grounds.
- Its documented host-code model matches the OpenQP solvent design better than using PySCF as a runtime dependency.
- It has a plausible path not only for energy, but eventually for force/Fock/KS contributions.
- The ddX core C/Fortran library built and its C example ran successfully on this macOS host.
- The C API exposes solution/cavity getters that support a Fortran/C-level OpenQP integration spike.
- The first OpenQP target remains `reference_scf` energy-only RHF/ROHF coupling.

### What did not work

- The PyPI `pyddx` macOS wheel does not import on this host, at least with the tested Python 3.9 and Python 3.11 venvs.
- I did not yet validate a live `pyddx` Python energy number locally because the import failed before the example could run.

### Surprises

- The upstream docs and example differ slightly in `Model` constructor shape: docs mention sphere charges as a mandatory argument, while the current example uses `Model("pcm", centres, rvdw, solvent_epsilon=...)`. The source/example should be treated as authoritative for current pyddx.
- ddX's Python example is multipole/sphere based, while OpenQP's SCF coupling will need a density-derived electrostatic potential and reaction-field matrix. That mapping still needs a focused source-level read before implementation.

### Recommendation for the real build

Proceed with ddX-first, but add one more narrow build-system/API spike before touching OpenQP SCF internals:

1. Add optional `ENABLE_DDX`/`ENABLE_SOLVENT` CMake plumbing and a small compile/link smoke test against `ddx.h`/`libddx`.
2. Build a minimal OpenQP-side adapter around the C API: allocate model, setup, solve, retrieve `x/s/xi`, and deallocate cleanly.
3. Identify the exact OpenQP AO-integral path for building the reaction-field matrix from ddX outputs.
4. Then wire the adapter into RHF/ROHF SCF for `runtype=energy` only.

Keep PCMSolver as fallback if ddX's Fock matrix contribution interface proves too hard to map onto OpenQP's AO integrals.
