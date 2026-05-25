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

This looks like a macOS/PyPI wheel compatibility issue rather than an API-design blocker. The next runtime test should use one of:

1. conda-forge `pyddx`, or
2. source build of ddX/pyddx, or
3. C/Fortran build directly, skipping the Python wheel.

The probe script in this spike records the expected point-charge example and fails gracefully if `pyddx` cannot import:

```bash
python3.11 spikes/001-ddx-api-probe/probe_pyddx_point_charges.py
```

## Verdict: PARTIAL

### What worked

- ddX is still the strongest backend candidate on API/maintenance grounds.
- Its documented host-code model matches the OpenQP solvent design better than using PySCF as a runtime dependency.
- It has a plausible path not only for energy, but eventually for force/Fock/KS contributions.
- The first OpenQP target remains `reference_scf` energy-only RHF/ROHF coupling.

### What did not work

- The PyPI `pyddx` macOS wheel does not import on this host, at least with the tested Python 3.9 and Python 3.11 venvs.
- I did not yet validate a live ddX energy number locally because the import failed before the example could run.

### Surprises

- The upstream docs and example differ slightly in `Model` constructor shape: docs mention sphere charges as a mandatory argument, while the current example uses `Model("pcm", centres, rvdw, solvent_epsilon=...)`. The source/example should be treated as authoritative for current pyddx.
- ddX's Python example is multipole/sphere based, while OpenQP's SCF coupling will need a density-derived electrostatic potential and reaction-field matrix. That mapping still needs a focused source-level read before implementation.

### Recommendation for the real build

Proceed with ddX-first, but add one more narrow build-system/API spike before touching OpenQP SCF internals:

1. Build/install ddX from source or conda-forge on this host.
2. Run the point-charge example and record the reference energy/force shape.
3. Inspect the Fortran/C interface for the Fock/KS contribution path and identify the exact arrays OpenQP must provide.
4. Only then add optional `ENABLE_DDX`/`ENABLE_SOLVENT` CMake plumbing and a small import/link smoke test.

Keep PCMSolver as fallback if ddX's Fock matrix contribution interface proves too hard to map onto OpenQP's AO integrals.
