# PCM solvent backend spike

This note records the recommended first steps for adding continuum-solvent support to OpenQP without claiming full excited-state PCM or PCM gradients prematurely.

## Recommendation

Use a two-track approach:

1. Use PySCF as the reference/prototyping engine for regression data.
2. Prefer ddX as the first optional OpenQP backend candidate, with PCMSolver kept as the classic PCM fallback candidate.

The first production scope should be energy-only SCF coupling for RHF/ROHF references. For MRSF-TDDFT, this means a PCM-solvated high-spin ROHF reference followed by the existing MRSF response calculation. It is not full state-specific or nonequilibrium excited-state PCM.

## Backend comparison

| Candidate | Role | Advantages | Concerns |
| --- | --- | --- | --- |
| ddX | Preferred backend spike | Active project; Fortran core; C/Fortran/Python interfaces; ddCOSMO, ddPCM, ddLPB; documents Fock/Kohn-Sham and force contributions for host codes | Need to map its host-code interface onto OpenQP density/Fock data; model differs from classic IEFPCM defaults |
| PCMSolver | Classic PCM fallback | Host-code API is simple: provide surface MEP, solve apparent surface charges, compute polarization energy | Less active; C++ dependency; OpenQP still needs Fock matrix and gradient coupling |
| PySCF solvent | Reference only | Fastest way to generate expected energies and check conventions for PCM/ddCOSMO/ddPCM/SMD | Not a good OpenQP runtime dependency because it is a full QC stack with its own molecule, basis, integral, and SCF objects |

## Initial user-facing input shape

The branch reserves a `[pcm]` section but intentionally blocks `enabled=true` until runtime coupling exists:

```ini
[pcm]
enabled=false
backend=ddx
mode=reference_scf
model=ddpcm
solvent=water
epsilon=78.3553
radii=uff
```

Planned modes:

- `reference_scf`: self-consistent PCM on RHF/ROHF reference; first target.
- `post_state_correction`: perturbative MRSF state-density polarization correction; later and off by default.
- `reference_scf_plus_post_state`: combined mode; later, after double-counting is defined and validated.

## First implementation scope

Implement only:

- `runtype=energy`
- RHF and ROHF references
- optional backend dependency
- self-consistent electrostatic PCM contribution to the SCF Fock/Kohn-Sham matrix
- PCM polarization energy bookkeeping in output/JSON metadata
- MRSF-TDDFT support only through a PCM-solvated ROHF reference

Do not implement in the first PR:

- analytic PCM gradients
- geometry optimization with PCM
- MRSF state-specific self-consistent PCM
- nonequilibrium/linear-response excited-state PCM
- SMD empirical terms

## SCF coupling sketch

For each SCF iteration:

1. Build the current total density. For ROHF, use `D_total = D_alpha + D_beta`; the scalar reaction field is not spin-dependent.
2. Compute the solute electrostatic potential on backend cavity/grid sites from nuclei and electronic density.
3. Ask the backend to solve the solvent reaction field/apparent surface charges.
4. Build the one-electron reaction-field matrix and add it to the alpha/beta Fock-like matrices.
5. Include `E_pcm = 1/2 * V_total dot q_pcm` or the backend-equivalent polarization-energy expression in total-energy bookkeeping.

## Validation targets

Generate reference data before enabling runtime support:

- H2O RHF/6-31g* in water
- formaldehyde BHHLYP/6-31g* in water
- a small ROHF triplet/radical reference in water
- one tiny MRSF-TDDFT case with and without the PCM-solvated ROHF reference

Record:

- gas-phase total energy
- solvated total energy
- PCM polarization energy if exposed
- backend, model, dielectric, cavity/radii settings
- SCF iteration count
- final dipole if available

## Current branch changes

This branch only adds input-level scaffolding, tests, optional backend link plumbing, and a disposable ddX API spike:

- `[pcm]` schema defaults in `pyoqp/oqp/molecule/oqpdata.py`
- input checker validation and guardrails in `pyoqp/oqp/utils/input_checker.py`
- optional `ENABLE_DDX` CMake link plumbing in `CMakeLists.txt`, `cmake/FindDDX.cmake`, and `source/CMakeLists.txt`
- OpenQP-owned ddX adapter API in `source/solvent_ddx_adapter.c` and `source/solvent_ddx_adapter.h`
- tests in `tests/test_pcm_scaffold.py` and `tests/test_ddx_cmake_scaffold.py`
- C link smoke test source in `tests/ddx_link_smoke.c`
- C adapter smoke test source in `tests/ddx_adapter_smoke.c`
- ddX API probe notes/script in `spikes/001-ddx-api-probe/`

The guardrail is deliberate: `pcm.enabled=true` currently produces an input-check error because no runtime Fock/energy coupling has been implemented yet.

## ddX probe outcome

The ddX API still looks like the best first backend candidate. Local `pip install pyddx` produced a macOS import failure in isolated Python environments (`_PyObject_ClearManagedDict` missing), but the core ddX C/Fortran library built from source with CMake/Ninja and its C `run_ddx` example ran successfully.

The branch now includes optional `ENABLE_DDX` CMake plumbing:

- `cmake/FindDDX.cmake` finds `ddx.h` and `libddx` from `DDX_ROOT`/`$DDX_ROOT`.
- `ENABLE_DDX=ON` requires the backend and links `oqp` against imported target `DDX::ddx`.
- A `oqp_ddx_link_smoke` CTest executable verifies compile/link/runtime access to `ddx_get_banner` without touching SCF internals.
- A `oqp_ddx_adapter_smoke` CTest executable exercises the OpenQP-owned adapter lifecycle: allocate model/electrostatics/state, build multipole RHS, solve, retrieve `x/s/xi/cavity`, validate energy/norms, and deallocate.

Verified locally with a ddX source build under `/tmp/ddx-openqp-smoke`:

```bash
cmake -S . -B /tmp/openqp-ddx-on-config -G Ninja \
  -DENABLE_DDX=ON -DDDX_ROOT=/tmp/ddx-openqp-smoke \
  -DENABLE_PYTHON=OFF -DUSE_LIBINT=OFF -DLINALG_LIB=none
cmake --build /tmp/openqp-ddx-on-config --target oqp_ddx_adapter_smoke
ctest --test-dir /tmp/openqp-ddx-on-config -R 'oqp_ddx_(link|adapter)_smoke' --output-on-failure
```

Result: both `oqp_ddx_link_smoke` and `oqp_ddx_adapter_smoke` passed. The next OpenQP step is to identify the AO-integral path for turning ddX reaction-field output into an RHF/ROHF Fock contribution.
