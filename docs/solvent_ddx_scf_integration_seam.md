# ddX SCF integration seam

This note identifies the first OpenQP code seam for turning ddX reaction-field output into an RHF/ROHF energy-only PCM contribution.

## Existing OpenQP pieces

### SCF Fock/energy builder

The main SCF loop in `source/scf.F90` calls:

```fortran
call calc_fock(basis, infos, molgrid, pfock, energy, mo_a, pdmat, mo_b, nschwz, fold, dold)
```

`calc_fock` is implemented in `source/scf_addons.F90`. It copies density matrices from tagarray or the SCF loop override, then calls `calc_jk_xc`.

`calc_jk_xc` is the narrowest first insertion point for energy-only PCM because it already has:

- packed AO density matrices `d(:, :)`
- packed AO Fock matrices `f(:, :)`
- number of spin blocks `nfocks`
- the one-electron core Hamiltonian `hcore`
- energy accumulator `E`

Today it does:

```fortran
call fock_jk(...)
do ii = 1, nfocks
  f(:,ii) = f(:,ii) + hcore
end do
...
E%ehf1 = sum_i Tr[D_i Hcore]
E%ehf  = 0.5 * (sum_i Tr[D_i F_i] + E%ehf1)
E%etot = E%ehf + E%nenergy
```

For `reference_scf` PCM, the solvent reaction-field matrix should be treated as an additional density-dependent one-electron operator:

```fortran
V_pcm[D] = reaction_field_matrix_from_ddx(D_total)
F_i      = F_i + V_pcm[D]
E_pcm    = 0.5 * Tr[D_total V_pcm[D]]    ! or backend-equivalent expression
```

For ROHF/UHF-style two-spin storage, the reaction field is scalar/electrostatic and should be added to both spin Fock blocks. The source density for ddX is the total density:

```fortran
D_total = D_alpha + D_beta
```

For RHF, OpenQP stores one density block; existing energy expressions indicate that block is already the density used in traces.

### Existing one-electron electrostatic-potential kernel

`source/integrals/int1.F90` already had a private kernel:

```fortran
int1_el_pot(basis, x, y, z, d, pot, tol)
```

This evaluates the electronic electrostatic potential from a packed AO density at arbitrary points. The older public `electrostatic_potential` wrapper multiplies the result by grid weights and is therefore appropriate for quadrature-style property integrations, not for ddX `phi_cav`.

The branch now adds a safe unweighted public wrapper:

```fortran
electrostatic_potential_unweighted(basis, x, y, z, d, pot, logtol)
```

This normalizes the density before calling `int1_el_pot`, restores density normalization afterwards, and deliberately does not multiply by weights. It is the intended OpenQP entry point for electronic `phi_cav` construction.

### Existing one-electron external-charge integral kernel

`source/integrals/int1.F90` already had a private kernel:

```fortran
int1_coul_ext_chg(h, basis, nat, x, y, z, chg, tol, chgtol)
```

This computes packed one-electron Coulomb integrals from external point charges:

```text
sum_k q_k <mu | 1/|r-r_k| | nu>
```

The branch now adds a normalized public wrapper:

```fortran
external_charge_potential(basis, v, x, y, z, chg, logtol, chgtol)
```

This wrapper zeroes `v`, calls the existing kernel, and applies `bas_norm_matrix` so it returns a packed AO matrix in the same normalization convention as `Hcore` and `Fock`.

This is the natural first path for converting ddX apparent charges/cavity-point charges into a PCM reaction-field matrix.

## Proposed first runtime implementation

1. Add a Fortran/C bridge wrapper around the C ddX adapter for the quantities needed in SCF.
2. At the start of a PCM-enabled `calc_jk_xc` call, form `D_total` from packed density blocks. The branch now has the dependency-light Python guard/helper `reference_scf_total_density()` for the intended convention: one RHF density block is used as-is, while ROHF/UHF-style alpha/beta packed blocks are summed elementwise before building ddX `phi_cav`.
3. Use `int1:electrostatic_potential_unweighted` to compute electronic potential on ddX/cavity sites.
4. Add the nuclear potential at those sites to form total `phi_cav`.
5. Use ddX to solve the reaction field and retrieve the point representation needed for OpenQP.
6. Use `external_charge_potential` to build packed `V_pcm` from reaction charges/sites.
7. Add `V_pcm` to all spin Fock blocks before diagonalization/convergence acceleration. The dependency-light helper `reference_scf_reaction_fock_updates()` now records this intended replication: one packed RHF block receives one update, while ROHF two-block storage receives the same scalar/electrostatic reaction-potential update on both alpha and beta Fock blocks. The branch also has a source-level Fortran handoff helper, `add_reference_pcm_reaction_field(f, reaction_potential, nfocks)`, that applies the same packed reaction field to each requested Fock block after the caller has validated the packed AO matrix contract.
8. Add a dedicated PCM energy term to `scf_energy_t` rather than hiding it inside `ehf1`; this avoids confusing final energy breakdowns.

## Deliberate non-goals for the first implementation

- No gradients/optimization.
- No state-specific or nonequilibrium MRSF solvent response.
- No SMD surface terms.
- No incremental-Fock shortcut for PCM until the full non-incremental energy path is validated.

## ddX mapping findings from source inspection

The relevant ddX state fields are documented in `src/ddx_core.f90`:

- `state%phi_cav`: electric potential at cavity points, dimension `ncav`; used to build the primal RHS.
- `state%psi`: representation of the solute density/potential in spherical harmonics, dimension `(nbasis, nsph)`; used as RHS for the adjoint linear system.
- `state%xs`: forward ddCOSMO solution, dimension `(nbasis, nsph)`; used by COSMO and PCM.
- `state%s`: adjoint ddCOSMO solution, dimension `(nbasis, nsph)`.
- `state%q`: effective total adjoint ddPCM solution, defined in ddX as
  `Q = S - 4*pi/(eps - 1) * Y`; this is only populated by the ddPCM adjoint solve.

For COSMO/PCM setup, ddX does:

```fortran
call cav_to_spherical(..., phi_cav, state%phi)
state%phi = -state%phi
state%phi_cav = phi_cav
state%psi = psi
state%rhs_done = .true.
state%adjoint_rhs_done = .true.
```

The ddX solvation energy is:

```fortran
esolv = 0.5 * dot(state%xs, state%psi)
```

The Python/C++ API confirms the expected host-code call shape for a general QM host:

```python
state = pyddx.State(model, psi, phi)
state.solve()
state.solve_adjoint()
energy = state.energy()
```

For the OpenQP SCF path this means:

1. `phi_cav` should be the unweighted total solute electrostatic potential at ddX cavity points.
2. `psi` should be ddX's spherical-harmonic representation of the same solute source. ddX exposes helpers for point multipoles (`ddx_multipole_psi`), but not yet a generic C helper that converts an arbitrary QM density potential into `psi` directly.
3. The forward solution `xs` is enough for the ddX energy expression, but ddX documentation says the adjoint solve is required for the Fock/Kohn-Sham operator contribution. Therefore the first SCF implementation should not assume `xs` alone is the reaction-field charge vector for `external_charge_potential`.
4. For ddPCM, `state%q` is the physically relevant effective adjoint quantity for derivatives/Fock-like response. The public C API routine `ddx_get_xi` projects ddPCM `state%q` onto cavity points through `ddproject_cav`; it is therefore better interpreted in OpenQP as a cavity-projected `state%q`/`q_cav` vector than as generic adjoint `s`.

## Current implementation status

The branch now includes two ddX adapter smoke paths:

1. `oqp_ddx_run_point_charge_smoke`: the original high-level ddX point-charge lifecycle using `ddx_ddrun`.
2. `oqp_ddx_run_explicit_pcm_smoke`: an explicit host-code PCM path that builds host-supplied `psi` and `phi_cav`, then calls `ddx_pcm_setup`, `ddx_pcm_solve`, `ddx_pcm_solve_adjoint`, `ddx_pcm_energy`, and `ddx_get_xi` to retrieve the cavity-projected `state%q` as `q_cav_norm` in the smoke result.
3. `oqp_ddx_run_explicit_pcm_reaction_field_smoke`: the same explicit PCM path, but it also copies caller-owned arrays `cavity_xyz[3*ncav]` and `q_cav[ncav]`. This matches the OpenQP-side `external_charge_potential(basis, v, x, y, z, chg, ...)` handoff shape: split `cavity_xyz` into x/y/z arrays and pass a sign/scale-adjusted candidate external-charge vector.
4. A finite-difference check perturbing `phi_cav[0]` shows the local derivative obeys `dE/dphi_cav[0] ≈ -0.5*q_cav[0]` for the smoke case. In the current smoke, direct use of `q_cav[0]` misses the derivative by about `0.0000971173`, while the `-0.5*q_cav[0]` relation is within about `0.0000000140`.

The explicit path is still a smoke/proof-of-seam, not production OpenQP SCF coupling. It verifies that OpenQP can drive the lower-level ddPCM setup/forward/adjoint API that will be needed after `psi` and `phi_cav` come from AO density and nuclear potentials.

### Reviewed payload and calc_fock handoff chain

The branch also has a dependency-light Python handoff chain for carrying a reviewed reference-SCF reaction field toward the SCF Fock builder without enabling runtime PCM:

1. `reference_scf_pcm_runtime_payload(density_blocks, reaction_potential)` validates the RHF/ROHF reference density and packed AO reaction potential, then stores only `OQP::pcm_reaction_potential`, matched `OQP::pcm_epcm`, `nbf`, `packed_ao_length`, `expected_packed_ao_length`, `packed_ao_shape_formula`, payload version, first-scope metadata, and `backend_validation_status="pending PySCF/ddX/reference cross-check"`.
2. `reference_scf_pcm_reaction_potential_from_payload(payload)` is the consumer gate. It rejects non-mapping inputs with `PCM runtime payload must be a mapping` before field-level validation. For mapping payloads, the shape metadata (`nbf`, `packed_ao_length`, `expected_packed_ao_length`, and `packed_ao_shape_formula`) is required, not optional defaults; the consumer recomputes the packed AO length and rejects stale shape metadata, missing `backend_validation_status`, and boolean numeric payload values before exposing `pcm_reaction_potential_in`.
3. `reference_scf_pcm_calc_fock_handoff(payload)` exposes only `calc_fock_kwargs={"pcm_reaction_potential_in": ...}` plus compact metadata for a future opt-in prototype `calc_fock(..., pcm_reaction_potential_in=...)` call.
4. `reference_scf_pcm_calc_fock_handoff_from_molecule(mol)` is the molecule-level no-runtime gate: an empty/restored payload returns empty `calc_fock_kwargs`, while a present payload must pass the reviewed consumer before any packed reaction potential is exposed.
5. `reference_scf_pcm_calc_fock_request(mol, incremental_fock=...)` is the first caller-facing request guard. It forwards reviewed payloads only to the non-incremental Fock path, returns a labeled disabled request when no payload is present, and fails fast with `reference PCM incremental Fock is not validated` if a restored reference-PCM payload would otherwise enter OpenQP's existing incremental-Fock shortcut.
6. `reference_scf_pcm_calc_fock_request_from_scf_state(mol, dens_old=..., f_old=...)` is the SCF-state shim for a future call site. It derives the incremental-Fock flag from actual old-buffer state (`dens_old is not None or f_old is not None`) before delegating to the request guard, so reviewed reference-PCM payloads are blocked whenever either old density or old Fock state is present, while no-payload molecules still return an explicit disabled request. Blocked reviewed payload diagnostics preserve the old-buffer provenance before exposing `pcm_reaction_potential_in`: `dens_old_present=true/false`, `f_old_present=true/false`, and ordered trigger labels such as `incremental_trigger_fields=dens_old,f_old`, `incremental_trigger_fields=dens_old`, or `incremental_trigger_fields=f_old`.
7. `reference_scf_pcm_calc_fock_call_site_bridge()` is the final dependency-light SCF call-site staging helper before native wiring. It preserves the same disabled/no-payload behavior for empty molecule payloads, forwards a reviewed payload only as `pcm_reaction_potential_in` on the non-incremental path, and records `forward_pcm_reaction_potential` so downstream prototype callers can audit whether a packed reaction field would be sent to `calc_fock` without enabling runtime PCM.

This is intentionally not a production solvent switch: runtime PCM remains disabled, with no state-specific or nonequilibrium MRSF solvent response, no MRSF-kernel solvent response, no analytic PCM gradients, and no solution-phase optimization claim.

## Consequence for the next implementation step

Before adding a production SCF hook, the remaining ddX/OpenQP mapping question is now narrower:

A. Preferred: for the first prototype AO reaction-field matrix, test `chg = -0.5*q_cav` as the finite-difference-consistent derivative of the ddPCM energy with respect to `phi_cav`; keep this explicitly marked as a candidate until cross-checked against a package/reference implementation.

B. If upstream ddX already has a public routine for the Fock/KS contribution that was missed, wrap that routine directly instead of reconstructing charges in OpenQP.

C. For a temporary COSMO-only prototype, test whether `xs` projected to cavity points can reproduce a known ddCOSMO/COSMO Fock contribution, but keep it explicitly marked as a prototype until validated.

The OpenQP-side AO matrix seam is now ready: once a physically correct point/cavity reaction representation is available, use `external_charge_potential` to build the packed AO matrix and add it to all spin Fock blocks.

## Open questions before coding the SCF hook

- What exact ddX quantity should be exposed for the Fock/Kohn-Sham operator: `q`, `qgrid`, `xi`, or a higher-level ddX routine?
- Whether ddX wants `psi` for arbitrary QM density from spherical projection of the potential, direct source-density projection, or a combination of electronic+nuclear multipoles.
- Whether OpenQP's public `electrostatic_potential_unweighted` wrapper should be extended to include nuclear potential assembly, or whether nuclear `phi_cav` should be built separately in the solvent adapter.
- Where solvent metadata and energy fields should be stored in OpenQP JSON/database output.
