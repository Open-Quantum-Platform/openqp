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
2. At the start of a PCM-enabled `calc_jk_xc` call, form `D_total` from packed density blocks.
3. Use `int1:electrostatic_potential` to compute electronic potential on ddX/cavity sites if ddX requires MEP input in cavity-point form.
4. Use ddX to solve the reaction field and retrieve the point representation needed for OpenQP.
5. Use `external_charge_potential` to build packed `V_pcm` from reaction charges/sites.
6. Add `V_pcm` to all spin Fock blocks before diagonalization/convergence acceleration.
7. Add a dedicated PCM energy term to `scf_energy_t` rather than hiding it inside `ehf1`; this avoids confusing final energy breakdowns.

## Deliberate non-goals for the first implementation

- No gradients/optimization.
- No state-specific or nonequilibrium MRSF solvent response.
- No SMD surface terms.
- No incremental-Fock shortcut for PCM until the full non-incremental energy path is validated.

## Open questions before coding the SCF hook

- Which ddX output should be mapped to reaction charges for `external_charge_potential`: forward solution `x`, projected adjoint/cavity vector `xi`, or another cavity representation from the ddX API? The point-charge smoke test proves lifecycle and retrieval, not yet the physical mapping for QM density.
- Whether ddX wants the electronic MEP weighted (`electrostatic_potential` currently multiplies by `wt`) or unweighted (`int1_el_pot` is private and returns unweighted values).
- Where solvent metadata and energy fields should be stored in OpenQP JSON/database output.
