# SOC-NAMD-QMMM production guidance

**Status:** implemented and benchmarked  
**Recommendation:** use `soc_basis=mch` for production SOC-NAMD-QMMM

This note records the current SOC-NAMD-QMMM implementation status, the measured
energy-conservation behavior on the equilibrated formaldehyde/water PBC test,
and the recommended input settings.

## Implemented methods

OpenQP currently supports the following MRSF-TDDFT surface-hopping modes:

| Mode | Main switches | Force model | Intended use |
|------|---------------|-------------|--------------|
| IC-NAMD | `soc=False` | exact active-state MRSF gradient | production IC dynamics |
| SOC option 1 | `soc=True`, `soc_basis=adiabatic` | SHARC-like spin-adiabatic weighted MCH gradient | comparison and diagnostics |
| SOC option 2 | option 1 + `soc_du_dt_corr=True` | option 1 plus finite-difference `dU/dt` correction | experimental diagnostics |
| SOC option 2b | option 1 + `soc_tdc_grad_corr=True` | option 1 plus TDC-projected NAC-gradient correction | experimental diagnostics |
| SOC option 3 | `soc=True`, `soc_basis=mch` | MCH-basis propagation with exact active-root gradient | preferred production SOC mode |

Option 3 is preferred because the nuclear force is the analytic gradient of the
active spin-pure MCH state.  The spin-adiabatic options hop on SOC-diagonal
states, but their gradients remain approximations unless all MCH derivative
couplings and SOC-gradient matrix elements are available.

## Recommended production input

For solvated formaldehyde SOC-NAMD-QMMM, use an equilibrated PBC water box and:

```ini
[md]
soc=True
soc_basis=mch
init_state=S1
thrshe=0.1
dt=0.5
substep=20
decoherence=edc
velocity=maxwell
seed=101
```

For QM/MM:

```ini
[input]
qmmm_flag=True

[qmmm]
embedding=electrostatic
cutoff=PME
qm_atoms=0-3
```

Use different `seed` values for independent trajectories.  Keep the solvent
box, force field, timestep, initial state, and QM region fixed when comparing
SOC options.

## Energy-conservation validation

The production validation used formaldehyde plus 338 TIP3P waters in a
22 Angstrom cubic PBC box.  The box was minimized and equilibrated at fixed box
size by restrained NVT followed by unrestrained NVT, then used for five
independent 100 fs trajectories per method with 0.5 fs nuclear timestep.

| Method | Final drift (kJ/mol) | Energy spread (kJ/mol) | Hops/trj |
|--------|----------------------|-------------------------|----------|
| Gas IC-FSSH | -0.19 +/- 1.95 | 4.89 +/- 0.59 | 0.2 +/- 0.4 |
| Gas SOC option 1 | +5.42 +/- 2.95 | 9.94 +/- 3.72 | 28.4 +/- 11.9 |
| Gas SOC option 2 | -3.22 +/- 34.41 | 38.25 +/- 18.30 | 23.2 +/- 8.2 |
| Gas SOC option 3 | +0.21 +/- 0.17 | 0.68 +/- 0.26 | 8.4 +/- 5.0 |
| IC-NAMD-QMMM | +0.53 +/- 3.54 | 5.86 +/- 1.44 | 1.0 +/- 0.7 |
| SOC-QMMM option 1 | +2.31 +/- 1.99 | 5.58 +/- 1.84 | 13.0 +/- 7.7 |
| SOC-QMMM option 2 | -5.02 +/- 11.43 | 21.16 +/- 4.41 | 8.0 +/- 5.0 |
| SOC-QMMM option 3 | +0.89 +/- 1.73 | 2.96 +/- 0.72 | 7.2 +/- 1.9 |

The option 3 MCH-basis path has the best energy conservation in both gas phase
and QM/MM/PBC.  Option 1 is reasonably stable for this 100 fs formaldehyde
window because the active SOC states stay mostly dominated by one or a few MCH
components, so the omitted derivative-of-mixing terms are small or partially
cancel.  That behavior is system dependent and should not be treated as a proof
that the spin-adiabatic weighted-gradient approximation is generally exact.

## SOC/ISC hop bookkeeping

For SOC-QMMM hops, the code rescales QM velocities for the electronic energy gap
inside the FSSH hop kernel.  In QM/MM, an accepted ISC hop can also change the
ESPF electrostatic energy because the active-state density and fitted charges
change.  Both spin-adiabatic SOC-QMMM and MCH-basis SOC-QMMM therefore recompute
the post-hop embedded gradient and total potential, then apply an additional
isotropic full-system velocity rescaling for

```text
Delta E_ESPF = (E_pot_old - E_pot_new) + Delta E_electronic_gap
```

This keeps the hop energy bookkeeping consistent with the state-dependent ESPF
charges.  When `ESPF_ROHF=1` is set, the ESPF charges are state independent and
this correction should be zero or negligible.

## Overlap tracking and reproducibility

The NAMD driver stores the previous-step geometry and electronic data in memory
and sets `properties.back_door=True` so overlap/state-tracking calculations use
the previous trajectory step rather than recomputing an unrelated reference.
Random numbers come from `numpy.random.default_rng(seed)`, so repeated runs are
reproducible when the same input, executable, thread behavior, and seed are used.

The old rectangular AO-overlap OpenMP race can appear as `returncode=-11` in a
first attempt.  A completed production matrix should be counted by completed
trajectory directories/logs, not by raw first-attempt status lines.  The current
five-trajectory validation completed all 40 requested trajectories; three failed
first attempts were rerun to completion and excluded from the statistics.

## Limitations

- `soc_basis=mch` gives exact active MCH gradients, but hops occur in the
  spin-pure MCH basis rather than between SOC-diagonal spin-adiabatic states.
- Spin-adiabatic option 1 is SHARC-like but uses the weighted-gradient
  approximation.  It omits the off-diagonal derivative-Hamiltonian terms.
- Option 2 and `soc_tdc_grad_corr` are approximations from time-domain
  information.  They are useful diagnostics but are not recommended as the
  default production path.
- The Wigner initial-condition workflow is not exposed as a one-command NAMD
  setup in this checkout.  Current benchmark trajectories use Maxwell velocity
  seeds from the same equilibrated geometry.
- Longer production studies should still monitor total energy, hop acceptance,
  active-state labels, and population traces for every trajectory.

## Regression coverage

The focused regression tests in `tests/test_soc_namd_qmmm_production.py` check:

- `soc_basis=mch` dispatches SOC-QMMM jobs to `NAMD_SOC_MCH_QMMM`.
- the adiabatic SOC-QMMM path remains available for comparison.
- SOC-QMMM hops include the ESPF energy-change velocity correction.
- MCH-QMMM recomputes the exact active-root embedded gradient after accepted
  hops.
- overlap tracking and seeded RNG contracts remain present.
