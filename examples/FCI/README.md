# FCI examples (`method=fci`)

Small **dense** Full Configuration Interaction energy calculations on top of a
closed-shell RHF reference, using OpenQP-native molecular orbitals and AO
two-electron integrals.

> **MVP scope.** This is a correctness baseline and architecture seed, *not* a
> production-scale FCI engine. It builds an explicit determinant-basis
> Hamiltonian and diagonalizes it densely, so cost and memory grow
> exponentially with the active space. Use it for small systems only.

## How to run

```bash
openqp --nompi examples/FCI/H2_RHF-FCI_ENERGY.inp
```

A minimal input:

```ini
[input]
runtype=energy
method=fci
basis=sto-3g
system=
   H   0.0   0.0   0.0
   H   0.0   0.0   0.74

[scf]
type=rhf
multiplicity=1

[fci]
nroot=1
max_det=1000
```

## `[fci]` options

| key                | default  | meaning |
|--------------------|----------|---------|
| `nroot`            | `1`      | number of CI roots (states) to return |
| `active_electrons` | `0`      | electrons in the active space (`0` = all) |
| `active_orbitals`  | `0`      | spatial orbitals in the active space (`0` = all after the frozen core) |
| `frozen_core`      | `0`      | doubly-occupied orbitals frozen into the core |
| `max_det`          | `5000`   | hard cap on the determinant count |
| `max_memory`       | `2048`   | memory budget in MiB (dense Hamiltonian / AO integrals / Davidson vectors) |
| `eig_tol`          | `1.0e-10`| eigensolver residual tolerance (dense sanity check / Davidson convergence) |
| `integral_backend` | `native` | integral source (only `native` is supported) |
| `integral_cutoff`  | `5.0e-11`| magnitude below which integral contributions are skipped |
| `solver`           | `auto`   | `auto` (dense, or Davidson once dense exceeds `max_memory`), `dense`, or `davidson` |
| `davidson_maxiter` | `100`    | maximum Davidson iterations |

An active space is defined as `frozen_core` doubly-occupied core orbitals plus
`active_orbitals` active orbitals; any higher virtual orbitals are dropped
(CASCI-style). Setting `active_electrons` infers `frozen_core` automatically.

## Spin diagnostics

Each FCI root is labeled with `<S^2>` and an inferred spin multiplicity in the
log. The values are also stored in `mol.data` as `OQP::FCI_S2` (`float64`) and
`OQP::FCI_MULTIPLICITY` (`int64`). For multi-root jobs, OpenQP prints a note
when the requested roots span more than one multiplicity. As with any
non-spin-adapted alpha/beta determinant solver, exactly degenerate roots can mix
spin characters; treat non-integer `<S^2>` values in a degenerate manifold as a
mixing diagnostic rather than a clean spin label.

## Energy/properties scope

`method=fci` reports the variational electronic FCI/CASCI energy only. D4
(`input.d4=True`) is rejected because dispersion correction is not part of the
FCI solver. Any requested `[properties] scf_prop` entries are RHF
reference-density diagnostics, not FCI-density properties; set `scf_prop` empty
for a pure FCI energy job.

## Included example inputs

| input | purpose | main `[fci]` coverage |
|-------|---------|------------------------|
| `H2_RHF-FCI_ENERGY.inp` | minimal H2/STO-3G full FCI energy | default full active space, `nroot=1`, `max_det` |
| `H2_RHF-FCI_DAVIDSON_3ROOT.inp` | H2 three-root calculation with spin labels | `nroot=3`, `solver=davidson`, `max_memory`, `eig_tol`, `integral_cutoff`, `davidson_maxiter` |
| `H2O_RHF-FCI_CAS2E2O_FROZEN_CORE.inp` | H2O CASCI-style frozen-core active-space calculation | `active_electrons`, `active_orbitals`, `frozen_core`, `solver=dense`, `max_memory` |

Each `.inp` file has a matching `.json` reference so `openqp --nompi --run_tests FCI`
checks the FCI examples through OpenQP's normal example-test harness.

## Validation

The dense solver has been checked against PySCF (FCI / CASCI) with **identical
MO integrals**, where it reproduces the reference energy to machine precision:

| system / space                  | PySCF reference (Eh) | dets | |Δ| vs PySCF |
|---------------------------------|----------------------|------|--------------|
| H2 / STO-3G, full FCI           | −1.137283834489      | 4    | ~9e-16 |
| H2O / STO-3G, full FCI          | −75.012647118992     | 441  | ~1e-12 |
| HeH⁺ / STO-3G, full FCI         | −2.826674836464      | —    | ~4e-16 |
| H2O / STO-3G, CAS(8e,6o)        | −75.012569053757     | 225  | ~1e-14 |
| H2O / STO-3G, CAS(2e,2o)        | −74.964311401165     | 4    | ~3e-14 |
| H2, three lowest roots          | −1.1373/−0.5308/−0.1684 | — | ~1e-16 |

End to end (OpenQP's own RHF orbitals and AO integrals), `H2_RHF-FCI_ENERGY.inp`
yields −1.137283835869487 Eh (stored in `H2_RHF-FCI_ENERGY.json`). The small
difference from the PySCF total above reflects SCF/integral thresholds, not the
FCI algorithm. These checks are reproduced by `tests/test_fci.py`.

## Current limitations

1. Closed-shell RHF singlet **energy** only (no gradients, no open shell, no UHF/ROHF).
2. Two solvers: a dense `numpy.linalg.eigh` path and a sparse-Hamiltonian block
   **Davidson** path (`solver=davidson`/`auto`) that avoids the `ndet**2` dense
   matrix. The matrix-element enumeration is still pure Python, so it is the
   current speed bottleneck — a fully matrix-free / string-based (Olsen) sigma
   build is the next scaling step.
3. Dense `nbf**4` AO ERI storage (no symmetry compression); guarded by `max_memory`.
4. Practical only for small active spaces (Davidson relaxes the memory wall but
   not the pure-Python build cost).
