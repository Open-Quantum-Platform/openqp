# Duschinsky transformation Python API

OpenQP provides a harmonic-coordinate transformation between two stationary
geometries.  For mass-weighted normal coordinates, the returned quantities
satisfy

\[
Q_{\mathrm{exc}} = J Q_{\mathrm{gnd}} + K .
\]

The implementation removes translations and rotations independently at the
two geometries, aligns the excited-state geometry to the ground-state frame by
a mass-weighted proper rotation, and then constructs the transformation in the
physical vibrational spaces.

## Two OpenQP `Molecule` objects

```python
from oqp.library import compute_duschinsky_from_molecules

result = compute_duschinsky_from_molecules(ground_mol, excited_mol)

J = result.J
K = result.K
```

The helper reads coordinates, masses, atomic identities, and Cartesian
Hessians through `get_system()`, `get_mass()`, `get_atoms()`, and `get_hess()`.
Both molecules must contain Hessians.  Either stored Hessian can instead be
supplied explicitly:

```python
result = compute_duschinsky_from_molecules(
    ground_mol,
    excited_mol,
    ground_hessian=ground_hessian,
    excited_hessian=excited_hessian,
)
```

The two molecules must have identical atom order and isotopic masses.  A
mismatch raises `ValueError` before normal-mode analysis.

## Explicit arrays

```python
from oqp.library import compute_duschinsky_from_arrays

result = compute_duschinsky_from_arrays(
    ground_geometry,       # shape (N, 3)
    ground_masses,         # shape (N,)
    ground_hessian,        # shape (3N, 3N)
    excited_geometry,      # shape (N, 3)
    excited_masses,        # shape (N,)
    excited_hessian,       # shape (3N, 3N)
    ground_atoms=ground_atomic_numbers,
    excited_atoms=excited_atomic_numbers,
)
```

`ground_atoms` and `excited_atoms` are optional, but they must be supplied
together.  Providing them enables an exact check of atomic identity and order;
the two mass arrays are always compared.  The coordinates and Hessians must use
one consistent unit system.  If length is in bohr and mass is in atomic mass
units, `K` is in `sqrt(atomic mass unit) * bohr`.

## Physical and numerical conventions

- A nonlinear molecule has `3N - 6` physical modes and a linear molecule has
  `3N - 5`.  A transformation between a linear and a nonlinear geometry is
  rejected because the mode counts differ.
- Negative projected force constants are retained.  The API does not replace
  them with positive values or convert them into signed frequencies.
- Mode phases are aligned reproducibly.  A rotation within a degenerate
  subspace is applied only when the same contiguous mode block is degenerate
  at both geometries; nondegenerate Duschinsky mixing is preserved.
- `result.orthogonality` reports singular values and two-sided residuals of
  `J`.  Because translations and rotations are projected at two different
  geometries, their physical vibrational subspaces need not coincide exactly.
  The implementation therefore reports a nonorthogonal `J` when appropriate
  and does not replace it by a nearest orthogonal matrix.
- `orthogonality_tolerance` can be supplied to either public helper when an
  application requires a fail-closed residual threshold.

## Scope boundary

The API returns normal modes, the Duschinsky matrix `J`, displacement `K`, the
rigid alignment, and numerical diagnostics.  It does **not** calculate
Franck--Condon factors, Herzberg--Teller contributions, transition-moment
derivatives, thermal populations, line shapes, or spectral intensities.  Those
quantities require additional electronic and vibrational data and a separate
spectral model.
