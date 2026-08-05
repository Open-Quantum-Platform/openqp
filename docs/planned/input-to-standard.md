# Implement the `input_to_standard` back-transform

## The gap

`input_to_standard` is written at `pyoqp/oqp/molecule/molecule.py:489` and copied
at `:641` and `:649`, and **consumed nowhere**. A grep over `pyoqp/`, `source/`
and `tests/` returns only those three sites.

Nothing back-transforms a gradient, a normal mode, or a velocity out of the
standard frame.

## Why it matters

This single missing piece is the concrete technical reason the reorientation
runtype allow-list exists (`molecule.py:445-447`), which keeps the optimizer,
Hessian, IRC, NEB and NAMD in the input frame and therefore **excludes them from
every integral-symmetry benefit**.

It is the prerequisite for widening `use_integral_symmetry` past a single-point
energy.

## Also fix while here

The non-converged branch at `:479` returns `False` **without restoring the
geometry** the loop already mutated through `update_system()`.

## Restart hazard to document

Reorientation is idempotent — the loop at `:456-472` iterates to the identity
frame — so restarting from a standard-frame `.oqp` is a no-op. But the stored
`input_to_standard` then becomes the identity, and the mapping back to the user's
original axes is **permanently lost**. Decide and document whether the original
frame is preserved across restarts.
