# The `input_to_standard` back-transform

Status: **scoped, not implemented.** The two bugs this draft originally bundled
with the feature are fixed in this PR; the transform itself remains open.

## What exists

`reorient_for_integral_symmetry` moves the molecule into the symmetry standard
frame and records the map that takes user input axes to it:

    r_std = (r_input - origin) @ rotation^T
    v_input = v_std @ rotation

It is written to `symmetry_metadata['integral_symmetry']['input_to_standard']`,
copied forward when the maps stage, and **read by nothing**. No gradient,
normal mode, dipole or velocity is transformed back.

## What this PR fixed instead

Two things were wrong today, and neither needed the transform.

**The geometry was left moved on every failure path.** `input_coords` was
captured in the reorientation loop and never used, so both failure exits
returned False with the molecule already rotated -- and *translated*, because
the standard frame is centred on the charge-weighted centroid. That moves every
molecule, including C1, where the rotation is exactly the identity and there is
no symmetry to exploit at all (measured: |R - I| = 0, |origin| = 2.8 bohr for a
CHFClBr-like geometry). Staging could also decline after the basis existed, for
reasons only visible then, leaving the same residue.

Both are now undone. The pre-reorient coordinates live in
`meta['_reorient_input_coords']`, deliberately outside the `integral_symmetry`
block -- every bail in `stage_integral_symmetry_maps` replaces that block
wholesale, so keeping them inside silently loses them. The caller pops the key
immediately, so it never reaches `save_data`.

**`prop` was on the allow-list.** It is the only admitted runtype that consumes
an *externally supplied* geometry expressed in the caller's frame:
`load_previous_data` writes `OQP::xyz_old` from it, and `get_basis_overlap`
overlaps that against the current -- by then rotated and translated -- structure.
The cross-geometry MO overlap, the phase/reorder alignment built on it, and
NACME are then all wrong. That is PR #319's finite-difference failure mode one
level up. The allow-list is now `('energy', 'grad')`. (`'properties'` was never
a `run_func` key, so it admitted nothing.)

## The remaining feature

Add `Molecule.to_input_frame()` implementing `v_input = v_std @ R` and
`r_input = r_std @ R + O`, then:

- Hook the vector back-transform into `Molecule.symmetrize_gradient`, which
  already sits on the gradient exit path and is already gated on
  `status == 'active'`. Apply the rotation **after** the symmetric projection --
  the projection is only valid in the standard frame.
- Restore the geometry at the job boundary, before the closing log dump.

### The constraint that makes this more than an inverse map

**The geometry cannot be restored before `save_data()` / `write_molden()`
without also transforming the MO coefficients.** The AO basis is atom-centred
and orientation-dependent: `OQP::VEC_MO_A/B` are expressed in the standard-frame
AO basis. Writing input-frame coordinates next to standard-frame MOs produces a
file that is internally inconsistent and silently wrong when read back as a
guess. Either transform both or neither.

### What this does and does not unlock

It **does** unlock the finite-difference case: a child that back-transforms to
its own input frame returns the gradient at the coordinates the parent handed
it. (Note PR #319 currently solves that by switching the reduction off in
children, which is the correct scope for a speed optimisation; the transform
would let them keep it.)

It does **not** by itself unlock the optimiser, IRC or NEB. Those displace the
geometry every step, so they need per-step re-detection and have to cope with
the point group *changing* between steps -- not just an inverse map. The draft
previously claimed otherwise.

## Prior art in this tree

`grd2.F90` produces a skeleton gradient that `symmetrize_gradient` projects, and
`scf_addons.F90` does the same for the skeleton XC matrix. The back-transform is
the same shape of operation applied one level further out, at the frame rather
than the group.
