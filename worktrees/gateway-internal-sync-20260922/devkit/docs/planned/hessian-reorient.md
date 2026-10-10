# Numerical-Hessian and IR/Raman child jobs silently reorient

## The defect

`pyoqp/oqp/library/single_point.py:1499` (Hessian) and `:1163` (property /
IR-Raman) spawn child jobs with `copy.deepcopy(self.mol.config)` — including the
`[symmetry]` section — and set `runtype='grad'` and `runtype='energy'`
respectively. Both of those pass the reorientation allow-list in
`molecule.py:445-447`, while the **parent** is `runtype='hess'` and is excluded.

So each displaced child re-detects its own point group and rotates into its own
standard frame, then returns a gradient (or dipole / polarizability) expressed in
that frame, while the parent assembles them in the input frame at `:1441-1445`
with no inverse rotation.

## Why it is nasty

Fully symmetry-broken (C1) displacements get the identity rotation and are fine.
The corruption therefore hits **exactly the displacements that retain a symmetry
element** — the ones carrying the most Hessian signal — and is intermittent
rather than uniform.

## Scope

Only runs that explicitly set the experimental `use_integral_symmetry`. Those get
wrong frequencies, wrong thermochemistry, wrong IR/Raman intensities, wrong
`oqp.init_hessian=numerical` TS searches, and a wrong TS Hessian in the IRC
(`liboqp.py:1449`).

## Fix

Force `use_integral_symmetry=False` in the child config at both sites.

This is strictly better than back-rotating the returned vectors: back-rotation
would fix the frame but leave the child seeded with a parent-frame MO guess
against a rotated geometry.

## Verification

One regression test: a C2v numerical Hessian with the flag true and false. Today
they differ; after the fix they must agree.
