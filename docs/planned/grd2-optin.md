# Give `grd2` the same per-caller petite opt-in that `int2` has

## Precise statement

`load_petite_shell_map` (`source/integrals/int2.F90:334-372`) does check the
global `OQP::sym_petite` flag, so `grd2` is off whenever symmetry is off. The
difference is that **`int2` additionally requires a per-driver-instance opt-in**
(`enable_petite`, used only by `fock_jk(..., petite=.true.)`), which is what
keeps response and CPHF builds out of the reduction.

`source/integrals/grd2.F90:253` has no such gate, so *any* `grd2_driver` call
inherits the reduction.

## Why that is a problem

`fock_deriv_contract` / `_os` (`source/modules/fock_deriv.F90:114`, `:266`) route
CPHF probe densities through it, and their only callers are in
`source/modules/hf_hessian.F90` (:179, :353, :360, :721, :914, :931, :1268,
:1280). Those perturbed densities `dPx` are **not totally symmetric**, so the
petite weighting is invalid for them.

## Status

**Latent today**: `hess` fails the reorientation runtype allow-list and
`stage_integral_symmetry_maps` requires `status == 'reoriented'`
(`pyoqp/oqp/molecule/molecule.py:511`), so the maps are never staged for it.

It becomes **live** the moment that gate is widened — which is exactly what the
integral-symmetry work would do.

Not verified yet: whether every `grd2` contraction reachable under
`runtype='prop'` uses a totally symmetric density. That needs checking too.

## What to do

Give `grd2` an explicit per-caller opt-in mirroring `int2`'s `enable_petite`.

This is **not a speedup** — it is the audit that must precede every other
integral-symmetry item.
