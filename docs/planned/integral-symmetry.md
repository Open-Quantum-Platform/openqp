# Petite-list integral reduction: measured, and currently worth nothing

## What was measured

Benzene, D2h, cc-pVTZ (264 basis functions), RHF, 8 OpenMP threads, same machine,
back to back:

| run | wall clock | energy |
|---|---|---|
| `[symmetry] enabled=false` | 40.1 s | -230.7789916457 |
| `enabled=true, use_integral_symmetry=true` | 40.5 s | -230.7789916457 |

Energies agree to all printed digits, so whatever runs is not wrong. But there is
**no speed benefit at all** — the symmetric run is marginally slower. The same
comparison at 6-31G* also gave identical energies and no measurable difference.

## Why this is surprising

The reduction is implemented. `source/integrals/int2.F90:724` computes
`q4 = petite_quartet_weight(...)` per shell quartet and `cycle`s when `q4 == 0`,
weighting the surviving orbit representative. For D2h (order 8) that should be a
large saving in the 2e build.

The point group is detected correctly (`d2h` in the log) and the flag is on
(`use integral symmetry: yes`).

## What to find out first

1. Does the petite path actually engage at run time, or does it fail-safe to C1?
   `stage_integral_symmetry_maps` (`pyoqp/oqp/molecule/molecule.py:499`) returns
   early unless `meta['integral_symmetry']['status'] == 'reoriented'` (:511), and
   `reorient_for_integral_symmetry` has its own gates. If the maps are never
   staged, `OQP::sym_petite` stays 0 and `int2.F90:306` returns before doing
   anything. Instrument `this%petite` and the skipped-quartet count and compare.
2. If it does engage, where does the saving go? Schwarz screening may already be
   removing most of what the petite list would remove, in which case the honest
   answer is that the reduction is redundant with screening for this class of
   system and the gain only appears elsewhere (tight thresholds, diffuse sets).
3. Only then decide about defaults.

## Blockers to enabling it by default

- `use_integral_symmetry` **reorients the molecule to the standard frame at load
  time**, so anything that persists or compares geometries (restarts, saved
  `.oqp`, gradients, NAMD trajectories, optimisation histories, reference data)
  has to be checked first.
- `input_to_standard` is written (`molecule.py:489`) and copied (:641, :649) but
  **consumed nowhere** — nothing back-transforms a gradient, mode or velocity out
  of the standard frame. That missing piece is why the runtype allow-list keeps
  the optimizer and Hessian in the input frame today.
- `source/integrals/grd2.F90:253` has no per-caller opt-in, unlike `int2.F90:291`
  (`enable_petite`). CPHF probe densities routed through
  `fock_deriv_contract` are not totally symmetric. Latent today because `hess`
  fails the runtype gate; live the moment that gate is widened.

## Suggested order

Instrument first (is it even on?), then fix whatever prevents the saving, then
the back-transform, then the `grd2` opt-in audit, and only then discuss defaults.
