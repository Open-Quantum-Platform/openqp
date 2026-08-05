# Port the irrep-coverage guard from `mrinivec` into `inivec` (TDA/RPA/SF)

## Why

PR #317 fixes a defect in the MRSF Davidson guess: a symmetry block that gets no
initial trial vector contributes **no roots at all** — absent, not unconverged —
so the survivors are renumbered and every state index silently shifts.

`inivec` (`source/tdhf_lib.F90:328-379`) has the same construction and **none of
the guard**: it is a plain energy-ordered seed, insertion sort into `itmp`, then
unit vectors, with no `pair_irrep` reference anywhere.

Callers:
- `source/modules/tdhf_energy.F90:199` — TDA/RPA
- `source/modules/tdhf_sf_energy.F90:286` — SF
- `source/modules/tdhf_mrsf_energy.F90:498` — quintet

## Where it bites hardest

**Plain TDA/RPA**, because `nvec = nstates` exactly (`tdhf_energy.F90:159`) with
no padding at all. An unseeded block is silently absent and every state index
below it shifts, propagating into excited-state gradients, optimisation,
MECP/MECI, the z-vector solve, NAMD and NMR.

Less severe for SF, where `nvec = min(max(2*nstates,5), mxvec)` already doubles
the guess.

## The table is already staged in the right layout

`stage_response_symmetry` accepts `td_type` in `tda/rpa/sf/mrsf`
(`pyoqp/oqp/molecule/molecule.py:684`), and `tdhf_energy.F90:199` calls with
`nocca = noccb = nocc`, giving `xvec_dim = nocc*(nbf-nocc)` — exactly the staged
`na*(nbf-na)`.

## Caveats to carry into the commit message

- Covering an irrep once guarantees only that block's **lowest** root. Two
  requested roots inside one block still need more vectors per block.
- The quintet path (`tdhf_mrsf_energy.F90:498`, `xvec_dim = noccb*nvira`) matches
  no staged table, so this is a no-op there unless `stage_response_symmetry`
  gains an `mrst==5` layout.
- `inivec` takes no `infos`, so all three call sites change.

## Verification

Re-run every TDA/RPA/SF example reference — the substitution is inert when all
irreps are already covered, so they must be **bit-identical**. Then CH2O 6-31G
nstate=3 against an nstate=20 reference, plus a C1-distorted CH2O as the
mechanism check (single block, must be unaffected).

Depends on #317 landing first (it carries the helper and the convention).
