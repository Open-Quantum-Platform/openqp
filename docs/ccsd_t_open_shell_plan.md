# Open-shell CCSD(T): implementation plan

Status: **not implemented**. The solver currently accepts an RHF reference
only (`source/modules/ccsd_t_energy.F90`, the `scftype /= 1` guard). This
document is the design for lifting that, written against the code as it
stands so the work can start without re-deriving any of it.

## Why not spin-orbital

The quickest correct route is a spin-orbital code: one amplitude set over
antisymmetrised `<pq||rs>`, and the closed-shell equations collapse to a
much shorter list. It is the wrong choice here on memory alone.

Spin-orbital dimensions are `2*nmo`, so the ladder integrals go as
`(2 nv)^4 = 16 nv^4`:

| nv  | spatial `nv^4` | spin-orbital `(2nv)^4` |
|-----|----------------|------------------------|
| 40  | 0.02 GB        | 0.31 GB                |
| 88  | 0.45 GB        | 7.2 GB                 |
| 150 | 3.8 GB         | 60.7 GB                |

That gives back the 10.6x peak reduction the packed AO store just bought,
and then some. Use spin-integrated blocks.

## Spin-integrated structure

Three amplitude sets and three integral spin-cases:

- `t1a(i,a)`, `t1b(I,A)`
- `t2aa(i,j,a,b)`, `t2bb(I,J,A,B)`, `t2ab(i,J,a,B)`
- integrals `(pq|rs)` in the `aa`, `bb`, `ab` spin-cases

Same-spin blocks carry antisymmetrised integrals `<ij||ab>`; the opposite-spin
block does not antisymmetrise. `t2ab` has no permutational symmetry between
its two index pairs, so it is a full `no_a x no_b x nv_a x nv_b` array.

Equations: Gauss, Stanton, Bartlett, *J. Chem. Phys.* **95**, 2623 (1991) for
UCCSD; Watts, Gauss, Bartlett, *J. Chem. Phys.* **98**, 8718 (1993) for the
open-shell (T). Use the UCCSD(T) form, not ROHF-CCSD(T)-with-RHF-equations.

## Reference handling

A ROHF reference has a non-diagonal Fock matrix in the occ-occ and vir-vir
blocks, so the canonical denominators are undefined until it is
semicanonicalised. That code already exists and is directly reusable:

    source/mp2_lib.F90 :: semicanonicalize(nbf, nocc, cmo, fock_packed, cmo_sc, e_sc)

It diagonalises the occ-occ and vir-vir Fock sub-blocks per spin and returns
rotated orbitals plus their energies. `mp2_correlation` shows the exact
calling pattern for both spins, including the RHF case where the beta
pointers alias the alpha ones. For UHF it is a no-op; for ROHF it produces
the standard semicanonical basis.

Note that after semicanonicalisation the occ-vir Fock block is generally
**non-zero** for ROHF. The closed-shell code drops the `f_ov` terms from the
singles equation and the disconnected triples because they vanish for a
canonical RHF reference; the open-shell code must put them back.

## Integration points

1. `source/modules/ccsd_t_energy.F90`
   - Replace the `scftype /= 1` abort with a branch. Keep the closed-shell
     path exactly as it is — it is validated and faster; only route
     `scftype == 2` (UHF) and `3` (ROHF) to the new solver.
   - `nelec_a` / `nelec_b` are on `infos%mol_prop`; `nfzc` applies per spin.
   - Read `OQP_VEC_MO_B` / `OQP_E_MO_B` / `OQP_FOCK_B` alongside the alpha
     records, aliasing to alpha when `scftype == 1`.

2. `source/cc_ao2mo.F90`
   - The packed AO store is spin-free and needs no change.
   - `cc_build_mo_blocks` currently takes one `cmo` and emits six blocks.
     Generalise to take `(cmo_a, cmo_b)` and emit the `aa`, `bb`, `ab`
     cases. The two half-transformations stay as they are; only the second
     half's `scatter_pair` needs a spin-case argument.

3. `source/cc_lib.F90`
   - Add `uccsd_iterate` and `utriples_correction` beside the existing
     routines. Do not try to make one routine serve both spin cases — the
     closed-shell equations are a genuinely different (and cheaper)
     factorisation, not a special case to be branched into.
   - The DIIS, the ladder blocking and the `(T)` triple decomposition all
     carry over unchanged in structure.

## Parallelisation

Reuse what is already there and already measured:

- `(T)` over `a >= b >= c` virtual triples, MPI round-robin + OpenMP dynamic.
  Open-shell has three triples spin-cases (`aaa`, `aab`, `abb`, `bbb`);
  parallelise within each, not across them, so the load stays balanced.
- The ladder blocked over the last virtual index, with the block capped by
  rank and thread count (`BLOCKS_PER_RANK`).
- BLAS pinned to one thread for the whole solver. This matters: leaving it
  threaded under the solver's own OpenMP measured ~2x slower on four cores.

## Validation

Extend `tests/test_ccsd_t_reference.py`, which already drives the native
path and compares against `tests/data/ccsd_t_pyscf_validation.json`.

Generate the references from PySCF `UCCSD` / `uccsd_t` — do not transcribe
them by hand. Suggested first cases, all small enough to iterate on:

| system | basis | mult | reference |
|--------|-------|------|-----------|
| OH     | cc-pVDZ | 2 | UHF |
| CH2    | cc-pVDZ | 3 | UHF and ROHF |
| O2     | cc-pVDZ | 3 | UHF |
| NH2    | 6-31G   | 2 | ROHF |

Two traps the closed-shell validation walked into, both worth avoiding here:

- **Check the SCF energy and `nbf` before comparing any correlation energy.**
  A correlation energy is defined relative to its reference determinant; a
  mismatched reference produces a CC difference that looks like a CC bug.
  The existing test asserts both first, and says so when it fails.
- **OpenQP uses 6d Cartesian shells for the Pople sets.** PySCF needs
  `gto.M(..., cart=True)` to be compared against them, or the two codes are
  running different basis sets (28 vs 30 functions for CO/6-31G*).

A UHF reference on a closed-shell system must reproduce the RHF-CCSD(T)
number exactly — that is the cheapest first regression, and it catches most
spin-case bookkeeping errors before any genuinely open-shell case is run.

## Sizing

Roughly 3x the closed-shell equation set. The closed-shell solver is 1126
lines in `cc_lib.F90`. Budget accordingly, and validate incrementally:
UCCSD singles+doubles against PySCF `UCCSD` first, with `(T)` switched off,
before starting the triples.
