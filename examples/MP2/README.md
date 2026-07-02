# MP2 ground-state correlation

Standalone second-order Møller–Plesset (MP2) correlation energy for RHF/UHF/ROHF
references, invoked with `method = mp2` in the `[input]` block. The SCF reference
is converged first (as for `method = hf`); the MP2 correlation is then added on
top and the total `E(HF+MP2)` is reported.

## Method

- The correlation energy is assembled in the spin-blocked form
  `E_MP2 = E_aa + E_bb + E_ab` on **semicanonicalized** orbitals (the occ–occ and
  vir–vir Fock blocks are diagonalized first, so the canonical amplitude
  denominators are well defined even for a ROHF reference).
- The two-electron integrals are contracted through the validated `int2_compute`
  driver via per-occupied-MO-pair Coulomb builds, so **no O(N⁴) MO-integral
  tensor is ever stored** — memory scales as O(N²).

## Example

`h2o_ump2_6-31g.inp` — H₂O / 6-31G, UHF reference (closed shell ⇒ UMP2 = RMP2).

Reference (PySCF `mp.UMP2`, see `h2o_ump2_6-31g.ref.txt`):

| quantity            | value (Ha)      |
|---------------------|-----------------|
| E(UHF)              | −75.9842900384  |
| E(MP2 correlation)  |  −0.1278307500  |
| E(UHF+UMP2 total)   | −76.1121207884  |

Run:

```bash
python3 -m oqp.pyoqp h2o_ump2_6-31g.inp
grep -A6 'Moller-Plesset' h2o_ump2_6-31g.log
```
