# MP2 ground-state correlation

Standalone second-order Møller–Plesset (MP2) correlation energy for RHF/UHF/ROHF
references, invoked with `method = mp2` in the `[input]` block. The SCF reference
is converged first (as for `method = hf`); the MP2 correlation is then added on
top and the total `E(HF+MP2)` is reported.

Optional spin-scaled variants are selected in `[mp2]`:

```ini
[mp2]
variant=scs-mp2
```

Supported presets are `mp2`, `scs-mp2`, `sos-mp2`, `os-mp2`, `ss-mp2`, and
`scs-mi-mp2`. Use `variant=custom` with `same_spin_scale` and
`opposite_spin_scale` for other literature parameterizations.

## Method

- The correlation energy is assembled in the spin-blocked form
  `E_MP2 = E_aa + E_bb + E_ab` on **semicanonicalized** orbitals (the occ–occ and
  vir–vir Fock blocks are diagonalized first, so the canonical amplitude
  denominators are well defined even for a ROHF reference).
- `E_aa + E_bb` and `E_ab` can be scaled independently, so OS-only, SS-only,
  SCS-MP2, SOS-MP2, and custom spin-component-scaled totals reuse the same
  kernel.
- The two-electron integrals are contracted through the validated `int2_compute`
  driver via per-occupied-MO-pair Coulomb builds, so **no O(N⁴) MO-integral
  tensor is ever stored**. The same-spin contraction keeps only one occupied
  block of virtual-virtual pair integrals at a time.

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
