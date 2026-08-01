# CASSCF examples (`method=casscf`)

This directory keeps the original zero-orbital-update CASSCF scaffold example.
Current native one-macroiteration CASSCF and SA-CASSCF smoke examples live in
`examples/WaveFunction`.

`[casscf] max_macro_iterations=0` performs zero orbital-update macro iterations
and dispatches to the same fixed-orbital native CASCI determinant-CI path used by
`examples/WaveFunction/H2_RHF-CASCI_ENERGY.inp`.

## How to run

```bash
openqp --nompi examples/CASSCF/H2_RHF-CASSCF_ZERO_MACRO_CASCI_EQUIVALENT.inp
```

The matching JSON reference stores the CASCI/FCI-equivalent H2/STO-3G
CAS(2e,2o) energy:

```json
{
  "energy": -1.137283835869487
}
```

## Notes

- This is a scaffold/API smoke test for the explicit zero-orbital-update path.
- It is CASCI-equivalent because the RHF orbitals are held fixed.
- Use `examples/WaveFunction/LiH_RHF-CASSCF_ONE_MACROITERATION.inp` for the
  native state-specific macroiteration smoke example.
- Use `examples/WaveFunction/LiH_RHF-SA-CASSCF_ONE_MACROITERATION.inp` for the
  native equal-weight two-state SA-CASSCF macroiteration smoke example.
