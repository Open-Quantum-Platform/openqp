# NACME

NACME calculations use two geometries and an MRSF-TDDFT response calculation.
The current input style uses `[input] system` and `[input] system2`.

```ini
[input]
runtype=nacme
method=tdhf
functional=bhhlyp
basis=6-31g
system=
  O   0.000000000   0.000000000  -0.041061554
  H  -0.533194329   0.533194329  -0.614469223
  H   0.533194329  -0.533194329  -0.614469223
system2=
  O   0.000000000   0.000000000  -0.031061554
  H  -0.543194329   0.543194329  -0.624469223
  H   0.543194329  -0.543194329  -0.624469223

[scf]
type=rohf
multiplicity=3

[tdhf]
type=mrsf
nstate=10
```

Runnable reference:
`examples/other/h2o_nacme_rohf_mrsf-s_6-31g_bhhlyp.inp`.

## Notes

- Keep `[tdhf] nstate` large enough to cover all states used in the coupling
  analysis.
- Use consistent atom ordering between `system` and `system2`.
- The `[nac]` section controls finite-difference and restart behavior for NAC
  workflows.
