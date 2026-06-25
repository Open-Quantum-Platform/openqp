# MRSF-TDDFT

MRSF-TDDFT uses `[input] method=tdhf`, `[tdhf] type=mrsf`, and usually a triplet
ROHF reference.

## Energy

```ini
[input]
runtype=energy
method=tdhf
functional=bhhlyp
basis=6-31g*

[scf]
type=rohf
multiplicity=3

[tdhf]
type=mrsf
nstate=3
```

Runnable reference:
`examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.inp`.

## Gradient

Use `runtype=grad` and select the response state through
`[properties] grad`.

```ini
[input]
runtype=grad
method=tdhf
functional=bhhlyp
basis=6-31g*

[scf]
type=rohf
multiplicity=3

[tdhf]
type=mrsf
nstate=3

[properties]
grad=3
```

Runnable reference:
`examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_GRADIENT.inp`.

## Notes

- `[tdhf] nstate` must include every state requested by gradients, NACME, SOC,
  or EKT analysis.
- For ordinary TDDFT, use `[tdhf] type=rpa` or `tda`.
- For spin-flip TDDFT without mixed-reference correction, use
  `[tdhf] type=sf`.
- UMRSF-TDDFT uses `[tdhf] type=umrsf` with a UHF reference and is currently an
  energy-only workflow.
