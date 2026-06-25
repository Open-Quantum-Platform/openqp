# HF and DFT

HF and DFT calculations use `[input] method=hf`. Set `[input] functional` only
for DFT; leaving it empty gives Hartree-Fock.

## Energy

```ini
[input]
system=
   O   0.000000000   0.000000000  -0.041061554
   H  -0.533194329   0.533194329  -0.614469223
   H   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g*
method=hf

[guess]
type=huckel

[scf]
type=rhf
multiplicity=1
```

Runnable reference: `examples/HF/H2O_RHF-HF_ENERGY.inp`.

## Gradient

For a DFT gradient, request `runtype=grad`, set a functional, and use
`[properties] grad=0` for the reference-state gradient:

```ini
[input]
runtype=grad
method=hf
functional=bhhlyp
basis=6-31g*

[properties]
grad=0
```

Runnable reference: `examples/DFT/H2O_RHF-DFT_GRADIENT.inp`.

## Hessian

HF/DFT Hessians are run with `runtype=hess` and controlled by `[hess]`.

Runnable reference: `examples/HESS/H2O_RHF-DFT_ANA_HESS.inp`.
