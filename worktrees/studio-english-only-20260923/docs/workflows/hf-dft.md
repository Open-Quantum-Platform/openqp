# HF and DFT

HF and DFT use a basis-only or functional/basis `.oqp` route. The legacy
keyword form uses `[input] method=hf`; setting `[input] functional` selects
DFT, while leaving it empty gives Hartree-Fock.

## Energy

`.oqp`:

```text
hf/6-31g*
geom="h2o.xyz"
```

A single-point energy is the default task, so no `energy` keyword is needed; `hf/6-31g* energy` means the same thing.

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_hf", silent=1)
job.molecule(geometry="water", charge=0, multiplicity=1)
job.theory.hf(basis="6-31g*")

mol = job.run()
print("SCF energy:", mol.get_scf_energy())
```

Legacy `.inp`:

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

Runnable `.oqp`:
[`examples/HF/H2O_RHF-HF_ENERGY.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/HF/H2O_RHF-HF_ENERGY.oqp).
The same-stem `.inp` file is retained for legacy use.

## Gradient

For a DFT ground-state gradient, `grad` and `grad(S0)` are equivalent. In
Python, the same choice is written as `job.workflow.gradient(state=0)`.

`.oqp`:

```text
dft/bhhlyp/6-31g* grad
geom="h2o.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_dft_grad", silent=1)
job.molecule(geometry="water", charge=0, multiplicity=1)
job.theory.dft(functional="bhhlyp", basis="6-31g*")
job.workflow.gradient(state=0)

mol = job.run()
gradient = mol.get_grad()
```

Legacy `.inp`:

```ini
[input]
runtype=grad
method=hf
functional=bhhlyp
basis=6-31g*

[properties]
grad=0
```

Runnable `.oqp`:
[`examples/DFT/H2O_RHF-DFT_GRADIENT.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/DFT/H2O_RHF-DFT_GRADIENT.oqp).
The same-stem `.inp` file is retained for legacy use.

## Hessian

HF/DFT Hessians are run with `runtype=hess` and controlled by `[hess]`.
Analytical and numerical Hessian examples are collected on the
[Hessian and Frequencies](hessian.md) workflow page.
