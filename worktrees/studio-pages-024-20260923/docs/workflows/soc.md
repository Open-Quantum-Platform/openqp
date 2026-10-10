# Spin-Orbit Coupling

OpenQP supports MRSF-TDDFT spin-orbit coupling through `runtype=soc`.
SOC mixes spin-free electronic states through relativistic one-electron terms
and, when requested, a mean-field treatment of the two-electron SOC operator.
This workflow is useful when intersystem crossing, phosphorescence, or
singlet-triplet mixing is important. See the
[References](../references.md#spin-orbit-coupling) page for OpenQP's
relativistic MRSF-TDDFT SOC method and the mean-field SOC operator background.

## `.oqp`

Use one route count for equal singlet and triplet spaces:

```text
mrsf(nstate=3)/bhhlyp/6-31g* soc
geom="h2o.xyz"
```

This requests `S0`--`S2` and `T0`--`T2`. For unequal spaces, omit route
`nstate` and specify both counts:

```text
mrsf/bhhlyp/6-31g* soc(ns=3,nt=5)
geom="h2o.xyz"
```

Do not combine the two forms. The canonical route also recognizes MRSF-TDHF;
the detailed all-electron example below is MRSF-TDDFT.

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_soc", silent=1)
job.molecule(geometry="water", charge=0)
job.theory.mrsf(functional="bhhlyp", basis="6-31G(2df,p)", nstate=12)
job.workflow.soc(soc_2e=1, scal_rel=2)

mol = job.run()
soc = mol.get_soc()
```

Legacy `.inp`:

```ini
[input]
runtype=soc
method=tdhf
functional=bhhlyp
basis=6-31G(2df,p)
soc_2e=1

[scf]
type=rohf
multiplicity=3
scal_rel=2

[tdhf]
type=mrsf
nstate=12
```

Runnable `.oqp` inputs:

- [`examples/SOC/H2O_BHHLYP_SOC.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/SOC/H2O_BHHLYP_SOC.oqp)
- [`examples/SOC/CH3Br-BHHLYP-SOC.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/SOC/CH3Br-BHHLYP-SOC.oqp)

Each has a same-stem legacy `.inp` companion.

## SOC Terms

`[input] soc_2e` controls the SOC Hamiltonian. The current MRSF-SOC driver forms
one-electron Breit-Pauli SOC matrix elements and can add a mean-field
two-electron SOC contribution:

| Value | Meaning |
| --- | --- |
| `0` | One-electron SOC terms only. |
| `1` | One-electron plus mean-field two-electron SOC terms. |

The all-electron workflow on this page uses MRSF with its automatic triplet
ROHF working reference and `[tdhf] type=mrsf` in legacy input.
The SOC driver computes both singlet and triplet response roots internally, so
Python scripts should use `job.workflow.soc(...)` after
`job.theory.mrsf(...)` rather than setting
`job.tdhf.multiplicity`.
The Python helper sets `[scf] scal_rel=2` by default for the SOC workflow; pass
`scal_rel=0`, `1`, or `2` to override it.

## Scalar Relativistic Correction

Scalar relativistic correction is separate from SOC. It changes the spin-free
one-electron core Hamiltonian and does not mix spin states by itself. In the
current OpenQP input surface, the exposed scalar-relativistic keyword is
`[scf] scal_rel`:

| Value | Meaning |
| --- | --- |
| `0` | No scalar relativistic correction. |
| `1` | First-order Douglas-Kroll-Hess correction. |
| `2` | First- and second-order Douglas-Kroll-Hess correction. |

Use `scal_rel` when the spin-free Hamiltonian should include DKH scalar
relativistic effects. Use `runtype=soc` and `soc_2e` when the calculation needs
spin-orbit coupling. X2C support is under implementation and should be
documented as a user-facing option only after an OpenQP keyword is exposed.
