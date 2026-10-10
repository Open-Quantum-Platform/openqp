# Hessian and Frequencies

Hessian workflows use `[input] runtype=hess` and are controlled by `[hess]`.
OpenQP supports native analytical Hessians for supported HF/DFT ground-state
cases and numerical Hessians for broader workflows.

Frequency, IR, Raman, and thermochemistry analysis are built from Hessian data
when the selected workflow produces the required derivatives.

The generated `.hess.json` and `.freq.molden` files contain portable normal
modes and, for supported bases, the SCF molecular orbitals needed to view MO
surfaces and vibrational animation together. See
[Orbital and Vibrational Output](orbital-vibrational-output.md).

In Python scripts, start from a compact HF, DFT, or MRSF-TDDFT theory setup and
then select the Hessian workflow with `job.workflow.hessian(...)`.

## Analytical HF/DFT Hessian

`.oqp`:

```text
dft/bhhlyp/6-31g* hess(S0,type=analytical,clean=true)
geom="h2o.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_dft_hess", silent=1)
job.molecule(geometry="water", charge=0, multiplicity=1)
job.theory.dft(functional="bhhlyp", basis="6-31g*")
job.workflow.hessian(type="analytical", state=0, clean=True)

mol = job.run()
hessian = mol.get_hess()
```

Legacy `.inp`:

```ini
[input]
runtype=hess
method=hf
functional=bhhlyp
basis=6-31g*

[scf]
type=rhf
multiplicity=1

[hess]
type=analytical
state=0
clean=True
```

Runnable `.oqp`:
[`examples/HESS/H2O_RHF-DFT_ANA_HESS.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/HESS/H2O_RHF-DFT_ANA_HESS.oqp).
The same-stem `.inp` file is retained for legacy use.

## Numerical HF/DFT Hessian

Omit `type=analytical` to use the numerical finite-difference path.

`.oqp`:

```text
dft/bhhlyp/6-31g* hess(S0,clean=true)
geom="h2o.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_dft_num_hess", silent=1)
job.molecule(geometry="water", charge=0, multiplicity=1)
job.theory.dft(functional="bhhlyp", basis="6-31g*")
job.workflow.hessian(state=0, clean=True)
mol = job.run()
```

Legacy `.inp`:

```ini
[input]
runtype=hess
method=hf
functional=bhhlyp
basis=6-31g*

[scf]
type=rhf
multiplicity=1

[hess]
state=0
clean=True
```

For a symmetric molecule, `symmetry_unique=True` can reduce the numerical
finite-difference work to one displaced atom per symmetry-equivalent atom
orbit. It is an opt-in feature and falls back to the full displacement set
with an explanatory note whenever the current geometry or tolerance does not
support a complete symmetry reconstruction.

Runnable `.oqp`:
[`examples/HESS/H2O_RHF-DFT_NUM_HESS.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/HESS/H2O_RHF-DFT_NUM_HESS.oqp).
The same-stem `.inp` file is retained for legacy use.

## Numerical MRSF-TDDFT Hessian

MRSF-TDDFT Hessians use the numerical path in the documented examples. The
MRSF state numbering follows the MRSF target-state list; `state=1` is the
lowest MRSF target state, which can be the multiconfigurational ground state.

`.oqp`:

```text
mrsf(nstate=2)/bhhlyp/6-31g* hess(S0,clean=true)
geom="h2o.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_mrsf_hess", silent=1)
job.molecule(geometry="water", charge=0)
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=2)
job.workflow.hessian(state=1, clean=True)

mol = job.run()
```

Legacy `.inp`:

```ini
[input]
runtype=hess
method=tdhf
functional=bhhlyp
basis=6-31g*

[scf]
type=rohf
multiplicity=3

[tdhf]
type=mrsf
nstate=2

[hess]
state=1
clean=True
```

Runnable `.oqp`:
[`examples/HESS/H2O_BHHLYP-MRSFTDDFT_NUM_HESS.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/HESS/H2O_BHHLYP-MRSFTDDFT_NUM_HESS.oqp).
The same-stem `.inp` file is retained for legacy use.

## Notes

- HF/DFT ground-state Hessians use `state=0`.
- TDHF/TDDFT Hessian target states use positive excited-state indices.
- SF-TDDFT and MRSF-TDDFT use spin-flip/MRSF target-state ordering, where
  state `1` is the lowest target state.
- `[hess] restart=True` can continue a numerical Hessian workflow where the
  corresponding temporary files are available.
- `[hess] read=True` accepts only a current versioned `.hess.json` whose
  electronic-model configuration, state, atoms, geometry, and isotopic masses
  match the present job. Referenced file contents and the OpenQP binary are not
  content-fingerprinted, so regenerate the sidecar after changing either one.
  Unsigned legacy sidecars must be regenerated with `read=False`.
- `[hess] clean=True` removes temporary Hessian files where supported.
