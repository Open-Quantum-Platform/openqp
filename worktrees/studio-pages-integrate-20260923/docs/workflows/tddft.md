# TDDFT and TDHF

TDDFT and TDHF response calculations use `[input] method=tdhf`. A non-empty
`[input] functional` selects TDDFT; leaving `functional` empty gives TDHF.
Closed-shell examples normally use an RHF reference with `[scf] multiplicity=1`.

Use MRSF-TDDFT when multiconfigurational ground-state character, conical
intersections, or spin-flip target-state ordering are central to the study. Use
ordinary TDDFT/TDHF for standard single-reference vertical excitations and
excited-state gradients.

In Python scripts, TDHF and TDDFT are theory choices. Select the calculation
type separately with `job.workflow.*` only when you need a non-energy workflow.

## Energy

`.oqp`:

```text
tddft(nstate=3)/b3lyp5/6-31g*
geom="h2o.xyz"
```

For TDHF, use the same layout with the route `tdhf(nstate=3)/6-31g*`.

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_tddft", silent=1)
job.molecule(geometry="water", charge=0, multiplicity=1)
job.theory.tddft(functional="b3lyp5", basis="6-31g*", nstate=3)

mol = job.run()
print("TD energies:", mol.get_td_energies())
```

Legacy `.inp`:

```ini
[input]
runtype=energy
method=tdhf
functional=b3lyp5
basis=6-31g*

[scf]
type=rhf
multiplicity=1

[tdhf]
nstate=3
```

Runnable `.oqp`:
[`examples/TDDFT/H2O_B3LYP5-TDDFT_ENERGY.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/TDDFT/H2O_B3LYP5-TDDFT_ENERGY.oqp).
The same-stem `.inp` file is retained for legacy use.

For TDHF, use the same structure without `[input] functional`:
[`examples/TDHF/H2O_TDHF_ENERGY.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/TDHF/H2O_TDHF_ENERGY.oqp).

## Gradient

In ordinary TDHF/TDDFT, `S1` means the first excited state.

`.oqp`:

```text
tddft/b3lyp5/6-31g* grad(S3)
geom="h2o.xyz"
```

The state label sets the root count: `grad(S3)` makes OpenQP solve at least
three roots, so the route needs `nstate` only to request more roots than the
driver implies (for example `tddft(nstate=6)/b3lyp5/6-31g* grad(S3)`). The
same rule applies to `energy(...)`, `opt(...)`, and the SF/MRSF labels
(`grad(root=3)`, `grad(S2)`).

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_tddft_grad", silent=1)
job.molecule(geometry="water", charge=0, multiplicity=1)
job.theory.tddft(functional="b3lyp5", basis="6-31g*", nstate=3)
job.workflow.gradient(state=3)

mol = job.run()
gradient = mol.get_grad()
```

Legacy `.inp`:

```ini
[input]
runtype=grad
method=tdhf
functional=b3lyp5
basis=6-31g*

[scf]
type=rhf
multiplicity=1

[tdhf]
nstate=3

[properties]
grad=3
```

Runnable `.oqp`:
[`examples/TDDFT/H2O_B3LYP5-TDDFT_GRADIENT.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/TDDFT/H2O_B3LYP5-TDDFT_GRADIENT.oqp).
The same-stem `.inp` file is retained for legacy use.

TDHF gradient example:
[`examples/TDHF/H2O_TDHF_GRADIENT.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/TDHF/H2O_TDHF_GRADIENT.oqp).

## Notes

- In sectioned `.inp` input, `[tdhf] nstate` must include the highest state
  requested by gradients or follow-up analysis; `.oqp` input widens the root
  count from the driver's state label automatically.
- Ordinary TDHF/TDDFT state numbering is excited-state numbering: `state=1`
  in Python, or `grad=1` in `[properties]`, means the first excited state.
- SF-TDDFT and MRSF-TDDFT use a different target-state ordering: state `1` is
  the lowest spin-flip/MRSF target state, which can be the multiconfigurational
  ground state.
- Use `[tdhf] type=rpa` or omit `type` for the default response model. Use
  `[tdhf] type=tda` when a TDA calculation is desired.
- For TDHF in Python, use `job.theory.tdhf(basis=..., nstate=...)`.
