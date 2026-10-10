# SF-TDDFT

Spin-flip TDDFT uses a high-spin open-shell reference and `[tdhf] type=sf`.
It is useful for accessing low-spin target states from a high-spin reference.
MRSF-TDDFT is the recommended OpenQP route when the mixed-reference correction
is needed to reduce spin contamination and to describe multiconfigurational
ground-state surfaces.

In Python scripts, SF-TDDFT is a theory choice. Select the calculation type
separately with `job.workflow.*` only when you need a non-energy workflow.

!!! note "`.oqp` uses response roots"

    The spin character of an SF root is not known before diagonalization.
    State-specific `.oqp` input must therefore use `root=N`, for example:

    ```text
    sf(nstate=3)/bhhlyp/6-31g* grad(root=1)
    geom="h2o.xyz"
    ```

    Do not write an `S`/`T` label or omit the root on an SF state-specific
    driver.

## Energy

`.oqp`:

```text
sf(nstate=3)/bhhlyp/6-31g*
geom="h2o.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_sf", silent=1)
job.molecule(geometry="water", charge=0)
job.theory.sf_tddft(functional="bhhlyp", basis="6-31g*", nstate=3)

mol = job.run()
print("SF-TDDFT energies:", mol.get_td_energies())
```

Legacy `.inp`:

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
type=sf
nstate=3
```

Runnable `.oqp`:
[`examples/SF-TDDFT/H2O_BHHLYP-SFTDDFT_ENERGY.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/SF-TDDFT/H2O_BHHLYP-SFTDDFT_ENERGY.oqp).
The same-stem `.inp` file is retained for legacy use.

## Gradient

In SF-TDDFT, root `1` means the lowest spin-flip target state, not the first
ordinary TDDFT excited state.

`.oqp`:

```text
sf/bhhlyp/6-31g* grad(root=3)
geom="h2o.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_sf_grad", silent=1)
job.molecule(geometry="water", charge=0)
job.theory.sf_tddft(functional="bhhlyp", basis="6-31g*", nstate=3)
job.workflow.gradient(state=3)

mol = job.run()
gradient = mol.get_grad()
```

Legacy `.inp`:

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
type=sf
nstate=3

[properties]
grad=3
```

Runnable `.oqp`:
[`examples/SF-TDDFT/H2O_BHHLYP-SFTDDFT_GRADIENT.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/SF-TDDFT/H2O_BHHLYP-SFTDDFT_GRADIENT.oqp).
The same-stem `.inp` file is retained for legacy use.

## Notes

- SF-TDDFT currently uses an ROHF high-spin reference in the documented
  production path.
- In sectioned `.inp` input, `[tdhf] nstate` must include the highest spin-flip
  state requested by a gradient or follow-up workflow; `.oqp` input widens the
  root count from `root=N` automatically.
- For the mixed-reference correction and OpenQP's main multistate workflow, use
  [MRSF-TDDFT](mrsf-tddft.md).
