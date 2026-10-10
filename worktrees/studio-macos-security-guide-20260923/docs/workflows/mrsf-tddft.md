# MRSF-TDDFT

!!! tip "Recommended `.oqp` input"

    OpenQP 1.3.0 includes the `.oqp` parser. Start with the physical calculation
    you want:

    ```text
    mrsf/bhhlyp/6-31g* opt
    geom="h2o.xyz"
    ```

    This optimizes `S0`. Write `opt(S1)` for the next singlet or `opt(T0)` for
    the lowest triplet. Do not enter a ROHF multiplicity: OpenQP selects the
    internal MRSF reference automatically. Add controls only when needed:

    ```text
    mrsf(nstate=3)/bhhlyp/6-31g* opt(S1,maxit=100) scf(conv=1e-08)
    geom="h2o.xyz"
    ```

    See the [`.oqp` Quick Start](../oqp-input.md#quick-start).

MRSF-TDDFT uses `[input] method=tdhf`, `[tdhf] type=mrsf`, and usually a triplet
ROHF reference. Although it is configured through the TDHF/TDDFT response
section, MRSF-TDDFT is used for multiconfigurational ground-state surfaces as
well as excited-state surfaces.

MRSF-TDDFT combines the two high-spin spin-flip reference components into a
mixed-reference response problem. This keeps the linear-response workflow close
to TDDFT while reducing the spin contamination that can obscure ordinary
spin-flip TDDFT roots. See the [References](../references.md#mrsf-tddft) page
for the original theory, analytic-gradient implementation, and recent accounts.

## Energy

The concise form defaults to the singlet manifold:

```text
mrsf(nstate=3)/bhhlyp/6-31g*
geom="h2o.xyz"
```

This calculates `S0`--`S2`. Select another physical manifold through the
driver, not through the route:

```text
mrsf(nstate=3)/bhhlyp/6-31g* energy(T0)
geom="h2o.xyz"
```

This calculates `T0`--`T2`. All-electron MRSF also accepts `energy(Q0)` for a
quintet manifold. Route parentheses accept `nstate` only.

For MRSF-TDHF, use the explicit basis-only route:

```text
mrsf-tdhf(nstate=3)/6-31g*
geom="h2o.xyz"
```

Related basis-only routes are `sf-tdhf(...)`, `umrsf-tdhf(...)`, and
`tda-tdhf(...)`; `cis(...)` is an accepted alias of `tda-tdhf(...)`.

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_mrsf", silent=1)
job.molecule(geometry="water", charge=0)
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=3)

mol = job.run()
results = mol.get_results()
print("Ground/reference energy:", results["energy"])
print("TD energies:", results["td_energies"])
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
type=mrsf
nstate=3
```

Runnable `.oqp`:
[`examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.oqp).
The same-stem `.inp` file is retained for legacy use.

For post-run NTOs, attachment/detachment densities, transition densities,
cube files, QCSchema export, FCIDUMP export, and external-code comparisons, see
[MRSF Analysis and Interoperability](mrsf-analysis.md).

## Gradient

In canonical input, select the physical state directly; an omitted state
defaults to `S0`:

```text
mrsf/bhhlyp/6-31g* grad(S2)
geom="h2o.xyz"
```

Here physical `S2` lowers to internal response root `3`, matching the
traditional and Python examples below.

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_mrsf_grad", silent=1)
job.molecule(geometry="water", charge=0)
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=3)
job.workflow.gradient(state=3)

mol = job.run()
gradient = mol.get_grad()
```

Legacy `.inp` uses `runtype=grad` and selects the internal response root
through `[properties] grad`:

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

Runnable `.oqp`:
[`examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_GRADIENT.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_GRADIENT.oqp).
The same-stem `.inp` file is retained for legacy use.

## Notes

- Logs distinguish the physical target (`S0`, `S1`, and so on) from the
  triplet ROHF reference and internal response-root index. The reference is an
  implementation requirement, not the spin or identity of the requested
  target state.
- MRSF states are zero-based within every physical spin manifold. `S0`, `T0`,
  and `Q0` lower to response root 1; `S1`, `T1`, and `Q1` lower to root 2.
- In Python, `job.theory.mrsf(...)` supplies the usual ROHF triplet
  reference for MRSF-TDDFT. Raw input files still show `[scf] multiplicity=3`
  explicitly because they are the direct OpenQP keyword form. In that form,
  `[tdhf] multiplicity=1`, `3`, or `5` is the physical response multiplicity,
  not the SCF reference multiplicity.
- MRSF state numbering follows the spin-flip/MRSF target-state list. `state=1`
  in Python, or `grad=1` in `[properties]`, means the lowest MRSF target state,
  which can be the multiconfigurational ground state. This differs from
  ordinary TDHF/TDDFT, where state `1` means the first excited state.
- In sectioned `.inp` input, `[tdhf] nstate` must include every state requested
  by gradients, NACME, SOC, or EKT analysis; `.oqp` input widens the root count
  from the driver's state labels automatically, so `mrsf/bhhlyp/6-31g* grad(S2)`
  solves three roots without a route `nstate`.
- For SOC, `mrsf(nstate=3)/... soc` requests `S0`--`S2` and `T0`--`T2`.
  Use `mrsf/... soc(ns=3,nt=5)` for unequal counts; `ns` and `nt` must appear
  together and cannot be combined with route `nstate`.
- NAMD uses the same zero-based labels. `namd(T0,soc=true)` means the first
  triplet MCH state and lowers to `[md] init_state=T0`; it must not be written
  as `T1` merely because the internal response root is 1.
- For ordinary TDDFT, see [TDDFT and TDHF](tddft.md).
- For spin-flip TDDFT without mixed-reference correction, use
  [SF-TDDFT](sf-tddft.md).
- UMRSF-TDDFT uses `[tdhf] type=umrsf` with a UHF reference and is currently an
  energy-only workflow.
