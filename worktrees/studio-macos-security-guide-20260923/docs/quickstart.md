# Quickstart

This quickstart runs a small MRSF-TDDFT calculation on water. The recommended
starting point is the readable `.oqp` input. The same calculation is then shown
with the Python API, followed by the legacy sectioned `.inp` format.

The `.oqp` parser is included in OpenQP 1.3.0 and is the recommended input
style for this quickstart. The [legacy `.inp` example](#legacy-inp-input)
remains available for compatibility.

## `.oqp` Input

Create `h2o.xyz`:

```text
3
water
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
```

Create `h2o_mrsf.oqp` beside it:

```text
mrsf(nstate=3)/bhhlyp/6-31g*
geom="h2o.xyz"
```

Read the file from top to bottom: use MRSF-TDDFT with BHHLYP/6-31G*, calculate
three states, and read the water geometry. No task keyword means a single-point
energy, and OpenQP selects the required high-spin working reference
automatically. Nothing else is needed: the initial guess, SCF convergence, and
grid settings are sensible defaults and are written only when they change.

Whitespace outside quotes and parentheses has no semantic effect, so the same
items may be split over several lines. Examples put the route and task on one
line and `geom` last for readability.

Run it:

```bash
openqp h2o_mrsf.oqp
```

For a gradient or geometry optimization, add `grad(S0)` or `opt(S0)` after the
route. For HF and DFT, the only ground-state surface is implicit, so
`grad` and `grad(S0)` are equivalent, as are `opt` and `opt(S0)`.

See [`.oqp` Input](oqp-input.md) for routes, physical state labels,
workflow controls, and more examples.

## Python Script

The same calculation can be set up from Python:

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

For MRSF-TDDFT, the Python theory helper supplies the required ROHF triplet
reference internally. HF and DFT scripts can set `multiplicity` directly in
`job.molecule(...)` when the molecular reference multiplicity is part of the
ordinary SCF setup.

OpenQP writes a log and structured output files in the working directory. For
more Python examples, see [Run OpenQP from Python](python-scripting.md).

## Legacy `.inp` Input

Existing sectioned inputs remain supported. The legacy spelling of the same
calculation is:

```ini
[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g*
functional=bhhlyp
method=tdhf

[guess]
type=huckel

[scf]
type=rohf
multiplicity=3

[tdhf]
type=mrsf
nstate=3
```

Save it as `h2o_mrsf.inp` and run `openqp h2o_mrsf.inp`. Use this format when
maintaining an existing input deck or when a legacy-only option is required;
new input examples in this manual lead with `.oqp`.

## Next Calculations

Use these `.oqp` inputs as nearby templates. Every linked `.oqp` example has
a same-stem legacy `.inp` companion.

| Goal | Recommended example |
| --- | --- |
| MRSF-TDDFT energy | [`examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.oqp) |
| RHF energy | [`examples/HF/H2O_RHF-HF_ENERGY.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/HF/H2O_RHF-HF_ENERGY.oqp) |
| MP2 energy | [`examples/MP2/h2o_ump2_6-31g.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/MP2/h2o_ump2_6-31g.oqp) |
| Analytic RHF-MP2 gradient | [`examples/MP2/h2o_rmp2_6-31g_grad.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/MP2/h2o_rmp2_6-31g_grad.oqp) |
| DFT gradient | [`examples/DFT/H2O_RHF-DFT_GRADIENT.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/DFT/H2O_RHF-DFT_GRADIENT.oqp) |
| Analytic HF/DFT Hessian | [`examples/HESS/H2O_RHF-DFT_ANA_HESS.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/HESS/H2O_RHF-DFT_ANA_HESS.oqp) |
| Native geometry optimization | [`examples/OPT/H2O_RHF-DFT_OPTIMIZE_OQP.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/OPT/H2O_RHF-DFT_OPTIMIZE_OQP.oqp) |
| SOC | [`examples/SOC/H2O_BHHLYP_SOC.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/SOC/H2O_BHHLYP_SOC.oqp) |
| PCM/ddX energy | [`examples/PCM/H2O_RHF-HF_DDPCM_ENERGY_ISPHER.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/PCM/H2O_RHF-HF_DDPCM_ENERGY_ISPHER.oqp) |
| NMR shielding | [`examples/NMR/H2O_RHF-NMR.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/NMR/H2O_RHF-NMR.oqp) |
