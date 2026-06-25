# Quickstart

Create `h2o.inp`:

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

Run it:

```bash
openqp h2o.inp
```

OpenQP writes a log and structured output files in the working directory. For
examples with reference outputs, start from `examples/HF`, `examples/DFT`, and
`examples/MRSF-TDDFT`.

## Next Calculations

Use these examples as nearby templates:

| Goal | Example |
| --- | --- |
| RHF energy | `examples/HF/H2O_RHF-HF_ENERGY.inp` |
| DFT gradient | `examples/DFT/H2O_RHF-DFT_GRADIENT.inp` |
| Analytic HF/DFT Hessian | `examples/HESS/H2O_RHF-DFT_ANA_HESS.inp` |
| MRSF-TDDFT energy | `examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.inp` |
| Native geometry optimization | `examples/OPT/H2O_RHF-DFT_OPTIMIZE_OQP.inp` |
| SOC | `examples/SOC/H2O_BHHLYP_SOC.inp` |
| PCM/ddX energy | `examples/PCM/H2O_RHF-HF_DDPCM_ENERGY_ISPHER.inp` |
| NMR shielding | `examples/NMR/H2O_RHF-NMR.inp` |
