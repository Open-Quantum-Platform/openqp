# NMR, IR, and Raman

## NMR Shielding

Request NMR shielding through `[properties] scf_prop=nmr`.

```ini
[input]
runtype=energy
method=hf
basis=sto-3g

[scf]
type=rhf
multiplicity=1

[properties]
scf_prop=nmr
nmr_gauge=cgo
```

Runnable reference: `examples/NMR/H2O_RHF-NMR.inp`.

`nmr_gauge` accepts:

| Value | Meaning |
| --- | --- |
| `cgo` | Common-gauge-origin shielding. |
| `giao` | Gauge-including atomic orbital shielding where supported. |

## IR and Raman

IR and Raman intensities are produced from supported Hessian/frequency
workflows. Start from the Hessian examples:

- `examples/HESS/H2O_RHF-DFT_ANA_HESS.inp`
- `examples/HESS/H2O_RHF-DFT_NUM_HESS.inp`
