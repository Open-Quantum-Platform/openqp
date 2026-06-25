# Spin-Orbit Coupling

OpenQP supports MRSF-TDDFT spin-orbit coupling through `runtype=soc`.

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

[tdhf]
type=mrsf
nstate=12
multiplicity=3
```

Runnable references:

- `examples/SOC/H2O_BHHLYP_SOC.inp`
- `examples/SOC/CH3Br-BHHLYP-SOC.inp`

## SOC Terms

`[input] soc_2e` controls the SOC Hamiltonian:

| Value | Meaning |
| --- | --- |
| `0` | One-electron SOC terms only. |
| `1` | One-electron plus mean-field two-electron SOC terms. |

SOC workflows require a triplet ROHF reference and `[tdhf] type=mrsf`.
