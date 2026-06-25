# MRSF-EKT

MRSF-EKT computes ionization-potential and electron-affinity channels from an
MRSF-TDDFT reference.

```ini
[input]
runtype=ekt
method=tdhf
functional=bhhlyp
basis=6-31g

[scf]
type=rohf
multiplicity=3

[tdhf]
type=mrsf
nstate=10

[ekt]
ip=True
ea=False
```

Runnable references:

- `examples/other/h2o_rohf_mrsf_ekt_ip_6-31g_bhhlyp.inp`
- `examples/other/h2o_rohf_mrsf_ekt_ea_6-31g_bhhlyp.inp`

At least one of `[ekt] ip` or `[ekt] ea` must be true. Set both to true to run
both channels in one job.
