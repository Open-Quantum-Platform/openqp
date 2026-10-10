# MRSF-EKT

MRSF-EKT computes ionization-potential and electron-affinity channels from an
MRSF-TDDFT reference.

`.oqp`:

```text
mrsf(nstate=10)/bhhlyp/6-31g ekt
geom="h2o.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_mrsf_ekt", silent=1)
job.molecule(geometry="water", charge=0)
job.theory.mrsf(functional="bhhlyp", basis="6-31g", nstate=10)
job.workflow.ekt(ip=True, ea=False)

mol = job.run()
```

Legacy `.inp`:

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

Runnable `.oqp` inputs:

- [`examples/other/h2o_rohf_mrsf_ekt_ip_6-31g_bhhlyp.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/other/h2o_rohf_mrsf_ekt_ip_6-31g_bhhlyp.oqp)
- [`examples/other/h2o_rohf_mrsf_ekt_ea_6-31g_bhhlyp.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/other/h2o_rohf_mrsf_ekt_ea_6-31g_bhhlyp.oqp)

Each has a same-stem legacy `.inp` companion.

At least one of `[ekt] ip` or `[ekt] ea` must be true. Set both to true to run
both channels in one job.

With `[guess] save_mol=True`, the full JSON contains state-specific Dyson AO
coefficients, labels, electron binding energies, and pole strengths. With
`[scf] save_molden=True`, OpenQP also writes a dedicated state-labeled Dyson
Molden file after the EKT calculation. See
[Orbital and Vibrational Output](orbital-vibrational-output.md).
