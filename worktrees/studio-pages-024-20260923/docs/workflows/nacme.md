# NACME

NACME calculations use two geometries and an MRSF-TDDFT response calculation.

`.oqp`:

```text
mrsf(nstate=10)/bhhlyp/6-31g nacme(S0,S1)
geom="h2o.xyz"
geom2="h2o_previous.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

system = """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
"""

system2 = """
O   0.000000000   0.000000000  -0.031061554
H  -0.543194329   0.543194329  -0.624469223
H   0.543194329  -0.543194329  -0.624469223
"""

job = OpenQP("h2o_nacme", silent=1)
job.molecule(system, system2, charge=0)
job.theory.mrsf(functional="bhhlyp", basis="6-31g", nstate=10)
job.workflow.nacme()

mol = job.run()
```

Legacy `.inp`:

```ini
[input]
runtype=nacme
method=tdhf
functional=bhhlyp
basis=6-31g
system=
  O   0.000000000   0.000000000  -0.041061554
  H  -0.533194329   0.533194329  -0.614469223
  H   0.533194329  -0.533194329  -0.614469223
system2=
  O   0.000000000   0.000000000  -0.031061554
  H  -0.543194329   0.543194329  -0.624469223
  H   0.543194329  -0.543194329  -0.624469223

[scf]
type=rohf
multiplicity=3

[tdhf]
type=mrsf
nstate=10
```

Runnable `.oqp`:
[`examples/other/h2o_nacme_rohf_mrsf-s_6-31g_bhhlyp.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/other/h2o_nacme_rohf_mrsf-s_6-31g_bhhlyp.oqp).
The same-stem `.inp` file is retained for legacy use.

## Notes

- Keep `[tdhf] nstate` large enough to cover all states used in the coupling
  analysis.
- Use consistent atom ordering between `system` and `system2`.
- The `[nac]` section controls finite-difference and restart behavior for NAC
  workflows.
