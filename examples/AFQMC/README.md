# OpenQP-AFQMC example

This directory contains a deliberately short end-to-end installation smoke
test. It uses a triplet ROHF/STO-3G CH2 reference, prepares the native OpenQP
Hamiltonian, and runs six AFQMC projection steps with 16 walkers.

Place AFQMC at `external/afqmc` (or in a sibling `afqmc` checkout), install
OpenQP with `pip install .`, and run:

```sh
openqp-afqmc CH2_ROHF_STO3G_AFQMC.oqp
```

Use `--prepare-only` to export `HAMILTONIAN`, `FCIDUMP`, `TRIAL`,
`OOORBDAT`, `AFQMC.json`, and `AFQMC.KEYWORDS` without launching the native
executable.
