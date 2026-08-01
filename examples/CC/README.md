# Coupled cluster (CCSD, CCSD(T))

Energy-only ground-state coupled cluster on a Hartree–Fock reference.
`method=ccsd` stops after the doubles; `method=ccsd(t)` adds the perturbative
triples. Every example here runs in about a second.

Each calculation is given twice, once as a legacy `.inp` deck and once as a
semantic `.oqp` line, so the two formats can be compared directly.

| Calculation | `.inp` | `.oqp` |
| --- | --- | --- |
| H2O CCSD(T)/6-31G, frozen core | `h2o_ccsd_t_6-31g.inp` | `h2o_ccsd_t_6-31g.oqp` |
| H2O CCSD/STO-3G (no triples) | `h2o_ccsd_sto-3g.inp` | `h2o_ccsd_sto-3g.oqp` |
| CH2 triplet CCSD(T)/STO-3G, UHF | `ch2_triplet_ccsd_t_sto-3g.inp` | `ch2_triplet_ccsd_t_sto-3g.oqp` |
| OH doublet CCSD(T)/STO-3G, ROHF + frozen core | `oh_doublet_ccsd_t_rohf_fzc_sto-3g.inp` | `oh_doublet_ccsd_t_rohf_fzc_sto-3g.oqp` |

The `.oqp` form takes the reference as a model option and the solver controls
in a `cc(...)` call:

```
ccsd_t/6-31g h2o.xyz energy() cc(nfzc=1)
ccsd_t(reference=uhf)/sto-3g ch2.xyz mult=3 energy()
```

`ccsd_t`, `ccsd-t` and `ccsdt` all select CCSD(T); `ccsd` selects plain CCSD.
The route is model/basis — coupled cluster takes no functional, because it
requires a Hartree–Fock reference.

## `[cc]` keywords

| Keyword | Default | Meaning |
| --- | --- | --- |
| `nfzc` | 0 | Frozen core orbitals, excluded from the correlation treatment |
| `conv` | 1e-7 | Convergence on the amplitude RMS and the correlation energy |
| `maxit` | 50 | Maximum CCSD iterations |
| `ndiis` | 8 | DIIS subspace size; 0 disables DIIS |

## References and scaling

`[scf] type` may be `rhf`, `uhf` or `rohf`.

The closed-shell (RHF) path is the fast one: spin-adapted, every O(N^6)/O(N^7)
contraction cast as a DGEMM, OpenMP threaded and MPI distributed. It holds the
AO integrals in memory, so the practical ceiling is a few hundred basis
functions; the module prints the storage it needs and refuses above 64 GB.

Open-shell references (UHF, ROHF) go through a spin-orbital solver that stores
the full (2*nmo)^4 antisymmetrised tensor — sixteen times the spatial one — and
is OpenMP-only, with no MPI decomposition. It is meant for small systems; the
module prints its projected peak and refuses above 32 GB. An ROHF reference is
semicanonicalised first, and with `nfzc > 0` the core is removed *before* that
rotation, so the correlated space is the span of the reference orbitals that
were kept.

Reference energies for both paths, checked against PySCF, live in
`tests/data/ccsd_t_pyscf_validation.json` and
`tests/data/ccsd_t_open_shell_validation.json`.

## Python API

```python
from oqp.openqp import OpenQP

job = (OpenQP(project="h2o_ccsd_t")
       .molecule("h2o.xyz", basis="6-31g")
       .ccsd_t(reference="rhf", nfzc=1))
job.run()
```

`job.ccsd(...)` is the same call without the triples.
