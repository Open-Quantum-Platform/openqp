# oqp.quantum — Quantum-computing bridge

Export an OpenQP mean-field calculation as a **second-quantized molecular
Hamiltonian / FCIDUMP**, the universal entry point for quantum-computing
electronic-structure workflows (Qiskit Nature, OpenFermion, Block2/DMRG, …).

```python
from oqp.quantum import from_openqp

# `mol` is an OpenQP Molecule after a converged HF/DFT single point.
ham = from_openqp(mol)                     # ERIs computed natively
ham.to_fcidump("molecule.FCIDUMP")
```

## What it provides

| Piece | Status |
|-------|--------|
| One-electron MO integrals `h_pq` (from `OQP::Hcore` + MOs) | ✅ |
| Core / nuclear-repulsion energy, `nelec`, `MS2` metadata    | ✅ |
| AO→MO 1- and 2-index/4-index transforms                     | ✅ pure NumPy, tested |
| FCIDUMP write + read (8-fold symmetry, chemist notation)    | ✅ pure NumPy, tested |
| Two-electron AO integrals via native `oqp.int2e`            | ✅ (needs a build) |
| Two-electron MO integrals `(pq\|rs)` → full FCIDUMP         | ✅ |

**Verified:** the RHF energy reconstructed from an exported FCIDUMP matches
OpenQP's SCF energy to ~1e-14 on H2O/STO-3G and H2O/6-31G* (with d functions).
See `examples/QUANTUM/` and `tests/test_quantum_fcidump.py`.

## Two-electron integrals (`oqp.int2e`)

The ERIs are produced by a Fortran getter `oqp.int2e(mol)` that drives the
existing two-electron engine with a "dump" consumer and stores the full
`nbf**4` AO tensor (chemist notation) in the `OQP::ERI_AO` tag.
`from_openqp` calls it automatically. Because the integrals come from the same
engine and AO basis as `OQP::Hcore` and the MO coefficients, the resulting
Hamiltonian is self-consistent.

This is the **conventional in-core path** — memory grows as `nbf**4`, so it
targets small active systems / quantum-computing experiments, not production
basis sets (the routine aborts above ~16 GB). Override the source with
`eri_ao=`/`eri_provider=`, or skip the two-body part with `compute_eri=False`.

> Note: the Fortran (`source/modules/int2e.F90`, `oqp.int2e` C symbol, the
> `OQP::ERI_AO` tag) requires recompiling OpenQP. The pure-Python transforms
> and FCIDUMP I/O are tested independently and do not need a build.

## Conventions

* Two-electron integrals use **chemist notation** `(pq|rs)` (FCIDUMP / PySCF
  convention), tensor index order `[p, q, r, s]`.
* Symmetric one-electron OpenQP matrices are stored packed-triangular and are
  unpacked with `unpack_triangular`.

## Modules

* `integrals.py` — `unpack_triangular`, `ao_to_mo_1body`, `ao_to_mo_2body`
  (no dependency on the compiled `oqp` extension).
* `fcidump.py` — `write_fcidump`, `read_fcidump`.
* `hamiltonian.py` — `MolecularHamiltonian`, `from_openqp`.

See `examples/QUANTUM/export_fcidump.py` and `tests/test_quantum_fcidump.py`.
