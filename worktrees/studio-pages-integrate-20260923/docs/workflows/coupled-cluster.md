# Coupled Cluster (CCSD, CCSD(T))

!!! note "Available in OpenQP 1.3.0"
    Energy-only CCSD and CCSD(T) are included in OpenQP 1.3.0.

OpenQP supports ground-state coupled cluster with `[input] method=ccsd` and
`[input] method=ccsd(t)`. The calculation converges a Hartree-Fock reference
first, then adds the coupled-cluster singles and doubles correlation energy;
`ccsd(t)` additionally evaluates the perturbative triples correction.

Coupled cluster is currently an energy-only post-SCF workflow:

- Use `runtype=energy`.
- Leave `[input] functional` empty. Coupled cluster on Kohn-Sham references is
  rejected.
- Use `[scf] type=rhf`, `uhf`, or `rohf` for the reference.
- Use `[cc]` for the frozen core and the solver controls.

## Energy

### `.oqp`

```text
ccsd_t/6-31g cc(nfzc=1)
geom="h2o.xyz"
```

`ccsd_t`, `ccsd-t`, and `ccsdt` all select CCSD(T); `ccsd` stops after the
doubles. The route is model/basis and takes no functional slot, because
coupled cluster requires a Hartree-Fock reference.

The `reference` route option accepts `rhf`, `rohf`, or `uhf` and lowers to
`[scf] type`. It is separate from the `[cc]` solver section:

```text
ccsd_t(reference=uhf)/sto-3g mult=3
geom="ch2.xyz"
```

### Python

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_ccsd_t", silent=1)
job.molecule(geometry="water", charge=0, multiplicity=1)
job.theory.ccsd_t(basis="6-31g", reference="rhf", nfzc=1, conv=1.0e-7)

mol = job.run()
print("CCSD(T) total energy:", mol.get_results()["energy"])
```

`job.theory.ccsd(...)` is the same helper without the triples correction.
Note that `geometry="water"` is a built-in geometry and does not reproduce the
numbers tabulated below, which come from the explicit coordinates in the `.inp`
deck.

### Legacy `.inp`

```ini
[input]
system=
 8   0.000000000   0.000000000  -0.041061554
 1  -0.533194329   0.533194329  -0.614469223
 1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
method=ccsd(t)
functional=
basis=6-31g

[guess]
type=huckel

[scf]
type=rhf
multiplicity=1
conv=1.0e-10

[cc]
nfzc=1
conv=1.0e-7
```

Runnable inputs:
[`examples/CC/`](https://github.com/karmachoi/openqp/tree/claude/ccsd-t-parallel-lnoacw/examples/CC)
holds each calculation as both a `.oqp` line and a same-stem legacy `.inp`
deck. The link points at the preview branch of PR
[#302](https://github.com/Open-Quantum-Platform/openqp/pull/302), because that
directory does not exist on `main` yet; it becomes a `main` link once the
implementation merges.

Running the `.inp` deck shown above (`examples/CC/h2o_ccsd_t_6-31g.inp`)
reports:

| Quantity | Value (Ha) |
| --- | --- |
| `E(CCSD, correlation)` | `-0.1334141216` |
| `E((T), correction)` | `-0.0009510610` |
| `E(CCSD(T), correlation)` | `-0.1343651826` |
| `E(CCSD(T), total)` | `-76.1186552139` |

## References and Cost

`[scf] type` may be `rhf`, `uhf`, or `rohf`, but the two paths behind them are
very different in cost.

**Closed shell (RHF).** A spin-adapted formulation with one closed-shell
amplitude set. Every O(N^6) and O(N^7) contraction is cast as a large DGEMM,
the particle-particle ladder is blocked over the last virtual index, and the
triples are distributed over `a >= b >= c` virtual triples with MPI across
ranks and OpenMP within each. This is the path to use for production work.

**Open shell (UHF, ROHF).** A spin-orbital solver. It stores the full
`(2*nmo)^4` antisymmetrised integral tensor -- sixteen times the spatial one --
and is OpenMP-threaded but not MPI-distributed, so every rank would repeat the
same work. Use it for small systems.

The closed-shell path can also hold the ladder integrals as Cholesky vectors
instead of an explicit `v^4` array, and can build those vectors straight from
recomputed AO integrals rather than from a stored list. Both are decided on
memory rather than on speed, because both cost arithmetic to save it -- see
[`cholesky`](../keywords/cc.md#cholesky) and
[`cholesky_direct`](../keywords/cc.md#cholesky_direct). There is no disk-based
mode: whichever representation is chosen is held in memory.

Neither path uses a fixed memory ceiling. Each prints the storage it projects
and refuses when that exceeds what the machine can actually give, measured at
run time from physical RAM, the kernel's `MemAvailable`, and the cgroup limit --
so the same binary sizes itself correctly on a laptop and on a large node, and
sees a SLURM allocation rather than the whole node. In practice the closed-shell
path reaches a few hundred basis functions, and the open-shell one far fewer.

## Frozen Core

`[cc] nfzc` removes the lowest `nfzc` orbitals from the correlation treatment.

For an ROHF reference the orbitals are semicanonicalised before the amplitude
equations are solved, because the ROHF Fock matrix is not diagonal in its
occupied-occupied and virtual-virtual blocks. The frozen core is removed
*before* that rotation: the occupied rotation mixes core with valence, and
alpha and beta are rotated separately, so freezing afterwards would correlate a
different subspace in each spin. Freezing first keeps the correlated space
equal to the span of the reference orbitals `nfzc+1..nbf`, which is what
frozen-core coupled cluster means elsewhere.

## Implementation Notes

The AO integrals come from the shared two-electron engine through a collecting
consumer, are stored packed (only the canonical eighth of the tensor), and are
transformed to the MO basis by two half transformations driven through the
symmetric pair index, so neither the AO integrals nor the intermediate is ever
held as a dense `nbf^4` tensor.

For an ROHF reference the occupied-virtual Fock block survives
semicanonicalisation, and it is carried through both the amplitude equations
and the correlation energy -- Brillouin's theorem does not hold there, so the
singles term contributes.

Both paths are validated end to end against PySCF; the reference energies live
in `tests/data/ccsd_t_pyscf_validation.json` and
`tests/data/ccsd_t_open_shell_validation.json` in the main repository.

Derivative workflows are not implemented for coupled cluster. Use HF/DFT or
response methods for gradients, Hessians, optimization, NACME, SOC, PCM, NMR,
IR, and Raman workflows.
