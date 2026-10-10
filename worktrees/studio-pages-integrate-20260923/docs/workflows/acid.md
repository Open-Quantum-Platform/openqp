# ACID Current-Density Maps

ACID — the anisotropy of the induced current density, spelled "anisotropy of
the current-induced density" in the title of the original 2000 proceedings and
"Anisotropy of the Induced Current Density" in the 2005 *Chemical Reviews* paper
that is usually cited — maps where electrons are delocalised. It is the scalar field that an isosurface plot of aromaticity is
drawn from, and OpenQP writes it as Gaussian cubes together with the induced
current itself, so the map shows not only *where* the delocalisation is but
*which way* the current runs.

## What it needs

ACID is built from the magnetically perturbed density that the GIAO NMR
shielding calculation produces, so it is requested as a modifier on the
shielding call and requires `gauge=giao`.

A common gauge origin is not an option here. CGO shielding is gauge-contaminated
badly enough to invalidate the map: benzene NICS(0) comes out near −105 ppm in
6-31G* and is still ~43 ppm away from the GIAO answer at 6-311++G**. The
input checker refuses `acid` with `nmr_gauge=cgo` rather than producing a map
that cannot be compared with anything.

## Requesting it

`.oqp`:

```text
hf/6-31g* nmr(gauge=giao,acid=true)
geom="benzene.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("benzene_acid", silent=1)
job.molecule(geometry="benzene", charge=0, multiplicity=1)
job.theory.hf(basis="6-31g*")
job.workflow.nmr(gauge="giao", acid=True)

mol = job.run()
```

Legacy `.inp`:

```ini
[input]
runtype=energy
method=hf
basis=6-31g*

[scf]
type=rhf
multiplicity=1

[properties]
scf_prop=nmr,acid
nmr_gauge=giao
```

Runnable example:
[`examples/NMR/H2O_RHF-GIAO-ACID.inp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/NMR/H2O_RHF-GIAO-ACID.inp).

In the legacy `.inp` form `nmr` must come before `acid` in `scf_prop`; the
concise and Python surfaces order them for you.

## What it writes

Four Gaussian cubes, named after the log file:

| File | Contents |
| --- | --- |
| `<job>_acid.cube` | the ACID scalar field |
| `<job>_jx.cube`, `<job>_jy.cube`, `<job>_jz.cube` | the induced current-density vector, for **B** along *z* |

The ACID scalar is a tensor invariant and does not depend on the field
direction; only the three vector cubes do. They are ordinary Gaussian cubes, so
they load in VMD, Multiwfn and PyMOL, and the scalar cube renders as an
isosurface in [OQP Studio](../studio/analysis.md) like any other. Drawing the
current vectors on that isosurface is a separate Studio change and is not in the
shipped viewer yet.

## Grid controls

| Keyword | Default | Meaning |
| --- | --- | --- |
| [`acid_spacing`](../keywords/properties.md#acid_spacing) | `0.2` | grid spacing, bohr |
| [`acid_padding`](../keywords/properties.md#acid_padding) | `5.0` | padding around the molecule, bohr |

The defaults resolve a ring current comfortably. Coarsen them for a quick look
at a large molecule: the cube size grows as the cube of the inverse spacing, and
0.5 bohr with 3.0 bohr of padding turns a water map from about 2.9 MB per cube
into about 67 kB.

```text
hf/sto-3g nmr(gauge=giao,acid=true,acid_spacing=0.5,acid_padding=3.0)
geom="h2o.xyz"
```

## Reading the map

The standard plotting isosurface is **0.05 a.u.**, following Herges and
Geuenich ([references](../references.md#acid-anisotropy-of-the-induced-current-density)). On that scale a saturated bond falls below the threshold while
unsaturated and conjugated bonds stay above it, which is what makes the value
comparable across molecules. Critical isosurface values — the isovalue at which
the surface between two atoms breaks — at RHF/6-31G* GIAO:

| Bond | Critical isosurface value | vs 0.05 |
| --- | --- | --- |
| ethane C–C | 0.0475 | below: the surface has already broken |
| butadiene C1–C2 | 0.0558 | above |
| butadiene C2–C3 (the conjugated formal single bond) | 0.0620 | above |
| ethylene C=C | 0.0623 | above |
| benzene aromatic C–C | 0.0675 | above |

The current vectors give the circulation sense: diatropic (aromatic) current
runs clockwise seen from the **+z** side for **B** along *z*, paratropic
(antiaromatic) counter-clockwise. Benzene and cyclobutadiene separate cleanly on
this, with NICS(1)<sub>zz</sub> of −32.1 and +68.5 ppm respectively.

## Limits

- The grid evaluator walks Cartesian components up to *f*, so a
  spherical-harmonic basis needs `ispher=false` and a Cartesian basis containing
  *g* or higher shells is not supported. Both are refused with an explicit
  message rather than producing a partial map. Pople bases resolve to Cartesian
  shells under the default `ispher=auto` — the examples above run as written,
  and the log records `AO angular type: Cartesian (6d/10f/15g)` — while
  correlation-consistent bases such as cc-pVDZ resolve to pure spherical
  harmonics and need `ispher=false`.
- Ground state only. State-specific excited-state ACID is not exposed.
- Under MPI the cubes are written by the world root only.
- The map belongs to the calculation that produced it. The magnetic response is
  tied to the geometry, the basis and the orbitals it was computed from, and it
  is not reusable across a change in any of them: a new SCF drops it, and a
  moved geometry or a replaced basis is refused with a message naming what
  changed. Run the shielding again for a current map — that is a feature, not a
  cache miss, since combining one geometry's response with another's
  coordinates produces a map that is wrong without looking wrong.
