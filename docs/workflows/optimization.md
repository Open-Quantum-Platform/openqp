# Geometry Optimization

The default optimizer backend is the native OpenQP optimizer:

```ini
[optimize]
lib=oqp
```

The native backend uses redundant internal coordinates, DLC/TRIC-style
coordinate handling, and restricted-step RFO/P-RFO style steps. It supports
minima, transition states, MECI/MECP/TCI, NEB, IRC, and MEP workflows where the
corresponding run type is implemented.

## Native Minimum Search

```ini
[input]
runtype=optimize
method=hf
functional=bhhlyp
basis=6-31g*

[scf]
type=rhf
multiplicity=1

[optimize]
lib=oqp
istate=0
maxit=10

[oqp]
coordsys=tric
trust=0.2
```

Runnable reference: `examples/OPT/H2O_RHF-DFT_OPTIMIZE_OQP.inp`.

## geomeTRIC Backend

Use geomeTRIC when you need a workflow or constraint style that is best handled
by that external optimizer:

```ini
[optimize]
lib=geometric

[geometric]
coordsys=tric
trust=0.1
constraints_file=my.constraints
```

Runnable references are available in `examples/OPT/*_GEOMETRIC.inp`.

## Crossing Points

MECI and related workflows select the state pair or triplet in `[optimize]`:

```ini
[input]
runtype=meci
method=tdhf

[tdhf]
type=mrsf
nstate=5

[optimize]
lib=oqp
istate=1
jstate=2
```

Runnable references are available in `examples/OPT`.
