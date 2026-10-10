# Geometry Optimization and Reaction Paths

Concise `.oqp` geometry drivers use the native OpenQP optimizer automatically.
Users select a physical state and the calculation they want; there is no
backend or `lib` keyword in this format.

```text
dft/pbe0/def2-svp opt
geom="h2o.xyz"
```

```text
mrsf(nstate=5)/bhhlyp/6-31g* meci(S0,S1)
geom="guess.xyz"
```

The native engine supports minima, transition states, two-state and BaekA
multistate MECI, MECP, MEP, IRC, and NEB. Minimum, crossing-point, and TS
optimizers use redundant internal, DLC, TRIC, or Cartesian coordinates as
appropriate, with restricted-step RFO/P-RFO optimization. When `coordsys` is
omitted or set to `auto`, these optimizers start in DLC. The DLC transformation
is accepted only when it spans the complete molecular vibrational space; rank
loss activates the Cartesian recovery route rather than silently removing a
vibrational direction. For molecular complexes, OpenQP adds interfragment
distances when needed to condition the internal-coordinate transformation.
MEP and IRC instead trace a mass-weighted path, and NEB optimizes its Cartesian
band with FIRE; those path drivers do not read `coordsys`.

## Native Minimum Search

The shortest canonical input is:

```text
dft/bhhlyp/6-31g* opt
geom="h2o.xyz"
```

`opt` defaults to `S0`. Add native controls directly when needed:

```text
dft/bhhlyp/6-31g* opt(coordsys=tric)
geom="h2o.xyz"
```

Python uses the native backend by default. Its explicit backend selector is
retained for compatibility with existing scripts:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_opt", silent=1)
job.molecule(geometry="water", charge=0, multiplicity=1)
job.theory.dft(functional="bhhlyp", basis="6-31g*")
job.workflow.optimize(istate=0, trust=0.2)
mol = job.run()
```

The equivalent legacy `.inp` spelling remains supported:

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

[oqp]
trust=0.2
trust_max=0.5
```

These examples use the automatic DLC default. Set `coordsys=tric`, `ric`, or
`cart` only when that coordinate representation is intentionally required.

Runnable `.oqp`:
[`examples/OPT/H2O_RHF-DFT_OPTIMIZE_OQP.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/OPT/H2O_RHF-DFT_OPTIMIZE_OQP.oqp).
The same-stem `.inp` file is retained for legacy use.

## Native Transition State and IRC

Native TS optimization uses P-RFO. `follow` selects the initial mode index, and
`hessian` selects how the starting Hessian is obtained:

```text
dft/pbe0/def2-svp ts(S0,hessian=numerical)
geom="ts_guess.xyz"
```

`hessian=model` is the inexpensive default. `numerical` or `analytical`
calculates a real Cartesian Hessian for the selected state before the first TS
step. Availability of an analytical Hessian still depends on the electronic
method and basis. For an isolated molecule in full-rank TRIC or Cartesian
coordinates, OpenQP removes whole-molecule translation/rotation zero-mode noise
from a real Hessian and restores positive model curvature in those rigid
directions before P-RFO mode selection; internal-only RIC/DLC modes are left
unchanged. The standalone engine leaves this projection off by default so an
external field can retain genuine lab-frame curvature; its `auto` coordinate
selection therefore uses full-rank TRIC rather than removing those physical
translations and rotations through DLC. Explicit `coordsys=dlc` remains
available when the standalone objective is known to be invariant. Active QM/MM
OpenQP geometry and reaction-path jobs are currently rejected in preflight
because their force backend is not connected to these optimizers; supported
QM/MM workflows remain energy, MD, and NAMD.

After locating a transition state, trace each native IRC branch; the forward
branch is the default direction, so only the backward run names it:

```text
dft/pbe0/def2-svp irc(S0,hessian=analytical)
geom="ts.xyz"
```

```text
dft/pbe0/def2-svp irc(S0,direction=backward,hessian=analytical)
geom="ts.xyz"
```

Native IRC projects mass-weighted translation and rotation modes, then requires
exactly one significant negative vibrational mode. It rejects a minimum (none)
or a higher-order saddle (more than one) before tracing the path.

Native MEP uses the same gradient stopping threshold without requiring a
transition-state Hessian:

```text
mrsf(nstate=5)/bhhlyp/6-31g* mep(S0)
geom="start.xyz"
```

## Native NEB

The reactant comes from `geom`; `product` supplies the second endpoint. Native
NEB can align the endpoints, relax them, run climbing-image NEB, test both
maximum and RMS force thresholds, and write the final band:

```text
dft/pbe0/def2-svp neb(S0,product="product.xyz",images=7,frms=0.001,output="reaction_path.xyz")
geom="reactant.xyz"
```

`climb`, `align`, and `opt_ends` are booleans. If `output` is omitted, OpenQP
writes `<project>_neb.xyz` in the log directory. The multi-frame XYZ file
contains every final image and records each image energy in Hartree.
When `climb=true`, set `climb_fmax >= fmax` so the climbing image activates
before the final convergence threshold can be satisfied.

## Crossing Points

Physical state labels replace internal state indices in `.oqp`:

```text
mrsf(nstate=5)/bhhlyp/6-31g* meci(S0,S1)
geom="guess.xyz"
```

```text
mrsf(nstate=5)/bhhlyp/6-31g* mecp(S0,T0)
geom="guess.xyz"
```

```text
mrsf(nstate=5)/bhhlyp/6-31g* meci(S0,S1,algorithm=baeka)
geom="guess.xyz"
```

```text
mrsf(nstate=5)/bhhlyp/6-31g* meci(S0,S1,S2)
geom="guess.xyz"
```

```text
mrsf(nstate=6)/bhhlyp/6-31g* meci(S0,S1,S2,S3)
geom="guess.xyz"
```

Two-state MECI and MECP both default to `algorithm=auglag`, an augmented
Lagrangian. A least-squares multiplier removes the mean-gradient component
along the gap direction, so what remains is a gap term and a projected mean
gradient that are orthogonal to each other. Orthogonal terms cannot cancel, so
a vanishing total gradient forces the energy gap to zero.

This matters because a plain penalty cannot do that. Its stationary condition
balances the mean gradient against the gap term, and since those two are not
orthogonal they cancel at a residual gap of order `1/weight`, leaving the
optimizer converged onto a point that is not a crossing. Raising the weight
shrinks the residual only as `1/weight` and never removes it. The legacy MECP
`algorithm=quad` behaves this way and is kept only for reproducing older runs.

Use [`gap_sigma`](../keywords/optimize.md#gap_sigma) to scale the gap term;
`gap_sigma=1` recovers the plain Bearpark gradient projection.

```text
mrsf(nstate=5)/bhhlyp/6-31g* mecp(S0,T0,algorithm=auglag)
geom="guess.xyz"
```

`algorithm=baeka` selects the Baek adaptive penalty-function method within the
ordinary MECI driver. It accepts two or more states from one spin manifold,
uses the independent adjacent energy gaps, and updates the penalty strength
additively. Existing `tci(S0,S1,S2,...)` routes remain supported as an
independent legacy three-state workflow; they are not aliases for BaekA and
retain their established multiplicative update. See
[BaekA Multistate MECI](baeka-multistate-meci.md) for the
method, production controls, and regression-example scope.

Traditional `.inp` and Python scripts may continue to use the internal
`istate`, `jstate`, `kstate`, `imult`, and `jmult` fields documented under
[`[optimize]`](../keywords/optimize.md).

## Native Frozen-Distance Constraints

Native minimum searches can freeze one or more initial atom-pair distances.
Atom indices are one-based, and multiple constraints are separated by
semicolons:

```text
dft/bhhlyp/3-21g opt(freeze="distance(1,2)",coordsys=dlc,trust=0.05,trust_max=0.05)
geom="hcn.xyz"
```

The distance is fixed at its value in the input geometry. The native optimizer
projects the gradient into the constraint tangent space and corrects every
trial geometry back to the requested distance. The traditional equivalent is:

```ini
[optimize]
lib=oqp

[oqp]
coordsys=dlc
trust=0.05
trust_max=0.05
freeze=distance(1,2)
```

Runnable regression:
[`HCN_RHF-DFT_CONSTRAINED_OQP.inp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/OPT/HCN_RHF-DFT_CONSTRAINED_OQP.inp).

## Optional Legacy geomeTRIC Compatibility

geomeTRIC is not part of the concise `.oqp` geometry grammar. It remains
available to traditional sectioned `.inp` files and the Python API for advanced
constraint types beyond the native frozen-distance control.
Install the optional dependency first:

```bash
pip install "openqp[geometric]"
```

The compatible Python spelling is:

```python
job.workflow.optimize(
    lib="geometric",
    istate=0,
    coordsys="tric",
    trust=0.1,
    constraints_file="my.constraints",
)
```

The legacy `.inp` spelling is:

```ini
[input]
runtype=optimize
method=hf
functional=pbe0
basis=def2-svp
system=reactant.xyz

[scf]
type=rhf
multiplicity=1

[optimize]
lib=geometric
istate=0

[geometric]
coordsys=tric
trust=0.1
constraints_file=my.constraints
```

No standard shipped regression now depends on geomeTRIC. General minimum,
frozen-distance, crossing-point, TS, IRC, MEP, and NEB examples use the native
engine.
