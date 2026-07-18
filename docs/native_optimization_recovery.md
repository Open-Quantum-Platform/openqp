# Native geometry and electronic recovery

OpenQP's concise `.oqp` format defaults to a vertical energy when no driver is
written:

```text
mrsf-tddftb(nstate=3)
dftb(model=dtcam-tb)
geom="molecule.xyz"
```

`opt`, `meci`, `mecp`, and `tci` use the native `lib=oqp` optimizer unless a
different backend is explicit. `lib=geometric` remains available:

```text
mrsf-tddftb(nstate=3)
opt(S1,lib=geometric,maxit=100)
dftb(model=dtcam-tb)
geom="molecule.xyz"
```

## Electronic recovery is local to one geometry

The electronic methods have different convergence machinery:

- MRSF-TDDFT starts each geometry with DIIS. If necessary it escalates through
  SOSCF to TRAH. The configured primary converger is restored before the next
  geometry.
- MRSF-TDDFTB starts each geometry with its configured SCC mixer (`auto` by
  default). A specifically diagnosed SCC exhaustion retries the current
  geometry with Broyden and then charge/orbital trust-region TRAH. The primary
  mixer is restored in a `finally` path, including failed calculations.

Parameter, ABI, gradient, and other non-SCC DFTB errors are not hidden by the
SCC ladder. If every electronic method fails at a trial geometry, the native
geometry optimizer rejects that geometry and continues from the best
successfully evaluated structure with the next geometry-recovery stage.

## Shared geometry recovery

The same bounded recovery structure is used for MRSF-TDDFT, MRSF-TDDFTB,
ordinary minima, penalty MECI/MECP, and two- or multistate BaekA:

1. Small and medium molecules start with DLC restricted-step RFO, a model
   Hessian, and Powell-damped BFGS updates.
2. Systems with at least 30 atoms start directly in DLC with `trust=0.02` and
   `trust_max=0.05`. C60 validation required 38 macroiterations with this
   profile, versus 54 for a TRIC-preconditioned handoff.
3. Oscillation, an iteration limit, or electronic failure triggers a restart
   from the best geometry with a fresh model Hessian and the alternate
   coordinate system.
4. A final bounded stage uses Cartesian coordinates and half the recovery
   trust (or returns to DLC if Cartesian was already used).

Explicit `coordsys`, `trust`, and `trust_max` values bypass the automatic
initial profile. `auto_recovery=false` disables the phase handoff and recovery
stages.

The rational/model-Hessian step follows the robust default strategy used by
DFTB+ geometry optimization. The electronic scanner behavior follows the same
principle as PySCF geometry wrappers: do not accept a gradient from an
unconverged reference, and apply a stronger electronic solver only to the
current geometry.

## MECI selection

`meci(S0,S1)` selects `algorithm=auto`: a short two-state penalty phase is
followed by BaekA if both the seam gap and geometry criteria are not satisfied.
BaekA supports two states as well as consecutive multistate intersections.
Three or more states select BaekA directly.
