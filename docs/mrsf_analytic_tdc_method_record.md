# MRSF analytic time-derivative coupling method record

Status: scientific specification; no production implementation is implied by this record.

## Scientific objective

Evaluate the same-spin time-derivative coupling needed for nonadiabatic
dynamics,

\[
\tau_{IJ}(t)=\left\langle\Psi_I\middle|\frac{\partial\Psi_J}{\partial t}\right\rangle
             =\dot{\mathbf R}^{\mathsf T}\mathbf d_{IJ},
\qquad
\mathbf d_{IJ}=\left\langle\Psi_I\middle|\nabla_{\mathbf R}\Psi_J\right\rangle,
\]

without first storing all `3N` Cartesian components of `d_IJ`.  The velocity is
contracted with every nuclear-derivative contribution to the existing analytic
MRSF state-pair Lagrangian.  This is an exact directional derivative of the
implemented analytic NAC expression, not a truncated response theory.

## Electronic-structure specification

- Reference: the ordinary two-SOMO MRSF equal-weight mixed reference
  `rho_0^MR = [rho_0^(M_S=+1) + rho_0^(M_S=-1)]/2`.
- Orbital partition: closed `C`, open `O={O1,O2}`, and virtual `V`.
- Response space: the complete `CO`, `OV`, `CV`, and `OO` MRSF topology,
  including both parent contributions and the established parent-resolved
  open-open symmetrization.
- Target sector: two real, same-spin MRSF adiabatic states `I` and `J`.
- Perturbation: the time-even nuclear displacement contracted with the real
  nuclear velocity at one geometry.
- Observable: the real antisymmetric matrix `tau_IJ=-tau_JI`, in inverse atomic
  time.  The diagonal is zero.
- Electronic Hamiltonian and MRSF configuration basis: unchanged from the
  validated analytic full-vector NAC implementation.

The singlet open-open combination remains `(L-R)/sqrt(2)` and the triplet
combination remains `(L+R)/sqrt(2)` after parent phase alignment and averaging.
No CO, OV, CV, OO, spin-pairing, orbital-response, overlap-derivative, Pulay, or
nuclear term may be omitted merely because the final observable is a scalar.

## Inherited identities

For a nondegenerate pair,

\[
\mathbf h_{IJ}=(E_J-E_I)\mathbf d_{IJ},\qquad
\tau_{IJ}=\frac{\dot{\mathbf R}^{\mathsf T}\mathbf h_{IJ}}{E_J-E_I}.
\]

The directional implementation must reproduce, to rounding error, a posteriori
contraction of the existing full analytic vector with the identical velocity.
It must also approach the phase-aligned overlap time-derivative coupling as the
nuclear time step tends to zero.  A finite-step overlap coupling represents an
interval average; comparisons therefore use either a midpoint evaluation or a
time-symmetric average of endpoint analytic values.

## Distinct roles in surface hopping

The scalar `tau_IJ` controls electronic population transfer and the stochastic
hop probability.  It does not specify the nuclear momentum change after a hop.
When directional momentum rescaling is selected, the complete analytic vector
defines that direction.  With mass-weighted momentum `p` and mass-weighted unit
coupling direction `n_IJ`, the accepted-hop correction has the form

\[
\mathbf p' = \mathbf p + \alpha\mathbf n_{IJ},
\qquad
\tfrac12|\mathbf p'|^2 + E_J = \tfrac12|\mathbf p|^2 + E_I,
\]

where the energy-conserving root for `alpha` is used.  At a nonzero gap,
`h_IJ=(E_J-E_I)d_IJ` gives the same line and is numerically preferable near a
conical intersection.  The sign of this line does not affect the quadratic
energy equation, whereas its orientation relative to the nuclear momentum does.
Consequently, NAC-directed and isotropic/full-velocity rescaling may change the
number of frustrated hops, the post-hop momentum distribution, recrossing, and
eventual product branching.  This is tested as a physical dynamics choice and
is not folded into the validation of the electronic hop probability.

## New assumptions and optional approximations

The velocity-contracted analytic TDC itself adds no electronic-structure
approximation.  The following optional dynamics policies are separate and must
be identified in output and restart records:

1. `analytic-directional`: use the exact directional analytic TDC at every
   electronic step.
2. `adaptive`: use overlap TDC away from a transition and exact directional
   analytic TDC inside a hysteretic transition region.  The entry/exit criteria
   may depend on the energy gap and dimensionless coupling action
   `abs(tau_IJ)*dt`; hard state-index or geometry-specific rules are forbidden.
3. `analytic-on-hop`: propagate with overlap TDC and request the complete
   analytic vector only for an accepted hop or an explicitly requested
   NAC-directed momentum rescaling.
4. `local-lvc`: reuse an analytically determined local `g/h` branching plane.
   This is a model approximation and is out of scope until the exact directional
   implementation and its molecular comparisons pass.

## Response reuse between adjacent geometries

The first acceleration path preserves the exact analytic result.  Cache each
converged unordered-pair ROHF adjoint Z-vector and use it only as the initial
guess for the next MINRES solve.  Since the vector is represented in the MO
rotation space, it must first be transported to the new geometry using the AO
cross-overlap and phase-aligned maximum-overlap/Procrustes transformations
within the closed, open, and virtual subspaces.  State-pair permutation and
signs follow the same electronic-state overlap tracking used by dynamics.

Two predictors are evaluated:

\[
z_{n+1}^{(0)}=\mathcal T_{n\to n+1}z_n,
\]

and, after two accepted steps,

\[
z_{n+1}^{(0)}=\mathcal T_{n\to n+1}
\left[z_n+\eta\left(z_n-\mathcal T_{n-1\to n}z_{n-1}\right)\right],
\]

where `eta` is the ratio of consecutive nuclear time steps and is bounded to
avoid an unstable extrapolation.  MINRES then corrects this predictor until the
same certified true residual as a cold solve is reached.  This warm-start mode
therefore changes cost but not the final NAC within the solver tolerance.  It is
accepted only when it saves Hessian actions after including the extra `H z0`
residual evaluation.

The cache is invalidated upon a change in basis dimension, orbital occupation,
MRSF reference identity, tracked state manifold, failed subspace-overlap
criterion, large geometry displacement, nonfinite predictor, or increased
initial residual relative to the zero/preconditioned guess.  Cache contents,
transport overlaps, predictor type, initial/final residuals, and iteration
counts are restart data and are written to the calculation record.

A separately named approximate mode, `transport_approx` or `linear_approx`,
uses the transported predictor as the current Z vector and omits both MINRES
and the trial Hessian action.  The first point, a failed orbital/state-overlap
test, an excessive nuclear displacement, and every user-selected periodic
refresh point retain the exact solve.  All explicit current-geometry NAC terms
and the analytic HF/XC contraction of the approximate Z vector are still
evaluated.  A residual value of `-1` in the log means that the Z equation was
deliberately not tested at that step; it must not be read as convergence.

This solver-replacement mode is admissible for production dynamics only after
a paired corrected calculation has calibrated errors in `d`, `h`, and
`v dot d` against the finite-time-step error of the overlap/TLF propagation at
the same nuclear time step.  Periodic exact refreshes and immediate fallback
upon tracking failure limit cumulative drift.  Reusing derivative integrals,
XC grid terms, or an untransported Z vector at a displaced geometry is not
permitted as an unlabelled approximation.

TD-Baeck-An is retained as an independent, phase-free energy-curvature
diagnostic.  It is not an external reference for the signed analytic coupling
and must not replace `h` in the immediate vicinity of a conical intersection.

## Validity boundaries and diagnostics

- At exact degeneracy `d_IJ` is gauge-dependent and singular; the finite object
  is `h_IJ`.  Dynamics comparisons must record the gap, electronic time step,
  and state-manifold projector.
- An abrupt change of the triplet reference, a small T1-T2 reference gap, or a
  large state-specific response norm marks a possible breakdown of the local
  response description.  Such points remain in the record with diagnostics.
- State permutation and phase are determined by the established overlap
  tracking.  The same signed gauge is applied before comparing analytic and
  overlap TDC matrices.
- Translation removal is not applied to a physical velocity before forming the
  observable unless the dynamics integrator has explicitly removed center-of-
  mass motion.  The convention is recorded.

## Required evidence before molecular production

1. Preserve the existing independent MRSF model-space, dense-matrix, sigma,
   Davidson, and full analytic NAC evidence; the new path must not alter those
   results.
2. Unit-vector directional fixtures: contract each Cartesian unit direction
   and recover every component of the full analytic vector.
3. Random normalized directions: compare direct directional and stored-vector
   contractions for singlet and triplet pairs, with maximum absolute error at
   most `1e-10` atomic units on small fixtures.
4. Verify zero diagonal, antisymmetry, finite values, reference interchange,
   active-orbital rotation invariance, and rigid-translation behavior.
5. Compare phase-aligned analytic and centered overlap TDC at `dt`, `dt/2`, and
   `dt/4`, showing the expected convergence rather than accepting one time step.
6. Reproduce the H2/FCI directional coupling and at least one external MRCISD
   mode projection.
7. Only after items 1-6 pass, compare paired NAMD trajectories and adaptive
   policies.  Use identical initial conditions, random streams, gradients,
   decoherence, electronic substeps, and velocity-rescaling rule.
8. For response reuse, compare cold and transported warm starts over smooth and
   near-intersection paths.  Require identical certified residuals and NACs,
   record Hessian actions, and force cache invalidation at constructed orbital,
   reference, and state-switching events.  An approximate predictor-only mode
   requires a separate error calibration and must never serve as evidence for
   the exact method.

## Performance evidence

Report center electronic time, state-pair response time, directional derivative
assembly, full-vector derivative assembly, overlap-TDC time, and peak memory
separately.  A speed claim requires warm repetitions and must not compare a
single scalar to a full vector without naming the different observables.

## Provenance

- Starting OpenQP commit: `e619e38380f46499e277af5242b307e2a7f573d0` (PR-313 descendant, including merged PR
  311 and the signed NACME comparison criterion).
- Development branch: `agent/analytic-tdc-namd-20260829`.
- Central research task: `CHC-20260829-44893D`.
- Build, compiler, ILP64 linear-algebra, input, binary, and calculation hashes
  are mandatory in every numerical evidence package.
