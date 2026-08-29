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
