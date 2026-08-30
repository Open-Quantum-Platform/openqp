# Analytic-NAC FSSH method record

Status: implemented on the method-development branch. The native kernel,
gauge/scale, signed-TDC, compiled-runtime, and first molecular integration gates
have passed. Production molecular benchmarks remain subject to the staged gates
below.

## Scope

The method supplies two established quantities from the resident analytic
MRSF-TDDFT NAC implementation to fewest-switches surface hopping (FSSH):

1. the time-derivative coupling
   \(\sigma_{IJ}(t)=\dot{\mathbf R}(t)\mathbin{\cdot}\mathbf d_{IJ}(t)\), and
2. the direction used to conserve energy after an accepted hop.

It does not change the MRSF reference, response space, Hamiltonian, state
energies, gradients, analytic NAC Lagrangian, or electronic propagation
equations. The electron number, closed/open/virtual orbital partition,
equal-weight two-SOMO mixed reference, target singlet sector, response-vector
metric, and all fermionic phases are inherited unchanged from the analytic
MRSF-TDDFT energy, gradient, and NAC implementations. The observable added here
is a scalar inverse-time coupling and its use in a classical momentum update.
All electronic quantities and nuclear variables are in atomic units.

## Inherited equations and conventions

The analytic derivative coupling is

\[
 \mathbf d_{IJ}=\langle\Psi_I|\nabla_{\mathbf R}\Psi_J\rangle,
 \qquad \mathbf d_{JI}=-\mathbf d_{IJ},
\]

as produced by the certified resident MRSF NAC Lagrangian. State tracking first
places consecutive response vectors in one overlap-defined gauge. In that gauge,

\[
 \sigma_{IJ}=\sum_{A\alpha}v_{A\alpha}d_{IJ,A\alpha},
 \qquad \sigma_{JI}=-\sigma_{IJ}.
\]

For comparison with a finite-time overlap, the centered analytic reference for
the interval \([t_{n-1},t_n]\) is the trapezoidal value

\[
 \bar\sigma_{IJ}^{\rm an}=
 \tfrac12\left[\mathbf v_{n-1}\cdot\mathbf d_{IJ}(t_{n-1})+
                      \mathbf v_n\cdot\mathbf d_{IJ}(t_n)\right].
\]

Production propagation may use the endpoint value at \(t_n\); the output must
identify endpoint and centered quantities separately.

For a proposed hop \(I\rightarrow J\), let
\(\Delta E=E_J-E_I\). The directional velocity update is

\[
 \mathbf v'_A=\mathbf v_A+\gamma\frac{\mathbf d_{IJ,A}}{M_A}.
\]

Energy conservation gives

\[
 \tfrac12 A\gamma^2+B\gamma+\Delta E=0,
 \quad
 A=\sum_A\frac{|\mathbf d_{IJ,A}|^2}{M_A},
 \quad
 B=\sum_A\mathbf v_A\cdot\mathbf d_{IJ,A}.
\]

The accepted root is the real root with the smallest absolute momentum change.
An uphill hop is frustrated when \(B^2-2A\Delta E<0\). Multiplying
\(\mathbf d_{IJ}\) by any nonzero scalar, including a state-gauge sign, must
leave the final velocity invariant. A zero or non-finite direction is rejected,
not silently replaced by isotropic rescaling.

## Hop-triggered analytic NAC (HT-NAC)

HT-NAC is a computational approximation to the use of the analytic NAC at
every nuclear step; it is not an approximation to the resident analytic NAC
evaluated at a selected geometry.  Consecutive-state overlaps provide the NPI
time-derivative coupling used for electronic propagation and for the FSSH hop
probability.  The native FSSH kernel advances the amplitudes and decoherence
once, draws one random number, and returns an uncommitted candidate
\(I\rightarrow J\).  Only when such a candidate passes the energy-gap gate is
the exact resident analytic \(\mathbf d_{IJ}\) evaluated and used in the
directional velocity update above.  A successful update commits the state
change; a negative discriminant produces a frustrated hop.

The public input is

```ini
[md]
tdc = npi
rescale = hop_analytic_nac
```

Thus, HT-NAC supplies the physically defined momentum-rescaling direction at
an attempted transition while avoiding an analytic NAC calculation on steps
with no candidate.  It does not replace NPI hop probabilities by analytic
\(\dot{\mathbf R}\cdot\mathbf d\), and it does not provide a continuous NAC
record away from candidates.  Exact-NAC fields are explicitly empty on those
steps so that a preceding candidate vector cannot be mistaken for a current
quantity.

For the same overlap history, time step, amplitudes, random stream, and nuclear
state, HT-NAC must make the same stochastic candidate and final directional
rescaling decision as `tdc=npi, rescale=analytic_nac`.  It must use exactly one
electronic propagation and one random draw per step.  Analytic-NAC failure,
zero/non-finite direction, and a frustrated uphill update are hard recorded
failures of the candidate; there is no silent isotropic fallback.  The current
implementation inherits the same-spin, no-SOC, and no-QM/MM validity boundary
of analytic directional rescaling.

This addition changes neither the two-SOMO mixed reference nor the CO, OV, CV,
and OO MRSF response topology.  The electron number, orbital occupations,
target spin, response metric, operator, units, root tracking, and fermionic
phase convention are inherited unchanged.  The only new approximation is the
set of nuclear time points at which the already-defined analytic observable is
evaluated.

## Assumptions and validity boundary

- Only same-spin singlet MRSF-TDDFT is supported by the current analytic NAC.
- SCF and response convergence thresholds are those enforced by the analytic
  NAC driver; near a degeneracy, \(10^{-10}\) is recommended.
- The overlap-defined state permutation and phase must be applied before the
  analytic NAC is evaluated. A signed \(\sigma\) comparison is invalid across
  a tracking discontinuity.
- The derivative coupling is ill-conditioned as the adiabatic gap approaches
  zero. The finite branching-plane vector
  \(\mathbf h_{IJ}=(E_J-E_I)\mathbf d_{IJ}\) may define the same direction,
  but the scalar propagation coupling still requires a well-defined gauge and
  time discretization.
- Translation/rotation removal is an analysis convention. The first production
  implementation uses the full analytic vector and records any projected audit;
  a projected rescaling must be introduced as a separately named option.
- SOC, QM/MM constraint projection, complex electronic gauges, and degenerate
  multistate subspace propagation are outside the first implementation.

## Required gates before molecular use

1. Independent double-precision tests of the quadratic update: energy
   conservation, unchanged perpendicular velocity, gauge/scale invariance,
   downhill hops, frustrated uphill hops, and zero/non-finite directions.
2. Complete signed antisymmetry and zero-diagonal tests for
   \(\mathbf v\cdot\mathbf d\).
3. Central coordinate finite-difference validation of the inherited analytic
   NAC at several steps, with overlap/projector state tracking (existing gate).
4. Centered analytic-TDC comparison with the overlap TDC at decreasing nuclear
   time steps; endpoint and quadrature errors must not be conflated with an
   electronic-structure error.
5. Exact source/binary/build provenance and fail-on-skip test accounting.
6. HT-NAC candidate tests: no exact evaluation without a candidate; exactly one
   exact evaluation, one random draw, and one amplitude/decoherence propagation
   with a candidate; accepted/frustrated energy conservation; and no stale NAC
   record on a subsequent no-candidate step.
7. Only after gates 1--6 pass: H2/FCI, ordinary-point molecules, CI branching
   planes, then NAMD product-branching demonstrations.

## Initial verification record

- All seven independent two-SOMO MRSF reference fixtures pass with the
  established reference density, CO/OV/CV/OO response topology, metric, and
  fermionic phases unchanged.
- Native directional-rescaling tests pass energy conservation, gauge-sign and
  nonzero-scale invariance, downhill, frustrated-uphill, zero-vector, and
  non-finite-vector cases. Analytic TDC tests pass zero diagonal and signed
  antisymmetry.
- A one-step H2O MRSF-BH&HLYP/6-31G* trajectory using `tdc=analytic` and
  `rescale=analytic_nac` completed through state tracking, gradient, resident
  analytic NAC, the Fortran FSSH call, dense trajectory/restart output, and the
  NVE gate. The centered analytic-versus-overlap TDC RMS and maximum differences
  were (2.3008\times10^{-9}) and (3.8559\times10^{-9}) atomic units;
  the one-step total-energy drift was (1.6555\times10^{-8}) hartree.
- A corresponding H2 trajectory exposed an existing zero-dimensional
  open-shell gradient BLAS call (`LDC=0`, `M=0`) before NAC evaluation. The
  failure is retained as an H2-gradient edge-case diagnostic rather than being
  excluded or attributed to directional rescaling.
