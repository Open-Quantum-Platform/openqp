# Periodic full-ESPF QM/MM: Ewald electrostatics

`pyoqp/oqp/library/qmmm_ewald.py` provides the lattice-summed electrostatics
of the periodic ESPF QM/MM model (Bonfrate, Ferre, Huix-Rotllant, JCTC 2024,
20, 4338, eqs 7-8) for orthorhombic cells with tin-foil boundary conditions,
using standard Ewald sums (real-space erfc within half the shortest box edge,
reciprocal sum truncated at 1e-9). The full-ESPF QM/MM driver
(`OpenQpQMMM`, `embedding=electrostatic`) uses it whenever the OpenMM cutoff
is periodic (`cutoff=PME`, `Ewald`, ...); non-periodic systems keep the direct
Coulomb sum. It replaces the earlier minimum-image real-space sum, which was
conservative but had no long-range part.

## What is computed

* `Phi^MM_A`: MM potential (all MM charges and all their images) at every QM
  centre including link hydrogens, plus its gradient with respect to the QM
  centre. The frontier-charge redistribution (`frontier_scheme`) is applied to
  the MM charge set before the sum, as in the non-periodic case.
* Force on every MM atom from the QM ESPF charges and their images.
* `psi_img(A,B)`: the Ewald pair potential between QM centre A and the images
  of centre B (in-cell 1/r removed; diagonal = self-image term
  `sum_k w_k - 2 beta/sqrt(pi)`) with its gradient.
* A uniform neutralising background `-pi Q/(beta^2 V)` for non-neutral charge
  sets (a cut bond leaves the MM set non-neutral). It is a constant, so it
  carries no force. The paper's uniform charge correction `q_corr` (its eqs
  25-26) is **not** applied.

## Self-consistent QM-image term

The QM charges interact with their own images, `E_img = 1/2 q^T psi_img q`,
and `q` depends on the field. The driver solves this from outside the SCF:
the embedded SCF runs in `phi_eff = Phi^MM + psi_img q` and is repeated until
the ESPF charges are stable to 1e-7 e (3-4 iterations; the previous MD step's
charges seed the loop). If the loop has not converged after 50 iterations the
run stops with a `RuntimeError` rather than continuing with an inconsistent
energy/force pair (`OpenQpQMMM.IMAGE_MAXITER` / `IMAGE_TOL`). The total energy is

    E = E_QM[phi_eff] + Z.phi_eff - 1/2 q^T psi_img q + E_MM

and since dE_QM/dphi_A = -q_A the charge-response terms cancel at
self-consistency, so the force is the embedded QM gradient at fixed `phi_eff`
+ the Ewald coupling force (QM and MM sides) + the image force at fixed charges
`F_A = -q_A sum_B q_B grad_A psi_img(A,B)`. The earlier route through
`OQP::POTQM` / `add_potqm_contributions` inside the Fock build is not used
(its energy bookkeeping was never verified and it had no force).

The NAMD driver (`NAMD_QMMM`) runs the reference-density loop only to seed
the field, then iterates SCF -> MRSF -> active-state gradient until the
**relaxed ESPF charges of the propagated state** reproduce the field they
were computed in (`IMAGE_TOL_ACTIVE` = 1e-4 e, the Z-vector precision;
typically three gradient evaluations per step, the last of which is reused
as the force). The field the SCF sees therefore belongs to the state whose
force is integrated, and the response-term cancellation above holds for that
state. For the dipeptide box this moves the ESPF charges by 0.4 e and the
S0 energy by 6e-3 Ha relative to a field built from the ROHF reference
charges, and the NAMD force residual against finite differences improves
from 3.4 to 1.2 kJ/mol/nm at the same energy-conservation level (0.10 kJ/mol
over 50 fs). At a surface hop the force of the new state is evaluated in the
field of the previous state for that one step; the next step re-iterates.
The spin-adiabatic SOC-NAMD state is a mixture of MCH states whose relaxed
charges are not available per iteration, so periodic SOC-NAMD QM/MM is not
offered (`NotImplementedError`; use `NoCutoff`).

## Validation

`tests/test_qmmm_ewald.py` (pure numpy) checks the potential against an
explicit lattice sum and every analytic derivative against finite
differences. Driver level, alanine dipeptide with a link atom in a 1.6 nm
TIP3P box (337 atoms, HF/6-31G, AMBER-14, PME, OpenMM double precision,
PME tolerance 1e-6): the total force agrees with the finite difference of the
total energy to 0.023 kJ/mol/nm (5e-7 Ha/bohr) on QM, boundary and water
atoms; the QM-image self-consistency converges in 3-4 iterations. The ground-
state NVE in the box conserves energy at the same level as a pure-MM run of
the box with the same OpenMM Verlet integrator (flexible water, 0.5 fs). The
`split` embedding is not force-consistent under PBC (residuals of 1e3-1e4
kJ/mol/nm) and must not be used for dynamics.

Note: `OpenQpQMMM` in config mode used to write the QM geometry with six
decimals (Angstrom); the 5e-7 A rounding was a 0.3-0.9 kJ/mol/nm floor in
every finite-difference force test. It now writes twelve decimals.

## When the branch is used

The Ewald branch is entered for the periodic OpenMM nonbonded methods only
(`[qmmm] cutoff = PME`, `Ewald`, `LJPME`, `CutoffPeriodic`); `NoCutoff` and
`CutoffNonPeriodic` are treated as a finite cluster with the direct Coulomb
sum even when the PDB carries a `CRYST1` record. A periodic method with no box
vectors is an error, and only orthorhombic boxes are accepted. The
tight-binding NAMD path (`method=dftb/xtb`) has no ESPF charge operator to
iterate on and raises `NotImplementedError` under a periodic cutoff.

## Keys

* `[qmmm] ewald_tol` -- OpenMM PME/Ewald error tolerance for the MM-MM
  systems (default: OpenMM's 5e-4). Use 1e-6 for NVE validation.
* `[qmmm] lj_switch` -- switch the Lennard-Jones interactions smoothly to
  zero over the last 15% of the cutoff (default false = plain truncation).
* `[qmmm] h_lj` -- CHARMM TIP3P Lennard-Jones parameters on MM hydrogens that
  have none (default false). Too weak on its own to stop the collapse below.
* `[qmmm] mm_charge_width` (Angstrom, default off) -- Gaussian-smeared MM
  charges in the QM-MM electrostatics: the pair potential 1/r becomes
  erf(mu r)/r with mu = 1/(sqrt(2) w), consistently in energy, QM and MM
  forces, direct sum and Ewald real-space part (MM-MM and the QM-image term
  are untouched). This is the erf damping (ERFMU) of the reference Tinker
  ESPF code; w = 0.7 A corresponds to its default mu = 1/A. Under PBC the
  damping correction -erfc(mu r)/r is summed for the minimum image inside
  the real-space cutoff rc = L_min/2 only, so the width is limited to
  w <= rc/(4.5 sqrt 2) (1.26 A for a 16 A box); larger values are rejected
  (`EwaldQMMM.check_damping`). Point-charge
  embedding otherwise lets a TIP3P hydrogen collapse onto a QM carbonyl
  oxygen (2.25 -> 1.45 A in 90 fs of MRSF dynamics) through a polarisation
  runaway of the ESPF charges, after which the MRSF Z-vector equations
  diverge; with w = 0.7 A the Z-vector converges at that geometry and the
  spectrum is restored. Also used by the NAMD driver.

With both set, a periodic MRSF-TDDFT FSSH trajectory of the solvated alanine
link-atom system conserves energy to 0.05 kJ/mol fluctuation and -0.2
kJ/mol/ps drift over 88 fs (0.33 kJ/mol and +11 kJ/mol/ps with the OpenMM
defaults); the residual drift in solution is the classical truncation, not the
QM/MM coupling. The MRSF Z-vector divergence
seen after 60-90 fs in this solvated system with point charges is the
polarisation runaway described under `mm_charge_width`.

## Not implemented

* Lattice-parameter derivatives (pressure / NPT).
* Particle-mesh (PME) evaluation; the sums scale as N_QM*N_MM + N*n_k.
* The paper's non-neutral-cell charge correction.
* Non-orthorhombic cells.
