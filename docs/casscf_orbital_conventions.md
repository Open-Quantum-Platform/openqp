# What OpenQP's committed CASSCF orbitals mean

The CASSCF energy is invariant to rotations *within* the inactive, active and
virtual blocks, so a converged run still has to choose which representative of
that gauge to publish. OpenQP canonicalizes: it diagonalizes the inactive,
active and virtual blocks of the closed+active mean-field (generalized) Fock, so
the inactive and virtual orbitals carry orbital energies and the active orbitals
stay natural.

## The convention: the density is always root 0

The Fock that is diagonalized is built from the **root-0 CI density with weight
1.0**, and that is true regardless of what the run optimized:

* `[casscf] root=1` (state-specific excited CASSCF) — still root 0;
* SA-CASSCF with any weights — still root 0, not the weighted density.

Both implementations do this on purpose and in step:
`CASSCF._canonicalize` (`pyoqp/oqp/library/casscf.py`) calls
`_solve_active(..., np.array([1.0]), [0])`, and the native driver
(`source/modules/casscf_driver.F90`, `subroutine canonicalize`) mirrors it. They
must move together; if only one changed, the default native converger and the
Python fallback would disagree about the orbitals they publish, which is worse
than a consistent choice either way.

## What that does and does not affect

**Does not affect any energy.** No converged CASSCF or SA-CASSCF energy depends
on it, and no committed reference moves.

**Does not affect OpenQP's own PT2.** This is the important part, because it is
the usual reason the convention would matter. The perturbation theories do not
inherit these orbitals — each one semicanonicalizes from the density its own
zeroth-order Hamiltonian is defined against:

* MS-CASPT2 with the Fock `H0` uses the multi-set construction: every root is
  solved in semicanonical orbitals built from **that root's own** density
  (`pyoqp/oqp/library/caspt2_dyall.py`).
* Dyall-`H0` CASPT2, XMS-CASPT2 and the QDPT family (MRMP2 / MCQDPT2 /
  XMCQDPT2) use the single-set state-averaged canonical Fock, which is the
  GAMESS `IFORB=1` convention.
* The analytic PT2 gradient does not use the semicanonical set at all: the
  pre-semicanonicalization reference orbitals are preserved as
  `OQP::CASPT2_REFERENCE_MO`, because its orbital constraint (RHF canonicality
  or CASSCF stationarity) holds for those and not for the rotated set.

**Does affect what you read out of a run.** For an excited-state or
state-averaged solution the reported inactive and virtual orbital energies are
ground-state-like, and so is anything you build yourself on the committed
`OQP::VEC_MO_A` — an orbital plot, a Molden export, or your own post-processing.
If you need orbitals canonical with respect to the state the run actually
targeted, construct them from the state's own density rather than assuming the
published set already is.

See issue #338 for the discussion.
