An SCF continuation that exhausts its converger sequence is retried once from a fresh Huckel guess. Successful runs do not trigger recovery. SCF convergence remains mandatory; energy conservation and state-continuity checks remain active. The unit regressions deliberately exercise failed and recovered SCF calls.

With `md.disc_rescale=true`, energy recovery first retries the configured interval
with finer nuclear substeps (2, 4, 8, ... up to `md.disc_substeps`). If all configured
refinements fail the energy criterion, a converged electronic solution and positive
target kinetic energy permit isotropic numerical rescaling. The next interval uses
the configured dt again, with the same energy checks and recovery available.

The log records the uncorrected energy change, kinetic energies before/after the
correction, velocity factor, relative kinetic-energy change, and correction count.
Numerical corrections are not physical hops and are not compared with the hop-energy
tolerance. Stepwise and cumulative NVE checks remain active after correction.
`disc_tol` triggers refinement; it is not an upper bound on the correction.
No universal 0.5 eV acceptance rule is imposed.
