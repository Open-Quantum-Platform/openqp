# Smaller-time-step continuation

Run `oqp source.oqp` and then `oqp continuation.oqp` from this directory in a
validated OpenQP environment. This short formaldehyde example first saves step 1
at 0.1 fs, then propagates two 0.05 fs steps to 0.2 fs in separate outputs.
It demonstrates the interface; it is not a converged photochemical calculation.

`nstep` is the absolute final step index, not the number of added steps. A
checkpoint at step 7, time 0.7 fs continued with dt=0.05 and nstep=9 therefore
ends at 0.8 fs. Time origin is retained in the new trajectory and checkpoint.
Both continuation paths are required, and `restart=true` must not be supplied.
Paths are resolved from the working directory; absolute paths are recommended.
All new output paths must be unused and must not alias either source file.
The committed trajectory prefix is verified; any failed trailing source record
is left untouched. The generated child restart input uses ordinary restart.

Only fixed-step, same-spin analytic-TDC NAMD is supported. Hamiltonian,
electronic solver, acceptance, and random-seed settings must remain unchanged;
only a strictly smaller dt is permitted. The full saved electronic/reference
state, acceleration, amplitudes, phase history, energy history, and random
counter are restored. SCF acceptance criteria are unchanged. Different time
discretizations can encounter different hopping events despite the preserved
random counter. This is an explicit local continuation, not an adaptive-step
controller, and does not guarantee reconnection to an earlier trajectory.
