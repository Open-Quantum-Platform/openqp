# Single-run NAC validation harness

Validates the analytical MRSF NAC (`NAC.analytical_nac`) against the
numerical finite-difference NAC **within one process**, so both use the
same converged reference SCF/Davidson state (removing the run-to-run
phase/orientation confound of comparing two separate top-level runs).

## Run
```
export LD_PRELOAD=<gcc>/lib64/libstdc++.so.6
export LD_LIBRARY_PATH=<gcc>/lib64:~/.local/lib
export OPENQP_ROOT=~/openqp-nac OMP_NUM_THREADS=2
python single_run_validation.py
```
IMPORTANT: `import oqp` BEFORE `numpy` (otherwise numpy's LP64 LAPACK is
interposed and the ROHF Huckel guess DSYEVD fails with INFO=-2^32).

## Result (H2O MRSF-BHHLYP/6-31G*, states 2/3)
- numerical |d| ~ 1.057 /Bohr, respects C2v.
- CI part (gradient polarization): cos(num) ~ +0.81, correct C2v, but
  |d| ~ 0.009 (~1% of total). Direction right, magnitude small.
- frozen overlap part: breaks C2v, cos ~ 0.11 -- NOT a bug in isolation:
  the frozen S^x-half and the (missing) CPHF U^x parts of <phi|d_x phi>
  individually break point-group symmetry; only their SUM is symmetric.
  => the interstate-transition-density CPHF (U^x . gamma^IJ) term is
  required and large; it cannot be validated piecewise.
