# Troubleshooting

## OpenQP Cannot Find Runtime Files

For normal installs, run:

```bash
pip install .
```

from the repository root, or install the package from PyPI. The installed
package should locate its native library and data files without `OPENQP_ROOT`.
Only set `OPENQP_ROOT` for custom layouts where Python is separated from the
OpenQP runtime tree.

## BLAS/LAPACK Integer ABI

The default build expects ILP64 BLAS/LAPACK. LP64
(`-DLINALG_LIB_INT64=OFF`) is supported only on macOS, mainly for native
Accelerate builds. On Linux, use ILP64 OpenBLAS, MKL, or the configured bundled
Netlib path.

## SCF Does Not Converge

Try, in increasing order of intervention:

```ini
[scf]
diis_type=vdiis
```

```ini
[scf]
converger_type=soscf
```

```ini
[scf]
converger_type=trah
```

For difficult cases, use `alternative_scf=trah`, `escalation=soscf,trah`, MOM,
or pFON as appropriate.

## Too Many Threads

Use either:

```bash
openqp input.inp --omp 16
```

or:

```ini
[input]
omp_threads=16
```

This controls OpenMP threads per process or MPI rank. MPI rank count remains
controlled by the launcher, for example `mpirun -np 4`.

## PCM Input Runs as Vacuum

PCM-enabled runs should fail clearly if the ddX runtime path is unavailable.
Use:

```ini
[pcm]
enabled=true
backend=ddx
mode=reference_scf
model=ddpcm
```

and keep the run type at `energy` for the current production PCM path.

## Documentation Looks Stale

The keyword reference should be checked against:

- `pyoqp/oqp/molecule/oqpdata.py`
- `pyoqp/oqp/utils/input_checker.py`

If the code changes a default, allowed value, or runtime scope, update the
matching keyword and workflow page in the same pull request.
