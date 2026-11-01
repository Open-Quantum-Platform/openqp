# Bounded XC response cache

Run the water MRSF gradient with a 256 MiB cache (the default):

```sh
OQP_XC_RESPONSE_CACHE_MB=256 OQP_XC_TIMING=1 openqp h2o.inp
```

Repeat in a separate output directory with `OQP_XC_RESPONSE_CACHE_MB=0` for
the uncached reference, or `1` to exercise partial storage. Numerical grid,
screening thresholds and precision remain unchanged. Missing blocks are
recomputed. The `[XCCACHE]` timing lines report retained bytes and block hits.

The cache stores fixed-reference LDA/GGA XC data for TDDFT Davidson, Z-vector
and CPHF iterations, and, when space remains, pruned AO values and derivatives.
The MRSF gradient also reuses AO spatial derivatives between its repeated grid
sweeps, including the moving-grid contribution. Run `h2o_tda.inp` or
`h2o_rpa.inp` to exercise the TDDFT Davidson and Z-vector paths. It is
released after the solve. The limit is per active solver/grid cache on each MPI rank, shared by its
OpenMP threads; it is not a limit on total process memory. Meta-GGA reuses AO blocks while recomputing the reference XC data.
The MRSF Davidson implementation has no repeated semilocal XC integration; its
response equations are unchanged.
