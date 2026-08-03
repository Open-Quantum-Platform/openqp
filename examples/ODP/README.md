# ODP umbrella sampling and WHAM

`H2CO_BHHLYP-MRSF-NAMD-ODP.inp` is a minimal native ODP umbrella window. For a
production PMF, copy it into several windows, assign a unique `window` and
`center` in `[odp]`, and give each run its own packed `trajectory_file`.

Analyze all windows with the public Python API example:

```console
python odp_wham.py --temperature 300 --output odp-wham.npz window-*.namd.trj
```

The supplied temperature is the canonical temperature assumed by WHAM. The
current ODP production path is NVE, so users must decide whether their ensemble
preparation and sampling justify that analysis; OpenQP records this distinction
and emits a warning for non-NVT trajectory metadata.

To exercise checkpoint continuation end to end, run:

```console
python odp_restart.py --workdir odp-restart-smoke
```

The example performs one native ODP step, extends the generated canonical
`restart=true` manifest to step two, runs that manifest, and verifies the final
checkpoint. Omit `--workdir` for an automatically removed temporary run.
