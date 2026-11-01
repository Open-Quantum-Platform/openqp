# ODP umbrella sampling and WHAM

`H2CO_BHHLYP-MRSF-NAMD-ODP.oqp` is a minimal native ODP umbrella window (NVE)
and `H2CO_BHHLYP-MRSF-NAMD-ODP-NVT.oqp` the same window sampled canonically.
For a production PMF, copy one into several windows, give each a unique
`window` and `center` in `odp(...)`, and give each run its own packed
`trajectory_file` in `md(...)`.

Analyze all windows with the public Python API example:

```console
python odp_wham.py --temperature 300 --output odp-wham.npz window-*.namd.trj
```

The supplied temperature is the canonical temperature assumed by WHAM. For a
free-energy window run each window canonically, `md(...,ensemble=nvt)`
(Langevin thermostat at `temperature`, default 300 K) as in `H2CO_BHHLYP-MRSF-NAMD-ODP-NVT.oqp`; WHAM then
requires the supplied temperature to equal the recorded thermostat
temperature. An NVE window is still accepted, but OpenQP records the ensemble
and warns that WHAM assumes canonical sampling.

ODP also runs in periodic QM/MM (`qmmm(cutoff=PME)`, orthorhombic cell). Each
CV is measured at its own minimum image: every atom pair a CV uses (an angle's
two arms from its vertex) is taken at its minimum-image separation,
independently of the other CVs, so CVs whose pairs form a ring around the cell
are measured correctly too. The run stops only when a pair is exactly half the
box length apart along an axis, where the minimum image is not unique.

To exercise checkpoint continuation end to end, run:

```console
python odp_restart.py --workdir odp-restart-smoke
```

The example performs one native ODP step, extends the generated canonical
`restart=true` manifest to step two, runs that manifest, and verifies the final
checkpoint. Omit `--workdir` for an automatically removed temporary run.
