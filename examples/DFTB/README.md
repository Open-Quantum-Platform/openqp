# DFTB examples

These inputs use the OpenQP-DFTB tight-binding backend (`method=dftb`). They are
**documentation examples**, not self-contained regression tests: they need the
optional external openqp-dftb library (`libopenqp_dftb_c`) and a DFTB parameter
set, neither of which OpenQP or its CI ships. `openqp --run_tests all` therefore
skips every `method=dftb` input; the operator presets themselves are covered by
the openqp-dftb test suite.

Run one explicitly once a compatible openqp-dftb is installed, meaning one
that implements the preset the input uses:

```bash
openqp examples/DFTB/CH2_MRSF-TDDFTB_DTCAM-GAP_ENERGY.inp
```

If the installed openqp-dftb does not bundle a parameter set, set
`[dftb] parameter_path` (an `.opdftb` file or an SKF directory) or export
`OPENQP_DFTB_PARAMETER_PATH`.

| Input | What it shows |
|---|---|
| `CH2_MRSF-TDDFTB_DTCAM-GAP_ENERGY.inp` / `.oqp` | MRSF-TDDFTB with `[dftb] model=dtcam-gap`, the DTCAM-TB operator refitted on singlet-triplet gaps (needs Open-Quantum-Platform/openqp-dftb#37) |
