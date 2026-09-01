# Harmonic FC example

After installing OpenQP, evaluate the explicit one-mode harmonic model with

```bash
python -m oqp.library.vibronic \
  examples/VIBRONIC/HARMONIC_FC.json \
  harmonic-fc-output.json
```

The example contains no electronic-structure calculation.  It exercises the
public multidimensional harmonic FC model, exact thermal population weights,
Gaussian broadening, and unit-sum spectrum normalization from explicit
frequencies, `J`, and `K`.

`MRSF_PROPERTY_FD.json` documents the validated request schema for target-state
MRSF normal-coordinate property derivatives.  It declares the one-based MRSF
root, the spin-adapted two-SOMO CO/OV/CV/OO representation, mode phases, units,
state-tracking thresholds, the fixed `h`, `h/2`, `h/4` convergence sequence,
and the explicitly truncated SOS Raman gate.  Use
`MRSFPropertyFDRequest.from_dict` to validate the record and
`run_openqp_mrsf_spectroscopy_fd` with a completed central OpenQP MRSF
calculation to execute it.
