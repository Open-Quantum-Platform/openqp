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
