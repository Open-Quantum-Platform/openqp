# MRSF-TDDFTB (OpenQP-DFTB backend) examples

Plain long-range-corrected and uncorrected MRSF-TDDFTB through the optional
openqp-dftb backend (`[input] method=dftb`), without the DTCAM-TB operator
preset (for that, see `../DTCAM-TB`).

Parameters: nothing to configure with a current openqp-dftb wheel -- it
ships the bundled OB2W0PT3 SKF set (official shell-resolved `spinw.txt`
included, so the spin-polarization W kernels are complete), resolved
automatically when `[dftb] parameter_path` is empty and
`OPENQP_DFTB_PARAMETER_PATH` is unset. Point either knob at another
`.opdftb`/SKF directory to override the bundled default.

- `H2CH2_MRSF-TDDFTB_GRADIENT.inp`: analytic S0 gradient of CH2 on the
  triplet ROKS reference (relaxed Z-vector path).
- `C2_LC-MRSF-TDDFTB_ERF-TUNED_ENERGY.inp`: LC-corrected reference
  (`lc_ground_state=true`) with the erf-tuned response operator
  (omega=0.25, cam_beta=1.2) that restores the MRSF bright/dark ordering.

State convention: physical S0 is singlet response root 1; S1/S2 are roots
2/3. The high-spin ROKS determinant is only the auxiliary reference.
