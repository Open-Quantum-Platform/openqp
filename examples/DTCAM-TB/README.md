# DTCAM-TB (MRSF-TDDFTB) examples

`[dftb] model=dtcam-tb` selects the complete published DTCAM-TB operator and
its fast numerical protocol. The vector is resolved inside openqp-dftb
(single source of truth), so these inputs cannot drift from the paper:

- reference (ROKS) LC: (alpha_R, beta_R, omega_R) = (0, 0.04, 0.30 a0^-1), erf kernel
- independent response LC: (alpha_X, beta_X, omega_X) = (0, 1.0125, 0.2625 a0^-1)
- official OB2 v1.1.0 spin W with strengths (s_R^W, s_X^W) = (1.00, 0.6375)
- SPC (CO-CO, OV-OV, CO-OV) = (1.025, 0.25, 0.2625)
- response-only uniform one-center (ss, sp, pp) = (0, 0, -0.0125 Eh)
- numerics: Broyden 0.35 (history 12, max step 1.0), SCC tol 1e-8, Davidson
  response, relaxed Z-vector analytic gradients (single pair-space solve)

A model preset is all-inclusive: combining it with individual `[dftb]`
operator keys is rejected by the input checker. Omit `model=` to tune the
operator manually.

Parameters: point `[dftb] parameter_path` (or `OPENQP_DFTB_PARAMETER_PATH`)
at the OB2W0PT3 SKF directory containing the official `spinw.txt`.

State convention: physical S0 is singlet response root 1; S1/S2 are roots
2/3. The high-spin ROKS determinant is only the auxiliary reference.

MECI note: the example uses `meci_search=hybrid` (penalty -> updated
branching plane). A Baeck-An-type coupling-free algorithm is not yet
available in OpenQP's MECI drivers.

Files:
- `BUTADIENE_DTCAMTB_ENERGY.inp`      vertical S0/S1/S2 energies
- `BUTADIENE_DTCAMTB_S1_GRADIENT.inp` analytic S1 gradient
- `BUTADIENE_DTCAMTB_S1_OPTIMIZE.inp` S1 (1^1Bu+) minimum optimization
- `ETHYLENE_DTCAMTB_S1S0_MECI.inp`    S1/S0 twisted-pyramidalized MECI
