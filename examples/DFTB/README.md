# DFTB response-family examples

Ground-state DFTB, conventional TD-DFTB, SF-TDDFTB, and MRSF-TDDFTB all use
the optional openqp-dftb backend (`[input] method=dftb`). The four minimal
semantic inputs show the method routes without unrelated controls:

- `H2_DFTB_ENERGY.oqp`: SCC-DFTB ground-state energy.
- `H2O_TDDFTB_ENERGY.oqp`: conventional singlet TD-DFTB response. The current
  openqp-dftb response kernel is TDA, so logs identify this as
  `TD-DFTB (TDA)`.
- `CH2_SF-TDDFTB_ENERGY.oqp`: spin-flip TD-DFTB roots from the triplet ROKS
  reference. Use `root=N` when selecting one SF state because its spin
  character is assigned after the calculation.
- `CH2_MRSF-TDDFTB_ENERGY.oqp`: physical singlet S0 and higher MRSF states
  from the internal triplet ROKS reference.

Parameters: nothing to configure with a current openqp-dftb wheel -- it
ships the bundled OB2W0PT3 SKF set (official shell-resolved `spinw.txt`
included, so the spin-polarization W kernels are complete), resolved
automatically when `[dftb] parameter_path` is empty and
`OPENQP_DFTB_PARAMETER_PATH` is unset. Point either knob at another
`.opdftb`/SKF directory to override the bundled default.

Conventional TD-DFTB currently supports singlet targets only. Use the
SF-TDDFTB or MRSF-TDDFTB route for triplet-state work rather than requesting
a conventional TD-DFTB triplet.

The directory also retains two more detailed MRSF inputs:

- `H2CH2_MRSF-TDDFTB_GRADIENT.inp` / `.oqp`: analytic S0 gradient of CH2 on the
  triplet ROKS reference (relaxed Z-vector path).
- `C2_LC-MRSF-TDDFTB_ERF-TUNED_ENERGY.inp` / `.oqp`: LC-corrected reference
  (`lc_ground_state=true`) with the erf-tuned response operator
  (omega=0.25, cam_beta=1.2) that restores the MRSF bright/dark ordering.

MRSF state convention: physical S0 is singlet response root 1; S1/S2 are
roots 2/3. The high-spin ROKS determinant is only the auxiliary reference.

## Open-shell ground-state reference (`[dftb] reference=`)

A plain `dftb`/`dftb0` ground SCC can request an open-shell reference — computing
just the reference energy, no MRSF response — via the `reference` key:

- `reference=roks` — restricted-open-shell (Guest–Saunders)
- `reference=cuks` — constrained-UKS (default open-shell operator; ≡ ROKS on a
  common-orbital DFTB reference)
- `reference=uks`  — genuine unrestricted (independent α/β SCC; the DFTB+ `udftb`
  analogue). UKS lies at or below ROKS (variational).

`unpaired` (default 2) sets the number of unpaired electrons. The energy prints its
decomposition (H0 / SCC-gamma / spin / repulsive), for term-by-term comparison with
DFTB+. Examples: `CH2_UKS-DFTB_ENERGY.oqp`, `CH2_ROKS-DFTB_ENERGY.oqp`,
`CH2_CUKS-DFTB_ENERGY.oqp` (CH₂ triplet: UKS ≈ −2.5458, ROKS = CUKS ≈ −2.5419 Ha).
For an LC open-shell reference add `model=ob2` (LC has no closed-shell RHF form).
