# Orbital and Vibrational Output

OpenQP can write molecular orbitals, vibrational normal modes, and
state-specific MRSF-EKT Dyson orbitals in portable JSON and Molden formats.
These files can be opened together in
[OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) so the same
calculation can display both MO surfaces and animated normal modes. Files are
processed locally by the browser.

## Portable JSON

Set [`[guess] save_mol=True`](../keywords/guess.md#save_mol) to retain the full
same-stem `.json` calculation file. When the AO basis uses an ordering supported
by Molden, the JSON includes:

- `basis_set`: contracted shells, Molden-normalized primitive coefficients,
  spherical/Cartesian convention, and AO ordering;
- `molecular_orbitals`: alpha and, where applicable, beta energies,
  occupancies, and AO coefficients;
- `dyson_orbitals`: state-specific IP/EA AO coefficients, labels, parent state,
  electron binding energy, and pole strength for MRSF-EKT calculations; and
- `frequency_modes`: frequencies in cm⁻¹ and mass-unweighted Cartesian normal
  mode displacements when a frequency calculation has completed.

Hessian workflows also write `<project>.hess.json`. It contains the Hessian,
frequencies, normal modes, IR intensities, Raman activities, atoms, geometry,
and—when available—the same portable basis and MO blocks. Consequently, either
the full `.json` or `.hess.json` can provide both the molecular orbitals and the
vibrational view.

High angular-momentum H/I shells and mixed Cartesian/spherical bases remain
valid OpenQP calculations, but their optional portable viewer blocks are
omitted when no unambiguous Molden AO ordering is available.

## Combined Molden MO and frequency file

A completed Hessian workflow writes `<project>.freq.molden`. When the basis and
SCF orbitals support Molden ordering, this is a combined file with:

- `[Atoms]`, `[GTO]`, spherical-harmonic markers, and `[MO]`;
- `[FREQ]` and standard one-value-per-mode `[INT]` IR intensities;
- optional `[RAMAN]` activities; and
- `[FR-COORD]` and `[FR-NORM-COORD]` normal-mode geometry and displacements.

The standard Molden sections keep the file usable by other Molden readers,
including applications that ignore the optional `[RAMAN]` extension. If the
basis cannot be represented in portable Molden ordering, OpenQP safely falls
back to a frequency-only Molden file.

## State-specific MRSF-EKT Dyson Molden

For MRSF-EKT, enable [`[scf] save_molden=True`](../keywords/scf.md#save_molden).
After the EKT kernel finishes, OpenQP writes a dedicated `dyson` Molden file.
Each orbital is labeled by channel and root, for example
`Dyson-IP-state-1` or `Dyson-EA-state-2`; its energy records the EKT eigenvalue
and its occupancy records the pole strength. OpenqpView exposes these labels as
individually selectable Dyson orbitals.

The H₂O IP example enables both JSON and Molden exports and includes generated
SCF and Dyson Molden results:

[`examples/other/h2o_rohf_mrsf_ekt_ip_6-31g_bhhlyp.oqp`](https://github.com/Open-Quantum-Platform/openqp/blob/main/examples/other/h2o_rohf_mrsf_ekt_ip_6-31g_bhhlyp.oqp)

## Minimal input settings

`.oqp` Hessian output with MO data:

```text
dft/bhhlyp/6-31g* hess(S0,type=analytical,clean=true) guess(save_mol=true)
geom="h2o.xyz"
```

Legacy `.inp` MRSF-EKT output:

```ini
[input]
runtype=ekt
method=tdhf
functional=bhhlyp
basis=6-31g

[guess]
save_mol=True

[scf]
type=rohf
multiplicity=3
save_molden=True

[tdhf]
type=mrsf
nstate=10

[ekt]
ip=True
ea=False
```

