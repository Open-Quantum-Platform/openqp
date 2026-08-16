# OpenQP text-log format

The OpenQP text log is both a human-readable calculation record and a
long-standing input to user analysis scripts. Its scientific labels, units,
state numbering, and stable markers therefore form a compatibility contract.

## Section order

Every top-level calculation uses the following order. Sections that do not
apply to a requested observable are omitted.

1. `RUN`
2. `INPUT AND REFERENCE`
3. `CONVERGENCE AND ITERATIONS`
4. `ENERGIES AND STATES`
5. `GRADIENTS AND PROPERTIES`
6. `TIMING AND TERMINATION`

The heading grammar is:

```text
   ========================================================================
   PyOQP LOG | INPUT AND REFERENCE
   PyOQP: Calculation request
   ========================================================================
```

The descriptive title remains method-specific. The category identifies the
kind of information that follows and is common to all methods.

## Common fields, units, and precision

Key/value records retain the stable `PyOQP` prefix and align the field label in
a 28-character column. Boolean values in new common fields use `yes` or `no`.

| Quantity | Log unit | Common precision |
| --- | --- | --- |
| Total electronic energy | Hartree | 10 digits after the decimal point |
| Nuclear gradient | Hartree/Bohr | 8 digits after the decimal point |
| Excitation energy | eV when present in a state table | 8 digits after the decimal point |
| Frequency | cm^-1 | 2 digits after the decimal point |
| IR intensity | km/mol | 6 digits after the decimal point |
| Wall time | days, hours, minutes, seconds | integer seconds |

Method-specific tolerances and residual norms retain scientific notation and
the precision of their solver. This is intentional: a residual is not an
energy and should not be formatted as one.

## Calculation-path coverage

The common section grammar surrounds the existing method-specific output for
all run types accepted by the main driver.

| Calculation family | Methods and run types | Method-specific records retained |
| --- | --- | --- |
| SCF references | HF and DFT with RHF, UHF, or ROHF; `energy`, `grad`, `hess`, `prop`, `data` | SCF iteration table, converger changes, orbital and population records |
| Post-SCF correlation | MP2 and spin-scaled variants; CCSD and CCSD(T) | Correlation components, iterations, reference and total energies |
| Linear response | TDHF/TDDFT, SF-TDDFT, MRSF-TDDFT, UMRSF-TDDFT | Davidson and Z-vector iterations, physical state labels, oscillator strengths |
| MRSF properties | `ekt`, `soc`, `nac`, `nacme`, `bp` | IP/EA roots, Dyson data, spin-orbit states, coupling matrices and vectors |
| Active-space methods | FCI, CASCI, CASSCF, SA-CASSCF | Active-space definition, CI roots and vectors, macroiteration convergence |
| Multireference perturbation theory | CASPT2, MS-CASPT2, XMS-CASPT2, NEVPT2, MRMP2, MCQDPT2, XMCQDPT2 | Variant-specific reference, state, and perturbative-energy information |
| Nuclear derivatives | Analytic and numerical gradients and Hessians; `grad`, `hess`, `thermo` | Displacement progress, gradients, normal modes, IR/Raman data, thermochemistry |
| Geometry and paths | `optimize`, `ts`, `meci`, `mecp`, `tci`, `mep`, `irc`, `neb` | Step number, energy, gradient criteria, state gaps, path termination reason |
| Tight-binding methods | DFTB, TD-DFTB, SF-TDDFTB, MRSF-TDDFTB, xTB | SCC, Davidson, Z-vector, state-spectrum, and backend diagnostics |
| Dynamics and embedding | `namd`, legacy `md`, QM/MM, SOC-NAMD | Time step, active state, hopping, energy conservation, embedding, termination |

MRSF state labels retain their established physical meaning. In particular,
the internal high-spin reference remains identified as an internal reference;
the common formatter does not reinterpret engine root numbers or alter spin,
energy, gradient, or coupling data.

## Compatibility policy

The following markers remain stable:

- `PyOQP started at`, `PyOQP build:`, and `PyOQP terminated at`;
- `PyOQP method:`, `PyOQP electronic energies`, and
  `PyOQP electronic gradients`;
- `PyOQP state`, CASSCF/FCI field names, and physical MRSF state labels;
- native `SCF`, Davidson, Z-vector, CPHF, CC, and displacement iteration
  markers;
- the legacy final-energy table column order, including additive DFTB columns.

The section heading and explicit unit records are additive. Energy values in
the common final table now use 10 rather than 8 digits after the decimal point,
matching the native SCF and active-space summaries. Numeric parsers should
continue to read the value as a floating-point number.

Do not parse separator width or spacing. Prefer the stable field markers above,
or JSON output when a structured representation is required. A future rename
or removal of a stable marker requires a documented migration period and a
regression test for both the old and new representations.

## Formatter ownership

`pyoqp/oqp/utils/log_format.py` owns the Python section grammar, key/value
alignment, common units, total-energy precision, and the Python equivalent of
the native module banner. `pyoqp/oqp/utils/state_labels.py` owns user-facing
method and state labels. Native solver tables remain close to the corresponding
Fortran solver because their columns describe method-specific numerical
quantities; the common Python sections place those tables in a consistent
calculation record without changing their values or iteration control.
