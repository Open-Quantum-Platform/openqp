# OpenQP text-log format

The OpenQP text log is both a human-readable calculation record and a
long-standing input to user analysis scripts. Its scientific labels, units,
state numbering, and stable markers therefore form a compatibility contract.

## Section order

Every top-level calculation uses the following categories. Scientific stage
categories retain the listed order, and sections that do not apply to a
requested observable are omitted.

1. `RUN`
2. `INPUT AND REFERENCE`
3. `CALCULATION PROGRESS`
4. `CONVERGENCE AND ITERATIONS`
5. `ENERGIES AND STATES`
6. `GRADIENTS AND PROPERTIES`
7. `TIMING AND TERMINATION`

`CALCULATION PROGRESS` contains phase transitions and status messages such as
entering a gradient, geometry step, path point, or dynamics step. It may recur
between method-stage records; it does not contain final scientific results or
replace the convergence, energy/state, or gradient/property categories.

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
a 28-character column. Boolean values in common fields use `yes` or `no`.

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

## Verbosity

One level sets how much detail the whole log carries, in the Python driver and
in the native solvers alike. Set it with `[input] verbose`. The older spelling
`[scf] verbose` is still read; the parser gives both the default `1`, so the one
an input changes wins, `[input]` first. `[dftb] print_level` left at its default
follows the same level.

| Level | Name | What it adds |
| --- | --- | --- |
| `0` | quiet | Section headings, the calculation request, SCF and solver convergence results, final energies, gradients and properties, warnings |
| `1` | normal (default) | SCF, TRAH (native and OpenTRAH), Davidson, Z-vector, GMRES and CC iteration tables; one convergence summary per CPHF solve; orbital energies; SCF energy components; DFT grid statistics |
| `2` | detailed | MO coefficients and the orbital table of an unconverged SCF; the primitive-by-primitive basis listing; solver diagnostics (per-right-hand-side CPHF residuals, DFT XC integration timings, Hessian response residuals and storage notes, the NAC overlap table, MOM reordering, NMR gates); notes on symmetry reductions skipped by design; the dispersion block when dispersion is off |
| `3` | debug | Developer dumps: PCM, spin-orbit, scalar-relativistic and MRSF debug output; OpenTRAH MINRES internals |

`runtype = md` and `namd` run at level 0 unless the input sets a level above 1.

At every level the module banners and step timings are written. The LibXC
header, the DFT grid description and each functional's description and
literature references are written once per run rather than at every SCF,
response, gradient or Hessian step, and again only when the grid or exchange
parameters change. Dispersion settings and dispersion-corrected energies appear
only when dispersion is requested.

## Compatibility policy

The following markers remain stable:

- `PyOQP started at`, `PyOQP build:`, and `PyOQP terminated at`;
- `PyOQP method:`, `PyOQP electronic energies`, and
  `PyOQP electronic gradients`;
- `PyOQP state`, CASSCF/FCI field names, and physical MRSF state labels;
- native `SCF`, Davidson, Z-vector, CPHF, CC, and displacement iteration
  markers at the default verbosity (level 0 omits the iteration tables);
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
