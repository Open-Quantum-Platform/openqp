# Regression References and the Test Registry

OpenQP's example suite (`openqp --run_tests all`) is also its regression
harness: each same-stem `.inp`/`.oqp` pair shares one committed `.json`
reference against which its numbers are checked. What gets checked is defined
by a single **registry**, so the suite stays honest as features are added.

## The registry (single source of truth)

`pyoqp/oqp/utils/regression.py` declares every regression quantity once as a
`RegKey` row — its JSON key, the runtypes/methods/properties it applies to,
whether it is required, where it lives (primary JSON or the `.hess.json`
sidecar), and flags for sign/phase ambiguity. That one registry drives three
things so they cannot drift apart:

- **comparison** — `check_ref` compares exactly the registered keys for the run;
- **lean references** — `save_data(lean=True)` keeps only registry keys (plus
  identity/metadata), dropping the large internal `OQP::` arrays;
- **the gate** — `openqp --validate_examples` fails if any example reference is
  missing a value its runtype requires.

## Adding a new regression value

This is the one place you touch:

1. Add a `RegKey(...)` row in `regression.py` (key, `runtypes`, `required`,
   `source`, and any `needs_excited` / `needs_prop` / `phase_invariant` gates).
2. If it is a newly computed quantity, emit it under the same key from
   `Molecule.get_results()` (and, if it originates in Fortran, store it to a
   tagarray as in `OQP::nmr_shielding`).
3. Regenerate the affected references (see below). `lean_keep` keeps the new key
   automatically.

From then on the value is compared, kept in lean references, and **enforced**:
an example missing it fails the gate.

## Generating and validating references

Use the supported commands — do not hand-edit references:

```bash
# (Re)write lean references for one or more inputs (registry keys only)
openqp --generate_reference path/to/example.inp [more.inp ...]

# Verify every reference carries the values its runtype requires
openqp --validate_examples            # defaults to $OPENQP_ROOT/share/examples
```

`--generate_reference` runs each input and writes only the registry keys next to
the `.inp`, so committed references stay small and consistent. It validates what
it wrote before returning.

## What is enforced in CI

- `openqp --validate_examples` — every reference is complete (blocking).
- `openqp --run_tests all` — the numbers still match (blocking; the step uses
  the default `auto` input selection and `set -o pipefail`, so a failing example
  fails the build). Within the standard suite, `auto` prefers eligible `.oqp`
  decks and adds a representative legacy `.inp` compatibility set without
  doubling the entire suite.
- An agent review flags a new computed value added without a `RegKey`, an
  example added without a reference, and code/manual drift.

Maintainers can select the syntax within a test scope when diagnosing input
handling:

```bash
openqp --run_tests all --input-format inp
openqp --run_tests all --input-format oqp
openqp --run_tests all --input-format both
```

With `all`, these selectors retain the standard suite's slow and
non-self-contained exclusions. Pass an explicit examples directory when every
matching file below that directory should be included. The runner gives every
calculation an isolated output directory, and paired legacy jobs also receive
distinct project and log names, so parallel `.inp`/`.oqp` artifacts cannot
overwrite one another.

## Coverage notes

- **Hessian** references compare `freqs`, `infrared_intensities` (IR) and
  `raman_activities` (Raman) from the `.hess.json` sidecar, not just the matrix.
- **NAC** compares the NACME derivative-coupling matrix by magnitude
  (sign/phase ambiguous between builds).
- **IRC** is regression-tested on energy only; the final IRC geometry is a
  non-stationary path point that is not reproducible across BLAS/compilers.
- **Properties** (`dipole`, `mulliken_charges`, `lowdin_charges`,
  `resp_charges`) are opt-in via `[properties] scf_prop` and tested only when
  requested.
