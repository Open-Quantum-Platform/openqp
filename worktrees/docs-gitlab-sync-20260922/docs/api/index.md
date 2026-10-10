# API Documentation

OpenQP exposes several layers of programmatic access. Most users should start
with the Python API because it covers both compact `OpenQP` scripts and direct
`Runner` execution of existing input files.

## API Layers

| Layer | Entry point | Audience | Stability |
| --- | --- | --- | --- |
| Command line | `openqp input.oqp` | Production runs, shell scripts, workflows | Recommended |
| Legacy command line | `openqp input.inp` | Existing sectioned decks and legacy-only controls | Stable |
| Python API | `oqp.openqp.OpenQP`, `oqp.pyoqp.Runner` | Scripts, notebooks, services, and automated workflows | Recommended |
| Molecule data | `runner.mol` | Result extraction and advanced workflows | Semi-stable |
| MRSF analysis and interoperability | `oqp.interop` | NTOs, attachment/detachment densities, cubes, QCSchema, FCIDUMP, and external comparisons | Recommended |
| Input checker | `oqp.utils.input_checker.check_input_values` | Input builders, web forms, preflight validation | Recommended |
| Convenience wrapper | `oqp.openqp.OPENQP` | Compact dotted-key Python inputs | Experimental |
| C/Fortran handle | `include/oqp.h`, `oqp_handle_t` | OpenQP developers and native extensions | Internal |

The import package is `oqp`, even though the source tree contains a `pyoqp`
directory and the installable distribution is named `OpenQP`.

## Recommended Starting Points

- Use [Run OpenQP from Python](../python-scripting.md) for complete script-based
  calculation examples.
- Use [Python API](python-runner.md) for high-level `OpenQP` scripts, direct
  `Runner` execution, and the supported Python helper functions.
- Use [Results and Molecule Data](results.md) to retrieve energies, gradients,
  orbitals, TDDFT arrays, SOC data, and serialized result dictionaries.
- Use [MRSF Analysis and Interoperability](../workflows/mrsf-analysis.md) for
  post-run MRSF excited-state analysis, exports, and external-code comparisons.
- Use [Input Validation](input-validation.md) when building forms, agents, or
  custom front ends that need actionable diagnostics before launching OpenQP.
- Use the [Developer Guide](../developers/index.md) when contributing new
  keywords, workflows, native kernels, data tags, or bindings.

## Generated References

The manual pages above are the curated API guide. The older generated references
remain useful for source browsing:

- [Generated Python API reference](https://open-quantum-platform.github.io/openqp/sphinx/index.html)
- [Generated Fortran API reference](https://open-quantum-platform.github.io/openqp/ford/index.html)

Generated references can expose internal modules that are not intended as stable
user-facing APIs. Prefer the manual pages for supported scripting patterns.
