# OpenQP Documentation

This repository contains the user manual for
[OpenQP](https://github.com/Open-Quantum-Platform/openqp), the Open Quantum
Platform.

The manual covers the OQP Studio desktop application, installation, input
files, Python usage, build options, capabilities, examples, keyword references,
and workflows including HF/DFT, MP2, TDHF/TDDFT, SF-TDDFT, MRSF-TDDFT,
PCM/ddX, SOC, NACME, EKT, Hessians, optimization, molecular symmetry, and
spectroscopy-related properties.

Manual site:
[https://open-quantum-platform.github.io/openqp-docs/](https://open-quantum-platform.github.io/openqp-docs/)

API guide:
[https://open-quantum-platform.github.io/openqp-docs/api/](https://open-quantum-platform.github.io/openqp-docs/api/)

## Preview Locally

```bash
pip install -r docs/requirements.txt
mkdocs serve
```

## Build

```bash
mkdocs build --strict
```

## Maintenance

Development is maintained in the Open Quantum Platform group's GitLab project.
Use GitLab branches and merge requests for new documentation changes; its
`mkdocs build --strict` pipeline must pass before merging.

GitHub remains the public mirror and hosts the published manual. Updates to
GitHub `main` automatically create a GitLab synchronization merge request when
needed. That request is merged only after GitLab CI succeeds. Conflicts remain
visible for manual resolution. Updates from GitLab to GitHub are made on demand,
after checking that `main` is suitable for public release. Automatic outbound
mirroring stays disabled; private branches and archived pull-request refs must
not be published. The existing GitHub Pages workflow publishes the manual when
GitHub `main` is updated.

Keep keyword and API pages aligned with the OpenQP input schema, checker, and
Python entry points:

- [`oqpdata.py`](https://github.com/Open-Quantum-Platform/openqp/blob/main/pyoqp/oqp/molecule/oqpdata.py)
- [`input_checker.py`](https://github.com/Open-Quantum-Platform/openqp/blob/main/pyoqp/oqp/utils/input_checker.py)
- [`single_point.py`](https://github.com/Open-Quantum-Platform/openqp/blob/main/pyoqp/oqp/library/single_point.py)
- [`pyoqp.py`](https://github.com/Open-Quantum-Platform/openqp/blob/main/pyoqp/oqp/pyoqp.py)
- [`openqp.py`](https://github.com/Open-Quantum-Platform/openqp/blob/main/pyoqp/oqp/openqp.py)
- [`molecule.py`](https://github.com/Open-Quantum-Platform/openqp/blob/main/pyoqp/oqp/molecule/molecule.py)
- [`oqp.h`](https://github.com/Open-Quantum-Platform/openqp/blob/main/include/oqp.h)
- [`CMakeLists.txt`](https://github.com/Open-Quantum-Platform/openqp/blob/main/CMakeLists.txt)
- [`pyproject.toml`](https://github.com/Open-Quantum-Platform/openqp/blob/main/pyproject.toml)

When OpenQP changes input defaults, allowed values, workflow support, or Python
runner behavior, update the matching page under `docs/keywords`,
`docs/workflows`, `docs/api`, `docs/python-scripting.md`, or
`docs/build-options.md`. Also update this README when a new workflow, keyword
section, or maintenance source file becomes part of the documented surface.
