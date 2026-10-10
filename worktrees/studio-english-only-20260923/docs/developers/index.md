# Developer Guide

This section is for contributors extending OpenQP itself. User-facing scripts
should start with the [Python API](../api/python-runner.md); contributor work
usually touches the input schema, input checker, Python workflow layer, native
kernels, tagarray data records, tests, and manual pages together.

## Where To Start

| Task | Start here |
| --- | --- |
| Add or change a user keyword | [Contributing to OpenQP](contributing.md), then the relevant keyword page. |
| Add a Python workflow helper | `pyoqp/oqp/openqp.py`, tests, and the Python API page. |
| Add a native kernel or binding | [Native Interface](native-interface.md). |
| Store data across Python and Fortran | [Data Tags and Tagarray](tagarray.md). |
| Add result accessors | [Results and Molecule Data](../api/results.md) plus the tagarray page. |

## Contributor Principles

- Keep ordinary user setup in `molecule`, `theory`, `workflow`, `control`, and
  `settings`.
- Add compatibility paths only when they protect existing inputs or scripts.
- Validate unsupported method/workflow combinations before dispatching kernels.
- Document new user-facing behavior in the manual, not only in generated
  Sphinx or FORD output.
- Add focused tests near the Python helper, input checker, or native module
  that owns the behavior.
