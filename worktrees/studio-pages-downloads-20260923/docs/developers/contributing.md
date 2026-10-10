# Contributing to OpenQP

Most feature work crosses several layers. A small keyword or workflow normally
needs code, validation, examples, and manual updates.

## Common Change Path

| Layer | Typical files |
| --- | --- |
| Input schema | `pyoqp/oqp/molecule/oqpdata.py` |
| Input parsing and validation | `pyoqp/oqp/utils/input_parser.py`, `pyoqp/oqp/utils/input_checker.py` |
| Python workflow orchestration | `pyoqp/oqp/library/`, `pyoqp/oqp/openqp.py` |
| Native kernels | `source/modules/`, `source/integrals/`, `source/dftlib/` |
| Native bridge | `include/oqp.h`, CMake source lists, CFFI loading |
| Cross-language data | `source/tagarray_driver.F90`, `pyoqp/oqp/molecule/oqpdata.py` |
| Examples | `examples/` in the OpenQP source repository |
| Manual | `openqp-docs/docs/` |
| Tests | `tests/` |

## Adding A User-Facing Keyword

1. Add the keyword to the schema with a safe default.
2. Map it into the native handle if native kernels need it.
3. Add input-checker rules for unsupported combinations.
4. Add or update a runnable example when the workflow is user-facing.
5. Document the keyword in the Keywords chapter.
6. Add focused tests for parsing, validation, and dispatch behavior.

## Adding A Workflow

1. Decide whether it is a real workflow or a theory/settings detail.
2. Add the Python workflow routing in the library layer.
3. Add a high-level helper only when it simplifies a common user path.
4. Validate theory compatibility before expensive work starts.
5. Store reusable outputs in tagarray records when Python/native layers share
   the data.
6. Add a workflow manual page with both input-file and Python examples.

## Validation Commands

Useful focused checks:

```bash
python3 -m unittest tests.test_openqp_api
python3 -m unittest discover -s tests
python3 -m mkdocs build --strict
```

Use the narrowest test that covers the change while developing, then broaden
coverage when the change touches shared dispatch, native bindings, or public
workflow behavior.
