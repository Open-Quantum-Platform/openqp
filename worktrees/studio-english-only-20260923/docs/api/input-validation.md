# Input Validation API

OpenQP validates parsed inputs before a calculation starts. The same diagnostic
engine can be called from Python to build input editors, web tools, notebooks,
or agent workflows that need fast feedback.

## Validating a Config

Use `check_input_values` when you already have a typed OpenQP configuration.

```python
from oqp.utils.input_checker import check_input_values

report = check_input_values(
    config,
    raise_error=False,
    emit=False,
)

if not report.ok:
    print(report.to_text())
```

```python
check_input_values(
    config,
    *,
    raise_error=True,
    emit=True,
    input_dir=None,
) -> CheckReport
```

| Argument | Meaning |
| --- | --- |
| `config` | Parsed OpenQP configuration dictionary. |
| `raise_error` | If `True`, validation errors raise `SystemExit(1)`. |
| `emit` | If `True`, diagnostics are printed as text. |
| `input_dir` | Directory used to resolve paths stored relative to the input file. |

## Parsing Then Validating

For raw sectioned dictionaries, use `OQPConfigParser` with
`OQP_CONFIG_SCHEMA`, then pass the typed config into the checker.

```python
from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
from oqp.utils.input_parser import OQPConfigParser
from oqp.utils.input_checker import check_input_values

raw = {
    "input": {
        "system": "\nH 0 0 0\nH 0 0 0.74",
        "basis": "6-31g*",
        "method": "hf",
        "runtype": "energy",
    },
    "scf": {
        "type": "rhf",
    },
}

parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA, allow_no_value=True)
parser.load_dict(raw)
config = parser.validate()

report = check_input_values(config, raise_error=False, emit=False)
print(report.to_text())
```

The parser applies schema defaults before validation, so the report sees the
same effective input that a real OpenQP run would see.

## Diagnostic Objects

`check_input_values` returns a `CheckReport`.

| Attribute or method | Meaning |
| --- | --- |
| `report.ok` | `True` when there are no `ERROR` diagnostics. |
| `report.errors` | List of error diagnostics. |
| `report.warnings` | List of warning diagnostics. |
| `report.diagnostics` | Full diagnostic list, including informational messages. |
| `report.to_text()` | Human-readable report. |

Each diagnostic contains:

| Field | Meaning |
| --- | --- |
| `severity` | Usually `ERROR`, `WARNING`, or an informational severity. |
| `path` | Keyword path such as `input.runtype` or `tdhf.type`. |
| `message` | Short description of the problem. |
| `value` | Current value, when useful. |
| `expected` | Expected value or allowed values, when useful. |
| `action` | Suggested fix. |
| `wiki` | Extra help text for common user-facing errors. |

## Front-End Pattern

For an input builder or web service, keep validation non-fatal and render each
diagnostic next to the relevant field.

```python
report = check_input_values(config, raise_error=False, emit=False)

errors_by_path = {}
for item in report.errors:
    errors_by_path.setdefault(item.path, []).append(item.action or item.message)
```

This pattern avoids launching a calculation when the input is inconsistent, and
it keeps the suggested fix close to the keyword that needs attention.
