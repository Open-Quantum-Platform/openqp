# `[tests]`

The `[tests]` section is an internal regression-harness control. It is accepted
by the input schema so example decks can describe an expected failure.

### `exception`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | OpenQP example/regression harness |

When true, the harness expects the calculation to raise rather than finish
normally. Production calculation inputs should omit this section.
