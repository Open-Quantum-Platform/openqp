# Developer Documentation

Developer documentation should live in this repository beside the code it
describes. Keep user-facing behavior in the manual and implementation details in
developer pages.

## Current Pages

- [Release Packaging](../release-packaging.md)

## Suggested Additions

- Native module and C API structure
- Adding a new input keyword
- Adding and validating examples
- Updating the keyword reference after schema changes
- Build-system and external-cache internals

When a code change modifies input validation, defaults, or workflow support, the
same pull request should update the relevant page under `docs/keywords` or
`docs/workflows`.
