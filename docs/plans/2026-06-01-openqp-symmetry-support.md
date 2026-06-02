# OpenQP symmetry support autopilot plan

Scope: metadata/label support only. Do not enable integral or response-space reductions by default. Preserve C1 behavior and GAMESS NOSYM=1 compatibility. Keep GPU/METC/libintx work out of this branch.

Current gates:
- [x] Parser/checker coverage for symmetry metadata flags.
- [x] Metadata persistence tests.
- [x] Backend-free one-electron block/projector tests.
- [x] Cleanly separate the symmetry work from any GPU branch context.
- [x] Run focused pytest for `tests/test_symmetry_parser_checker.py`, `tests/test_symmetry_metadata.py`, and `tests/test_symmetry_one_electron_blocks.py`.
- [x] Remove or document any blocker under `/Volumes/External_Storage/Hermes-Agent/openqp-symmetry-autopilot/`.

Next increment rule: take exactly one unchecked item, run RED/GREEN or focused verification, then stop with a compact final report.
