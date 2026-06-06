# OpenQP symmetry support autopilot plan

Scope: metadata/label support only. Do not enable integral or response-space reductions by default. Preserve C1 behavior and GAMESS NOSYM=1 compatibility. Keep GPU/METC/libintx work out of this branch.

Current gates:
- [x] Parser/checker coverage for symmetry metadata flags.
- [x] Metadata persistence tests.
- [x] Backend-free one-electron block/projector tests.
- [x] Cleanly separate the symmetry work from any GPU branch context.
- [x] Run focused pytest for `tests/test_symmetry_parser_checker.py`, `tests/test_symmetry_metadata.py`, and `tests/test_symmetry_one_electron_blocks.py`.
- [x] Remove or document any blocker under `/Volumes/External_Storage/Hermes-Agent/openqp-symmetry-autopilot/`.

Detection/labeling gates (increment 2):
- [x] Geometry-based Schoenflies point-group detection (`symmetry_detect.py`): atom/linear/cubic/axial trees, covered by `tests/test_symmetry_point_group.py` (Kh, Coov/Dooh, C1/Ci/Cs/C2/C2h/C2v, D2h..D6h, D2d/D3d/D3h/D4h, Td/Oh, noisy-geometry tolerance).
- [x] Largest abelian (D2h-family) subgroup resolution with standard orientation and per-operation atom permutations, verified op-by-op.
- [x] Character tables for the eight D2h-family groups (Mulliken labels).
- [x] Symmetry-adapted transform (SALC) builder for Cartesian s/p/d shells with irrep labels; feeds the existing one-electron block diagnostics (`tests/test_symmetry_salc_mo_labels.py`).
- [x] MO irrep labeling via <m|S O|m> characters with 'mixed' fallback.
- [x] Molecule wiring: `initialize_symmetry_metadata` runs detection when enabled/auto (non-fatal, metadata-only), resolves 'auto' requests, honors `strict` mismatch as an error. Reductions remain hard-off.

Remaining (future increments):
- [ ] Label MOs/states/modes in the actual run outputs (log/json) using `assign_mo_irreps`.
- [ ] Pure spherical-harmonic shells (ISPHER) in the SALC builder.
- [ ] Integral-side symmetry reductions behind `use_integral_symmetry` (off by default until production-ready).
- [ ] Response-side reductions behind `use_response_symmetry`.

Next increment rule: take exactly one unchecked item, run RED/GREEN or focused verification, then stop with a compact final report.
