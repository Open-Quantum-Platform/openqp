# QM/MM ↔ tight-binding reconciliation (full-ESPF covalent boundary for `method=dftb` and `method=xtb`)

This note records the merge of `origin/codex/soc-namd-options` (the correct
full-ESPF covalent-boundary QM/MM) into the TB-adapter branch
(`feat/openqp-xtb-adapter`: DFTB adapter PR #266 + the `method=xtb` adapter + the
`oqp.utils.tb_backends` dispatch generalization). The goal is a single tree where
the full-ESPF QM/MM driver internals coexist with the TB dispatch so that
`method=dftb` **and** `method=xtb` QM/MM both conserve energy at covalent
(hydrogen link-atom) boundaries.

Merge-base: `55befcb8f` (PR #253). `soc-namd-options` was 83 commits ahead in
QM/MM; the TB branch was 13 commits ahead. The merge was performed
`soc-namd-options → feat/xtb-qmmm-fullespf`.

## What full-ESPF brings

The released ESPF path embedded covalent boundaries with **frontier-charge
redistribution** (`frontier_scheme = none|rcd|rc|z1`, PR #258): the M1 host
charge could be deleted/redistributed onto virtual midpoint charges. `soc-namd-
options` replaced that with a cleaner, gradient-exact **full-field ESPF** scheme
and removed the redistribution machinery entirely:

- `OpenQpQMMM.__init__` builds a per-QM-centre MM-exclusion structure
  (`_build_qm_mm_exclusions`) instead of a frontier-host redistribution map. The
  embedding potential (`_full_field_potmm`) and the coupling force
  (`_coupling_forces`) use the *same* per-centre view of the MM charges, so the
  analytic gradient stays exact across a covalent cut.
- The full-ESPF native path adds the **nuclear-MM interaction** `Σ_A Z_A φ_A`
  to complete the QM–MM electrostatic energy (commit `fac2be89c`, "nuclear-MM
  energy → FD-exact covalent-boundary forces"), and drops the frontier
  redistribution / virtual-site scatter.
- `espf_full` is now selected by `Embedding ∈ {espf, espf_full, electrostatic}`
  and is the default.
- PBC / min-image support (`_box_lengths_bohr`, `_min_image`) is retained.

The `redistribute_frontier_charges` / `assemble_embedding_sites` helpers,
`_frontier_hosts` / `_embedding_sites`, the `frontier_scheme` schema key and API
parameter, and the frontier tests/examples were all removed on `soc-namd-options`;
this merge takes that removal.

## TB dispatch contract (unchanged by full-ESPF)

For `method=dftb`/`xtb` the openqp-dftb/openqp-xtb library folds the per-QM-centre
MM potential (Hartree/e) **directly into the SCC Hamiltonian** via
`mol.dftb_external_potential = POTMM`:

- the returned state energy is the **complete** embedded QM energy (it already
  contains `E_ext = Σ_A q_A φ_A` with `q_A` the *net* atomic charge), so the
  driver must **not** add the native nuclear `Σ_A Z_A φ_A` term for TB;
- the returned analytic gradient is `d(E_embedded)/dR_QM` at fixed potential (the
  Pulay/charge-response coupling is already inside it), so the native
  `espf_op_corr` hcore mutation and `grad_esp_qmmm(_excited)` / `OQP::ESPF_GRAD`
  additions are **skipped**;
- the adapter publishes relaxed net atomic charges of the active state into
  `OQP::partial_charges`; the driver's backend-agnostic `_coupling_forces` reads
  those for the classical `dφ/dR` MM-field forces (the DFTB analog of
  `form_esp_charges`).

Only full-ESPF (`espf_full`) or mechanical (`potmm is None`) embedding is
supported for TB; the legacy `split` scheme is rejected (it would double-count
the coupling already inside the embedded TB energy).

## TB-dispatch call sites re-applied on the full-ESPF driver

Dispatch helpers live in `pyoqp/oqp/utils/tb_backends.py`:
`is_tb_method`, `tb_section_name`, `tb_config`, `make_tb_adapter`.

`pyoqp/oqp/library/qmmm_driver.py` (`OpenQpQMMM`):

- import `is_tb_method`;
- `_forces_qm`, **mol mode** — after `_update_mol_positions()`, before the native
  `SinglePoint`/`espf_op_corr` path: `if is_tb_method(...): return
  self._forces_qm_dftb(self.mol, potmm)`;
- `_forces_qm`, **config mode** — after `self.op = OPENQP(...)`, same early
  dispatch on `self.op.mol`;
- new method `_forces_qm_dftb(mol, potmm)`: sets `mol.dftb_external_potential`,
  runs `Gradient(mol).gradient()` (dispatches to the TB adapter), reads the active
  state energy/gradient and `OQP::partial_charges`; guards `espf_full` (else
  `NotImplementedError`).

`pyoqp/oqp/library/namd.py`:

- import `is_tb_method, make_tb_adapter, tb_section_name`;
- `NAMD_QMMM._electronic_structure` — TB branch sets `mol.dftb_external_potential`
  and runs the embedded SCF/excitation without the ESPF hcore mutation
  (covalent-boundary electronic step);
- `NAMD_QMMM._forces_qm` — TB early-return of the raw gradient (skips
  `grad_esp_qmmm`/`OQP::ESPF_GRAD`);
- `NAMD_QMMM._total_force` — TB energy uses `mol.energies[active]` with **no**
  nuclear `Σ Z_A φ_A` term;
- `NAMD_SOC`, `NAMD_SOC_MCH`, `NAMD_SOC_QMMM`, `NAMD_SOC_MCH_QMMM` — TB SOC arms:
  `tb_section_name(...)['target_multiplicity']` bookkeeping, `_dftb_soc_tags`
  (one-centre SOC + numpy `eigh`) in place of `oqp.soc_mrsf`, and
  `_dftb_spatial_overlap` in place of `oqp.get_states_overlap`;
- module helpers `_dftb_soc_tags`, `_dftb_spatial_overlap` (TB, dftb+xtb via
  `make_tb_adapter`).

`pyoqp/oqp/library/runfunc.py` — `compute_soc` TB branch (`xtb_soc`/`dftb_soc`);
`compute_namd` (soc's QM/MM NAMD dispatch) retained.

`pyoqp/oqp/library/single_point.py` — TB dispatch at `energy`, `reference`,
`excitation`, `gradient`, and the `dftb_overlap`/`dftb_states_overlap` state-
overlap paths (all via `is_tb_method`/`make_tb_adapter`); soc's
`import oqp.utils.qmmm` and `scf_grad` tweaks retained.

`pyoqp/oqp/utils/input_checker.py` — NAMD gate allows
`method ∈ {tdhf, dftb, xtb}`; the `method=dftb/xtb` QM/MM validator (rejects
`split`, requires full-ESPF electrostatic, gates SOC-QMMM off as "not yet wired")
retained from the TB branch and now describes the full-ESPF scheme.

`pyoqp/oqp/molecule/oqpdata.py` — `[dftb]`/`[xtb]` schema sections retained;
`frontier_scheme` removed.

`pyoqp/oqp/openqp.py` — `[dftb]`/`[xtb]` builder API retained; `frontier_scheme`
argument removed from `qmmm(...)`.

`pyoqp/oqp/utils/oqp_tester.py` — openqp-dftb/openqp-xtb "backend not found →
SKIPPED" handling and the `method=dftb`/`method=xtb` example-skip predicate
retained.

## Conflict resolution per file

| File | Type | Resolution |
| --- | --- | --- |
| `pyoqp/oqp/library/namd.py` | add/add | Take soc (full-ESPF), re-apply the TB delta (`674aab601..HEAD`), which applied cleanly (soc namd.py ≡ released namd.py apart from the 2-line `frontier_scheme` removal). |
| `pyoqp/oqp/library/qmmm_driver.py` | add/add | Take soc (full-ESPF); re-apply the TB delta — 3/4 hunks clean, the `tb_backends` import added manually (soc had edited the adjacent `qmmm_connectivity` import). |
| `pyoqp/oqp/library/qmmm_connectivity.py` | add/add | Take soc (frontier helpers removed). TB delta never touched this file. |
| `pyoqp/oqp/library/qmmm_md.py` | add/add | Take soc (removes `frontier_scheme` + the internal `runtype→energy` forcing). TB dispatch is in the driver, not here. |
| `pyoqp/oqp/utils/qmmm.py` | (auto) | Identical on both sides; auto-merged. |
| `pyoqp/oqp/utils/input_checker.py` | content | Take HEAD for the NAMD `method ∈ {tdhf,dftb,xtb}` gate; rest auto-merged (soc QM/MM checks + TB validator both present). |
| `pyoqp/oqp/molecule/oqpdata.py` | content | Take HEAD `timestep` (`float`, safer than soc's `int` — the MD/NAMD code casts to float and `int` truncates sub-fs steps); take soc (remove `frontier_scheme`). `[dftb]`/`[xtb]` sections auto-merged. |
| `pyoqp/oqp/openqp.py` | content | Take soc (remove `frontier_scheme` param/docstring/update). `[dftb]`/`[xtb]` API auto-merged. |
| `pyoqp/oqp/utils/oqp_tester.py` | content | Take HEAD (openqp-dftb/xtb SKIPPED handling + skip predicate). |
| `pyoqp/oqp/library/runfunc.py` | (auto) | Auto-merged: soc `compute_namd`/`QMMMOpt` + TB `compute_soc` branch both present. |
| `pyoqp/oqp/library/single_point.py` | (auto) | Auto-merged: soc `import qmmm`/`scf_grad` + TB dispatch both present. |
| `include/oqp.h`, `source/tagarray_driver.F90` | (auto) | Auto-merged; soc additions verified present. |
| `tests/test_openqp_api.py` | content | Take soc (drop `frontier_scheme` schema + `test_qmmm_frontier_scheme_sets_section_key`); TB tests auto-merged. |
| `tests/test_rohf_status_and_interface.py` | content | Take HEAD (the `tb_backends` stub bootstrap the shared `try` block depends on). |
| `README.md`, `examples/QMMM/README.md` | content/add | Take HEAD tutorial links in README; take soc for the examples README (drop the `frontier_scheme` section). |
| `tests/test_qmmm_frontier.py`, `tests/test_qmmm_frontier_openmm.py`, `examples/QMMM/ala-dipeptide_BHHLYP-QMMM-MD-RCD.inp` | stale | **Removed** — HEAD-only artifacts of the removed `frontier_scheme` feature; they import `redistribute_frontier_charges`/`assemble_embedding_sites` (deleted on soc) and would fail to load. Absent on soc. |

Adapter files `pyoqp/oqp/library/openqp_dftb.py`, `openqp_xtb.py`,
`pyoqp/oqp/utils/tb_backends.py` are TB-branch-only (no conflict) and kept as-is.

## Static verification (this agent)

- Every changed `.py` file parses (`ast.parse`).
- No conflict markers remain in the tree.
- Dependency-free unit tests (no compiled `liboqp`): the DFTB/XTB schema-hook
  tests, the runfunc/interface dispatch test, and the QM/MM builder API test pass
  — **54 passed, 20 skipped** (skips are OpenMM/liboqp-gated).
- Full `tests/` run vs the pre-merge HEAD baseline: **identical set of 51
  pre-existing failures** (all require the compiled `oqp` extension, which is not
  built in this environment) — the merge introduces **no new failures**. The
  drop in passing count equals exactly the removed frontier tests.

## Remaining runtime-validation checklist (needs the integrated runtime)

The compiled `_oqp`/`liboqp` extension plus the external openqp-dftb / openqp-xtb
libraries are **not** built in this agent's environment, so the numeric
covalent-boundary behaviour was not exercised. Once the integrated runtime is up:

1. **Covalent-boundary energy-vs-force FD, `method=dftb`** — a QM/MM partition
   that cuts a covalent bond (H link-atom cap), full-ESPF electrostatic
   embedding, `runtype=grad`: verify `F ≈ −dE/dR` by finite differences on both
   QM and frontier/host MM atoms (this is the original DFTB QM/MM bug this merge
   targets).
2. **Same FD check for `method=xtb`.**
3. **Whole-molecule QM/MM regression** (no link atom) for dftb and xtb — confirm
   full-ESPF still matches the pre-merge energies/forces.
4. **Ground-state QM/MM MD** (`runtype=md`, `QMMM_MD`) energy conservation for
   dftb and xtb across a covalent boundary.
5. **NAMD-QMMM** (`runtype=namd`, non-SOC, `embedding=electrostatic`) for dftb —
   the `NAMD_QMMM._electronic_structure`/`_forces_qm`/`_total_force` TB arms.

### Honest flags for the runtime tests

- **`_forces_qm_dftb` is the released method name/body applied onto soc's driver.**
  Its `espf_full` guard and its reliance on `_coupling_forces` reading
  `OQP::partial_charges` are consistent with soc's `_coupling_forces` (verified
  by reading both), but the *numerical* interaction of the TB net-charge coupling
  with soc's **new full-field exclusion structure** (`_build_qm_mm_exclusions`,
  which in practice applies *no* exclusions — see its NOTE) has not been run.
  Item 1/2 above is the authoritative check.
- **SOC-QMMM for TB.** `namd.py` carries TB SOC-QMMM arms
  (`NAMD_SOC_QMMM`/`NAMD_SOC_MCH_QMMM`), but `input_checker` still gates
  `dftb/xtb` SOC-QMMM off as "not yet wired". The code path exists but is not
  reachable through the checker; do not assume it is validated.
- **`qmmm_md.py`** dropped soc's removal of the internal `runtype→energy`
  forcing block. The TB QM/MM MD path routes through
  `OpenQpQMMM._forces_qm` (which has the TB dispatch), so this is expected to be
  fine, but a `runtype=md` dftb/xtb smoke run (item 4) should confirm the QM
  sub-run accepts the config.
