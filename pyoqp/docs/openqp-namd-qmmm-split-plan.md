# Splitting NAMD + QM/MM into `openqp-namd-qmmm`

Status: **plan** · this document motivates the incremental work; **step 1 (compile-time
gating) lands with this commit**. Target external package:
`github.com/Open-Quantum-Platform/openqp-namd-qmmm`.

## 1. Goal and shape

Extract the nonadiabatic-dynamics (NAMD/SOC-NAMD) and ESPF QM/MM features into **one
external Python package**, `openqp-namd-qmmm`, that sits on top of an unchanged compiled
core. The package has two internal submodules:

```
openqp-namd-qmmm/                    (import: openqp_namd_qmmm)
├── dynamics/   ← namd.py            (NAMD, NAMD_SOC, NAMD_SOC_MCH, *_QMMM)
└── qmmm/       ← qmmm_driver.py, qmmm_md.py, qmmm_connectivity.py, MMBackend
```

`openmm` is an **optional extra**: `pip install openqp-namd-qmmm` gives gas-phase
NAMD/SOC-NAMD; `pip install openqp-namd-qmmm[openmm]` adds ESPF QM/MM. The
`dynamics → qmmm` dependency is one-directional and lazy (only the `*_QMMM` classes
import it), so a later split into a standalone `openqp-qmmm` wheel stays cheap.

## 2. Why complete native separation is impossible (and undesirable)

The ESPF QM/MM coupling is **not** a self-contained module — it is fused into the
one-electron / SCF / gradient path in compiled Fortran, gated at runtime by
`control%qmmm_flag`:

| Native hook | used by (core) | Layer |
|---|---|---|
| `qmmm.F90` (`qmmm_mod`) | `int1e.F90:36` (adds `Hqmmm`→`hcore`), `scf.F90:71`, `hf_gradient.F90:154`, `tdhf_mrsf_z_vector.F90:1093` | 1e Hamiltonian / SCF / gradient / z-vector |
| `int1.F90` `omp_qmmm`, `grd1.F90` | 1e integral / gradient assembly | 1e integrals |
| `soc_mrsf.F90` | `runtype=soc` single point **and** SOC-NAMD; pulls the 2e engine | shared electronic structure |
| `namd.F90` (`mrsf_namd_hop`) | NAMD-only, but still links `types`/`c_interop`/`tagarray` | compiled in `liboqp` |

So the split is a **Python-layer extraction on an unchanged core**. Physically moving the
Fortran out of the monolithic `liboqp` buys nothing (it would still link `liboqp`) and
`soc_mrsf.F90`/`qmmm.F90` cannot leave the SCF path at all. The compute-intensive work in
an *ab-initio* NAMD step is ~100% electronic structure (SCF/MRSF/gradient/SOC/NAC/ESPF —
already Fortran) plus the `mrsf_namd_hop` kernel (already Fortran); the Python layer only
sequences these on tiny arrays, so it stays Python with no measurable cost.

## 3. Module inventory

**Move → `openqp-namd-qmmm`** (pure Python; imports *from* core only, never the reverse):

| Module | → submodule | Contents |
|---|---|---|
| `pyoqp/oqp/library/namd.py` | `dynamics/driver.py` | the 6 NAMD classes, velocity-Verlet, adaptive dt, TDC, SHARC local-diabatization, EDC, phase tracking, SHAKE/RATTLE |
| `compute_namd` (`runfunc.py:26-53`) | `dynamics/runfunc.py` | `runtype=namd` dispatch |
| `[md]` schema group (`oqpdata.py`, 24 keys) | `dynamics/schema.py` | NAMD-only config |
| `qmmm_connectivity.py` | `qmmm/connectivity.py` | link atoms, RCD/RC/Z1 frontier charges (numpy) |
| `qmmm_driver.py` | `qmmm/driver.py` + `qmmm/backends/openmm.py` | ESPF coupling, becomes `OpenMMBackend` |
| `qmmm_md.py` | `qmmm/md.py` | ground-state `QMMM_MD` (`runtype=md`) |
| `[qmmm]` schema keys | `qmmm/schema.py` | QM/MM config |
| tests / examples / `soc_namd_plan.md` | `tests/ examples/ docs/` | NAMD + SOC-NAMD-QMMM |

**Stays in core, consumed as the stable API:** native `oqp.*` kernels (`mrsf_namd_hop`,
`soc_mrsf`, `espf_op_corr`, `form_esp_charges`, `grad_esp_qmmm[_excited]`,
`get_states_overlap`); the five `single_point.py` drivers (`SinglePoint`, `Gradient`,
`LastStep`, `BasisOverlap`, `NACME`); the `Molecule`/`OQPData` accessor surface;
`dump_log`.

**Delete (dead legacy):** `libopenmm.py` and `runfunc.compute_md` (calls undefined
`qmmm_md.run_md()`).

**The one coupling to break:** the eager `from oqp.library.runfunc import (…, compute_namd)`
at `pyoqp.py:114`. Everything else is already lazy. Replaced by lazy/plugin registration so
`runtype=namd`/`md` without the package installed fails with a clean *"pip install
openqp-namd-qmmm"* message.

## 4. Compile-time gating (this PR — step 1)

Two opt-in CMake options, **default `ON`** so a stock build is byte-for-byte unchanged; they
exist so the external package can require a core built with the kernels, and so a lean
energy/opt-only core can drop them.

| Option | Default | Gates |
|---|---|---|
| `ENABLE_QMMM` | `ON` | compiles `qmmm.F90`; defines `OQP_ENABLE_QMMM`; guards the in-core ESPF hooks |
| `ENABLE_NAMD` | `ON` | compiles `namd.F90` (`mrsf_namd_hop`); **no in-core guards needed** (nothing does `use namd_mod`) |

SOC is intentionally **not** gated: `soc_mrsf.F90` backs `runtype=soc` (a standalone single
point) as well as SOC-NAMD, so it always compiles. (A separate `ENABLE_SOC` is possible but
a footgun — omitted unless a SOC-free core becomes a hard requirement.)

**Mechanism** (mirrors the existing `OQP_HAVE_OPENTRAH` / `OQP_ENABLE_DDX` pattern):
`option()` in the top-level `CMakeLists.txt`; `list(REMOVE_ITEM …/qmmm.F90|namd.F90)` from
the module GLOB when OFF (`source/modules/CMakeLists.txt`); `target_compile_definitions(oqp
PRIVATE OQP_ENABLE_QMMM|OQP_ENABLE_NAMD)` (`source/CMakeLists.txt`). Six `#ifdef
OQP_ENABLE_QMMM` guards wrap every `qmmm_mod` reference:

- `scf.F90:71` (`use`) and the `add_potqm_contributions` `select case` (the `hcore =
  hcore_bk + dhcore` line stays outside — `dhcore` is `source=0.0_dp`, so it's a correct no-op)
- `int1e.F90:36`, `hf_gradient.F90:154` (`use` lines; their calls are already commented out)
- `tdhf_mrsf_z_vector.F90:1092-1093` (`use`) and the `if(qmmm_flag)` excited-ESPF block

The `qmmm_flag`/`soc_2e` bind(c) struct fields and the tagarray tag-name constants are
**never** guarded (C-ABI / registry invariants).

**Build lines**

```bash
pip install .                                    # stock: QMMM+NAMD in, SOC in (default)
pip install . -C cmake.define.ENABLE_QMMM=OFF -C cmake.define.ENABLE_NAMD=OFF   # lean core
```

| Build | `qmmm.F90` | `namd.F90` | `runtype=energy/grad/soc` | `runtype=md`/`namd` |
|---|---|---|---|---|
| Default (`ON`/`ON`) | ✓ | ✓ | ✓ | ✓ (with package) |
| `OFF`/`OFF` | ✗ | ✗ | ✓ | package import probe fails cleanly |

The package adds an import-time capability probe (`hasattr(oqp.lib, "mrsf_namd_hop")`, …) so a
mismatched core errors early with the fix instead of a deep symbol failure mid-run.

## 5. General MM-backend interface (step 2)

**Load-bearing decision: all ESP/ESPF/coupling/frontier math stays in the generic layer; a
backend supplies only `(positions, MM charges, MM energy, MM forces, box)`.** This is already
how `qmmm_driver.py` is factored (pure NumPy + `oqp.espf_op_corr`, no OpenMM), so every other
engine gets correct QM-MM electrostatics for free.

```python
class MMBackend(ABC):
    def build_system(topology, forcefield, qm_atoms, cutoff_method, *, rigid_water=False)
    def set_positions(x_nm)
    def get_mm_charges() -> (q, xyz_bohr, idx)      # the entire ESP-embedding contract
    def set_qm_charges(qm_atoms, charges_or_none)   # ESP charges (split) or zero (full-ESPF)
    def mm_energy_forces() -> (E_hartree, F_hartree_per_bohr)
    def periodic_box_bohr() -> ndarray|None
    def get_natoms/get_masses/iter_bonds/element_of(...)
    # optional: supports_pme(), potmm_from_engine(qm_xyz), get_constraints_bohr()
```

- `OpenMMBackend` = today's `OpenQpQMMM` re-expressed = reference impl (behind the `openmm` extra).
- `NativeCoulombBackend` (pure NumPy) = first non-OpenMM backend, validates the interface.
- PME is the only engine-specific piece: the split-scheme POTMM/POTQM diffs OpenMM Ewald
  energies. **Full-ESPF real-space is the portable default**; PME is an optional
  `potmm_from_engine()` override that only OpenMM implements.
- Selection via a **`[qmmm] mm_engine = "openmm"` config key** through a small dict registry
  (parallel to `cutoff`/`embedding`/`frontier_scheme`), with an entry-point fall-through for
  out-of-tree backends. Amber/GROMACS/LAMMPS/Tinker/ASE all feasible for the full-ESPF path.
- Two roles kept separate: **MM-as-force-provider** (what `MMBackend` covers; portable) vs
  **MM-as-propagator** (`qmmm_md.py` drives MD inside OpenMM — stays an OpenMM-only extra).

## 6. Migration order — core changes land first

1. **Compile-time gating** (this PR), defaults `ON` → symbol-identical to today, validated by
   the full test suite.
2. **MM-backend refactor in-place**: introduce `MMBackend`, wrap `OpenQpQMMM` as
   `OpenMMBackend`, add `NativeCoulombBackend`, rewire `namd.py`/`qmmm_md.py` to
   `make_backend(...)`.
3. **Break the Python seam**: replace the eager `compute_namd` import with
   lazy/entry-point dispatch; add the package import probe.
4. **Extract** `dynamics/` + `qmmm/` into `openqp-namd-qmmm` with `openmm` as an extra.

Rationale: the guards and the backend seam are the two places a bug silently corrupts
scientific results (double-counted electrostatics, a missing kernel), so they are proven in
the monorepo before a package boundary can hide a regression.

## 7. Risks / open questions

- **cffi dlopen vs. hard-linked `_oqp`.** The `hasattr` probe is correct for the lazy dlopen
  path (`oqp/__init__.py:42`). If the package ships against the compiled `_oqp` extension,
  unconditional prototypes in `include/oqp.h` for absent symbols would fail at link time
  instead — decide whether those prototypes need `#ifdef` gating too.
- **Split scheme is OpenMM-only.** Consider having the config validator reject
  `mm_engine != openmm` combined with the split scheme (`espf_full=False`) to prevent a
  silently-wrong run.
- **ABI invariants.** `types.F90` `qmmm_flag`/`soc_2e` and the `tagarray_driver.F90` tag
  constants must never be `#ifdef`'d — guarding them desyncs the C struct/registry.
- **GLOB staleness.** `REMOVE_ITEM` needs a reconfigure; options are cache vars so toggling
  re-triggers it, but a stale build dir may need `--fresh`.
