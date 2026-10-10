# Python API

OpenQP exposes two Python layers.

| Layer | Import | Use |
| --- | --- | --- |
| High-level wrapper | `from oqp.openqp import OpenQP` | OpenQP-native scripts with molecule, workflow, and section-style keyword setup. |
| Runner | `from oqp.pyoqp import Runner` | Existing input files, explicit sectioned dictionaries, and front ends that already build OpenQP input data. |

Both layers use the same OpenQP input schema and the same native calculation
engine. `OpenQP` builds a sectioned dictionary and passes it to `Runner`.

For a user-level guide with complete scripts, see
[Run OpenQP from Python](../python-scripting.md). This page keeps the
user-facing Python helpers in one place; developer-facing native details are in
the [Developer Guide](../developers/index.md).

## High-Level OpenQP Wrapper

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_mrsf", silent=1)
job.molecule(geometry="water", charge=0)
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=3)

mol = job.run()
print(mol.get_results()["td_energies"])
```

### Constructor

```python
OpenQP(
    project="oqp_project",
    log=None,
    silent=0,
    usempi=True,
    config=None,
    **sections,
)
```

| Argument | Meaning |
| --- | --- |
| `project` | Project name used in logs and output files. It can be passed positionally, as in `OpenQP("h2o")`. |
| `log` | Log file path. Defaults to `<project>.log`. |
| `silent` | `0` prints parsed input and normal messages; `1` suppresses them. |
| `usempi` | Low-level runtime flag. New high-level scripts normally set this through `job.control(usempi=...)`. |
| `config` | Optional sectioned OpenQP input dictionary. |
| `**sections` | Optional section dictionaries, for example `input={...}` or `scf={...}`. |

The constructor is intentionally runtime-focused. Molecule and method setup use
OpenQP-native methods after construction.

The high-level API is organized around seven top-level calls:

| Call | Role |
| --- | --- |
| `job.molecule(...)` | Molecular identity: geometry, charge, multiplicity, and optional second geometry. |
| `job.d4(...)` | Enable DFT-D4 with functional defaults or a complete explicit rational-damping parameter set. |
| `job.theory.<model>(...)` | Quantum theory: method, functional, ordinary basis, reference, and response-state count. |
| `job.workflow.*(...)` | Calculation type: energy, gradient, Hessian, optimization, SOC, NACME, EKT, PCM, NMR, and related workflows. |
| `job.control(...)` | Hardware/runtime controls such as MPI and OpenMP threads. |
| `job.settings.*(...)` | Specialized detailed settings such as atom-wise basis and direct section overrides. |
| `job.run()` | Execute the calculation. |

### Molecule Setup

```python
job.molecule("H 0 0 0; H 0 0 0.74", charge=0, multiplicity=1)
job.molecule(geometry="water", charge=0, multiplicity=1)
job.molecule(geometry="benzene", source="pubchem", charge=0, multiplicity=1)
```

```python
job.molecule([
    ("H", 0.0, 0.0, 0.0),
    ("H", 0.0, 0.0, 0.74),
])
```

| Method | Returns | Use |
| --- | --- | --- |
| `molecule(system=None, system2=None, charge=None, multiplicity=None, unit="Angstrom", geometry=None, geometry2=None, source="auto", timeout=10, **kwargs)` | `OpenQP` | Writes molecular data into `[input]`, including explicit `system` text, optional `system2`, named `geometry`, charge, and the reference multiplicity. |
| `d4(enabled=True, s6=..., s8=..., s9=..., a1=..., a2=..., alp=...)` | `OpenQP` | Enables DFT-D4 and optionally writes the complete explicit rational-damping set into `[d4]`. Omit all six damping values to use functional defaults. |

Inline coordinates are Angstrom by default. Use `unit="Bohr"` for Bohr input
coordinates. `geometry=...` accepts built-in small molecules first and can fall
back to PubChem when `source="auto"` or `source="pubchem"` is used.
For compatibility with earlier scripts, `molecule(...)` can still accept
`basis=...`; new scripts should put basis in `theory(...)`.

The molecular charge set by `job.molecule(..., charge=...)` is passed to the
DFT-D4 charge model. For example:

```python
job = OpenQP("water_cation_d4")
job.molecule(geometry="water", charge=1, multiplicity=2)
job.theory.dft(functional="pbe", basis="6-31g*", reference="uhf")
job.d4(
    s6=1.0, s8=0.95948085, s9=1.0,
    a1=0.38574991, a2=4.80688534, alp=16.0,
)
```

### Theory, Workflow, and Control

```python
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=3)
job.workflow.gradient(state=3)
job.control(omp_threads=8, usempi=False)

dft_job = OpenQP("water_opt").molecule(geometry="water")
dft_job.theory.dft(functional="pbe0", basis="6-31g*")
dft_job.workflow.optimize(istate=0, coordsys="tric", trust=0.2)

mp2_job = OpenQP("water_mp2").molecule(geometry="water")
mp2_job.theory.mp2(basis="6-31g", reference="rhf", variant="scs-mp2")
mp2_job.workflow.gradient(state=0)
```

| Method | Returns | Use |
| --- | --- | --- |
| `theory.hf(basis=None, reference="rhf", multiplicity=None, **scf_keywords)` | `OpenQP` | Selects Hartree-Fock reference-SCF theory. |
| `theory.dft(functional, basis=None, reference="rhf", multiplicity=None, **scf_keywords)` | `OpenQP` | Selects Kohn-Sham DFT. The functional is part of the theory call. |
| `theory.mp2(reference="rhf", runtype=None, multiplicity=None, basis=None, variant=None, same_spin_scale=None, opposite_spin_scale=None, **scf_keywords)` | `OpenQP` | Selects MP2 with an HF reference. RHF supports analytic gradients and gradient-driven geometry runtypes; UHF/ROHF are energy-only. Use `variant` for named spin scaling, or `variant="custom"` with explicit same- and opposite-spin scales. |
| `theory.ccsd(reference="rhf", runtype=None, multiplicity=None, basis=None, nfzc=None, conv=None, maxit=None, ndiis=None, cholesky=None, cholesky_tol=None, cholesky_direct=None, triples=False, **scf_keywords)` | `OpenQP` | Selects energy-only CCSD with an HF reference. `nfzc`, `conv`, `maxit`, `ndiis`, `cholesky`, `cholesky_tol`, and `cholesky_direct` route to `[cc]`; every other keyword routes to `[scf]`. Passing a non-empty `functional`, or any `runtype` other than `energy`, raises. |
| `theory.ccsd_t(**kwargs)` | `OpenQP` | `theory.ccsd(...)` with `triples=True`, adding the perturbative triples correction. Takes the same arguments. |
| `fci(nroot=1, frozen_core=None, runtype=None, basis=None, reference="rhf", **keywords)` | `OpenQP` | Selects determinant-space FCI. |
| `casci(active_electrons=None, active_orbitals=None, frozen_core=None, nroot=1, runtype=None, basis=None, reference="rhf", **keywords)` | `OpenQP` | Selects fixed-orbital CASCI. |
| `casscf(active_electrons=None, active_orbitals=None, frozen_core=None, nroot=1, root=None, converger=None, hessian=None, max_macro_iterations=None, runtype=None, basis=None, reference="rhf", **keywords)` | `OpenQP` | Selects state-specific CASSCF. |
| `sa_casscf(active_electrons=None, active_orbitals=None, frozen_core=None, nstate=2, weights=None, target_roots=None, state=None, gradient_state=None, runtype=None, basis=None, reference="rhf", **keywords)` | `OpenQP` | Selects state-averaged CASSCF. For gradient run types, `gradient_state="averaged"` differentiates the weighted objective and an integer differentiates that physical CI root through the analytic Z-vector response. `state` remains a positional compatibility alias. |
| `caspt2(active_electrons=None, active_orbitals=None, frozen_core=None, nroot=None, variant=None, h0=None, ipea_shift=None, imaginary_shift=None, level_shift=None, runtype=None, basis=None, reference="rhf", **keywords)` | `OpenQP` | Selects CASPT2 or a multistate variant. |
| `nevpt2(active_electrons=None, active_orbitals=None, frozen_core=None, nroot=1, contraction=None, runtype=None, basis=None, reference="rhf", **keywords)` | `OpenQP` | Selects the Dyall-Hamiltonian PT2 path; `contraction="strong"` requests SC-NEVPT2. |
| `qdpt2(active_electrons=None, active_orbitals=None, frozen_core=None, nroot=None, variant=None, edshft=None, runtype=None, basis=None, reference="rhf", **keywords)` | `OpenQP` | Selects the QDPT family. |
| `theory.tdhf(basis=None, nstate=3, reference="rhf", multiplicity=1, **tdhf_keywords)` | `OpenQP` | Selects TDHF response theory. |
| `theory.tddft(functional, basis=None, nstate=3, reference="rhf", multiplicity=1, **tdhf_keywords)` | `OpenQP` | Selects TDDFT response theory. |
| `theory.sf_tddft(functional, basis=None, nstate=3, reference="rohf", multiplicity=3, **tdhf_keywords)` | `OpenQP` | Selects spin-flip TDDFT. `theory.sf(...)` is an alias. |
| `theory.mrsf(functional=None, basis=None, nstate=3, reference="rohf", **tdhf_keywords)` | `OpenQP` | Selects MRSF-TDDFT and supplies the usual triplet ROHF reference implicitly. `theory.mrsf_tddft(...)` is an alias. |
| `theory(method, functional=None, basis=None, nstate=3, reference=None, **keywords)` | `OpenQP` | Backward-compatible string dispatcher for existing scripts. In addition to HF/DFT, correlated, and response names, it accepts `"fci"`, `"casci"`, `"casscf"`, `"sa-casscf"`, `"caspt2"`, `"nevpt2"`, and `"qdpt2"` and dispatches them to the helpers above. `"cc"` is an alias for `"ccsd"`; `"ccsd-t"` and `"ccsdt"` for `"ccsd(t)"`. |
| `control(omp_threads=None, usempi=None, **kwargs)` | `OpenQP` | Sets hardware/runtime controls such as `[input] omp_threads` and the runtime-only MPI flag. |
| `workflow.gradient(state=None, **kwargs)` | `OpenQP` | Selects `runtype=grad` and stores the gradient state in `[properties] grad`. |
| `workflow.hessian(**kwargs)` | `OpenQP` | Selects `runtype=hess` and stores Hessian controls in `[hess]`. |
| `workflow.optimize(**kwargs)` | `OpenQP` | Selects `runtype=optimize`; native controls route to `[oqp]` by default, while explicit compatibility backends may route to `[geometric]` or legacy SciPy controls. |
| `workflow.meci(**kwargs)` | `OpenQP` | Selects `runtype=meci`. Two states default to `algorithm="auglag"`; use `states=[...]` with `algorithm="baeka"` for the BaekA two-or-more-state algorithm. Public controls include `gap_sigma` for `auglag`, and `sigma`, `alpha`, `delta_beta`, `beta_schedule`, and `gap`. The same workflow style is available for `mecp`, `mep`, `ts`, `irc`, and `neb`; `workflow.mecp(algorithm=...)` selects the MECP objective and accepts `gap_sigma`, and `workflow.tci(...)` remains a compatibility path. |
| `workflow.nacme(**kwargs)` | `OpenQP` | Selects `runtype=nacme` and requires MRSF-TDDFT. |
| `workflow.ekt(ip=False, ea=False, **kwargs)` | `OpenQP` | Selects `runtype=ekt`, requires MRSF-TDDFT, and requires IP, EA, or both. |
| `workflow.soc(soc_2e=1, scal_rel=2, **tdhf_keywords)` | `OpenQP` | Selects `runtype=soc` for an already configured MRSF-TDDFT theory, sets DKH2 scalar relativity by default, and rejects non-MRSF theories. |
| `workflow.pcm(**kwargs)` | `OpenQP` | Selects the current energy-only PCM/ddX path and requires HF/DFT reference-SCF, RHF/ROHF, `backend="ddx"`, and `mode="reference_scf"`. |
| `workflow.nmr(gauge="cgo", **kwargs)` | `OpenQP` | Requests NMR shielding and requires HF/DFT reference-SCF. CGO is RHF-only; use GIAO for open-shell references. |

The Python workflow API retains explicit `lib="geometric"` and `lib="scipy"`
selection for compatibility. This is intentionally broader than concise
`.oqp`, whose geometry drivers always use the native engine and reject backend
selectors. Native `oqp` handles frozen-distance minima; install
`openqp[geometric]` only before explicitly selecting geomeTRIC.

Plain energy calculations do not need a workflow call. Use `job.workflow.<name>(...)`
only when selecting a non-energy workflow or setting workflow-specific controls.
`job.control(...)`, the string form `job.theory("mrsf-tddft", ...)`, and the
older compact `job.hf(...)`, `job.dft(...)`, `job.mrsf(...)`, and `job.soc(...)`
helpers remain available for existing scripts.

MP2 clears any DFT functional and accepts named or custom spin-scaling
variants. It defaults to `runtype=energy`; an RHF calculation may instead use
`workflow.gradient`, `workflow.optimize`, `workflow.ts`, `workflow.mep`, or
`workflow.irc`:

```python
job.theory.mp2(basis="6-31g", reference="uhf", variant="scs-mp2", conv=1.0e-10)
job.theory("mp2", basis="6-31g", variant="custom",
           same_spin_scale=0.5, opposite_spin_scale=1.1)

grad_job = OpenQP("water_rmp2_grad").molecule(geometry="water")
grad_job.theory.mp2(basis="6-31g", reference="rhf", variant="sos-mp2")
grad_job.workflow.gradient(state=0)
```

The Python gradient helper uses `state=...` because users choose a molecular
state, even though the input-file keyword remains `[properties] grad`. For
HF/DFT and RHF MP2, `state=0` means the reference ground state. For ordinary TDHF/TDDFT,
`state=1` is the first excited state. For SF-TDDFT and MRSF-TDDFT, `state=1`
is the lowest spin-flip/MRSF target state, which can be the multiconfigurational
ground state. Existing scripts that use `grad=...` still work.

### Detailed Settings

```python
job.settings.scf(conv=1.0e-8)
job.settings.tdhf(target=2)

job.workflow.optimize(istate=0, coordsys="tric", trust=0.2)

job.settings.tdhf.nstate = 5

job.settings.basis(["LANL2DZ", "6-31g*"])
job.settings.basis(c1="cc-pvdz", h1="6-31g*")

job.set(**{"input.method": "tdhf", "tdhf.type": "mrsf"})
job.update({"scf": {"type": "rohf", "multiplicity": 3}})
```

| Method | Returns | Use |
| --- | --- | --- |
| `section(name, **kwargs)` | `OpenQP` | Updates one OpenQP input section. |
| `settings.basis(basis=None, **tags)` | `OpenQP` | Sets atom-wise basis assignment. Ordinary global basis selection belongs in `theory(..., basis=...)`. |
| `settings.<section>(**kwargs)` | `OpenQP` | Updates one detailed OpenQP input section, for example `settings.scf(conv=...)`. |
| `set(**kwargs)` | `OpenQP` | Updates dotted OpenQP keywords or section dictionaries. |
| `update(config=None, **kwargs)` | `OpenQP` | Merges a sectioned dictionary plus optional section overrides. |
| `to_input_dict()` | `dict` | Returns the sectioned dictionary that will be passed to `Runner`. |
| `run(run_type=None)` | `Molecule` | Builds `Runner`, executes the calculation, stores `job.runner` and `job.mol`, and returns the `Molecule`. |

For optimization workflows, `job.workflow.optimize(...)` routes ordinary optimization
keywords to `[optimize]`, while backend options such as `coordsys`, `trust`, and
`constraints_file` are sent to the selected backend section. The lower-level
`job.optimize(...)` section helper remains available for existing scripts.
The native backend is the default; `constraints_file` is relevant to the
optional legacy geomeTRIC path.

`oqp.openqp.OQP` is an alias for `OpenQP`.

### PySCF Conversion

`OpenQP.from_pyscf(mol, **kwargs)` is an explicit compatibility bridge. It
reads `atom`, `basis`, `charge`, `spin`, and `unit` attributes from a PySCF-like
object, translates them into OpenQP sections, and returns an `OpenQP` job.

```python
job = OpenQP.from_pyscf(pyscf_mol, project="mixed_workflow")
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=5)
mol = job.run()
```

After conversion, use normal OpenQP workflow and section calls.

## Runner

`oqp.pyoqp.Runner` loads the same configuration used by input files, runs the
selected workflow, and keeps the resulting `Molecule` object available as
`runner.mol`.

```python
from oqp.pyoqp import Runner

runner = Runner(
    project="water_mrsf",
    input_file="water_mrsf.inp",
    log="water_mrsf.log",
    silent=0,
    usempi=True,
)
runner.run()

summary = runner.results()
print(summary["energy"])
```

### Signature

```python
Runner(
    project=None,
    input_file=None,
    log=None,
    input_dict=None,
    silent=0,
    usempi=True,
)
```

| Argument | Type | Meaning |
| --- | --- | --- |
| `project` | `str` or `None` | Project name used in logs and output files. |
| `input_file` | `str` or `None` | Path to an OpenQP input file. |
| `log` | `str` or `None` | Log file path. |
| `input_dict` | `dict` or `None` | Sectioned input dictionary. If supplied, it is used instead of `input_file`. |
| `silent` | `int` | `0` prints parsed input and normal messages; `1` suppresses them. |
| `usempi` | `bool` | Enables MPI-aware behavior when the runtime supports it. |

`Runner` validates the parsed input before the calculation starts. A validation
failure prints an input-check report and exits before expensive kernels run.

### File-Based Runs

File-based runs are closest to command-line OpenQP:

```python
from pathlib import Path
from oqp.pyoqp import Runner

input_file = Path("examples/HF/H2O_RHF-HF_ENERGY.inp").resolve()
project = input_file.stem

runner = Runner(
    project=project,
    input_file=str(input_file),
    log=f"{project}.log",
)
runner.run()
```

Use this pattern when you want reproducible inputs, logs, JSON restart files,
and output files on disk.

### In-Memory Runs

For notebooks, tests, services, or agents, pass a sectioned dictionary through
`input_dict`.

```python
from oqp.pyoqp import Runner

config = {
    "input": {
        "system": "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74",
        "basis": "6-31g*",
        "method": "hf",
        "runtype": "energy",
    },
    "scf": {
        "type": "rhf",
        "multiplicity": "1",
    },
}

runner = Runner(
    project="h2",
    input_dict=config,
    log="h2.log",
    silent=1,
)
runner.run()
print(runner.results()["energy"])
```

The dictionary is organized as `{section: {keyword: value}}`, matching the input
file sections. Values can be strings; the OpenQP parser converts them using the
same schema used by file inputs.

### MRSF-TDDFT Example

MRSF-TDDFT normally starts from an open-shell reference. This compact example
uses an ROHF triplet reference and asks for MRSF-TDDFT states.

```python
from oqp.pyoqp import Runner

config = {
    "input": {
        "system": "\nO 0.000000 0.000000 -0.041062\nH -0.533194 0.533194 -0.614469\nH 0.533194 -0.533194 -0.614469",
        "basis": "6-31g*",
        "method": "tdhf",
        "runtype": "energy",
    },
    "scf": {
        "type": "rohf",
        "multiplicity": "3",
    },
    "tdhf": {
        "type": "mrsf",
        "nstate": "3",
    },
}

runner = Runner(project="water_mrsf", input_dict=config, log="water_mrsf.log")
runner.run()

mol = runner.mol
print(mol.get_results()["td_energies"])
```

### Runtime Methods

| Method | Returns | Use |
| --- | --- | --- |
| `run(test_mod=False)` | `None` | Executes the calculation selected by `[input] runtype`. |
| `results()` | `dict` | Returns a Python summary with atoms, coordinates, energies, gradients, NAC, SOC, and raw data tags. |
| `reload()` | `None` | Reloads guess data when the selected guess mode asks for JSON data. |
| `back_door(data)` | `None` | Supplies previous-state data for advanced internal workflows. |
| `test()` | `(message, diff)` | Compares against reference data in OpenQP test mode. |

## Legacy Dotted-Keyword Wrapper

`oqp.openqp.OPENQP` accepts dotted keyword names such as `input.system` and
`scf.type`, normalizes inline geometries, and then constructs a `Runner`.

```python
from oqp.openqp import OPENQP

op = OPENQP({
    "input.system": "H 0 0 0; H 0 0 0.74",
    "input.basis": "6-31g*",
    "input.method": "hf",
    "input.runtype": "energy",
    "scf.type": "rhf",
})
mol = op.run()
```

This wrapper is retained for existing scripts. New scripts should prefer
`OpenQP` for OpenQP-native scripting or `Runner` for direct input-file and
section-dictionary execution.
