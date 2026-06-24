"""
SimulationManager — Unified Python interface for multi-scale simulations in OpenQP.

Supports four modes:
    qp        : Pure QM single-point (OpenQP)
    qmmm      : QM/MM (OpenQP + OpenMM via QMMM_MD)
    namd      : Non-adiabatic MD (PyRAI2MD + OpenQP)
    namd_qmmm : NAMD with QM/MM potential (PyRAI2MD + OpenQP + OpenMM)

Minimal usage
-------------
    manager = SimulationManager.from_template(
        "template_namd.json",
        folder="Nac_run",
        title="test1-1",
        xyz="test1-1.xyz",
    )
    manager.show_config()
    manager.run()
"""

import json
import os
from copy import deepcopy
from pathlib import Path


# ====================================================================
#  Element symbol → atomic number lookup
# ====================================================================
ELEMENT_TO_Z = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
    "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
}


# ====================================================================
#  Utility
# ====================================================================
def _resolve_content(value):
    """If *value* is a path to an existing file, read its text.
    Otherwise treat as a raw string."""
    if isinstance(value, (str, Path)):
        s = str(value)
        if "\n" not in s and len(s) < 260:
            p = Path(s)
            try:
                if p.is_file():
                    return p.read_text(encoding="utf-8"), True
            except OSError:
                pass
    return str(value), False


# ====================================================================
#  Dependency checking
# ====================================================================
def _check_package(name):
    """Check if a Python package is importable. Returns (available, version_or_error)."""
    try:
        mod = __import__(name)
        version = getattr(mod, "__version__", "installed")
        return True, version
    except ImportError:
        return False, "not installed"


def _check_all_dependencies():
    """Check all packages and return a dict of results."""
    packages = {
        "oqp":        "OpenQP      (QM engine)",
        "openmm":     "OpenMM      (MM engine)",
        "PyRAI2MD":   "PyRAI2MD    (NAMD driver)",
    }
    results = {}
    for pkg, label in packages.items():
        available, info = _check_package(pkg)
        results[pkg] = {"label": label, "available": available, "info": info}
    return results


# Which packages each mode requires
_MODE_DEPENDENCIES = {
    "qp":        ["oqp"],
    "qmmm":      ["oqp", "openmm"],
    "namd":      ["oqp", "PyRAI2MD"],
    "namd_qmmm": ["oqp", "openmm", "PyRAI2MD"],
}


def _xyz_to_system_block(xyz):
    """Convert XYZ string/file → OpenQP system= lines with atomic numbers."""
    text, _ = _resolve_content(xyz)
    system_lines = []
    for line in text.strip().splitlines():
        parts = line.split()
        if len(parts) >= 4:
            try:
                float(parts[1]); float(parts[2]); float(parts[3])
            except ValueError:
                continue
            z = ELEMENT_TO_Z.get(parts[0], parts[0])
            system_lines.append(f"   {z}   {parts[1]}   {parts[2]}   {parts[3]}")
    return "\n".join(system_lines)


# ====================================================================
#  Input builders
# ====================================================================
def _sections_to_inp(sections, xyz=None):
    """Convert a dict of sections → OpenQP .inp / .openqp file content."""
    lines = []
    for section_name, params in sections.items():
        lines.append(f"[{section_name}]")
        if section_name == "input" and xyz is not None:
            system_block = _xyz_to_system_block(xyz)
            lines.append(f"system=\n{system_block}")
        for key, value in params.items():
            if isinstance(value, bool):
                value = str(value).lower()
            lines.append(f"{key}={value}")
        lines.append("")
    return "\n".join(lines)


def _namd_dict_to_control(namd_params, title):
    """
    Convert the namd dict → PyRAI2MD control file content.

    The 'openqp' key in &openqp section is auto-filled from
    the $OPENQP_ROOT environment variable if not set explicitly.
    """
    section_order = ["control", "molecule", "openqp", "md"]
    section_labels = {
        "control": "&CONTROL",
        "molecule": "&MOLECULE",
        "openqp": "&openqp",
        "md": "&MD",
    }

    # resolve openqp path from environment if not set
    openqp_root = os.environ.get("OPENQP_ROOT", "")

    lines = []
    for sec in section_order:
        if sec not in namd_params:
            continue

        lines.append(section_labels[sec])

        # inject title as first entry in &CONTROL
        if sec == "control":
            lines.append(f"title    {title}")

        params = namd_params[sec]
        for key, value in params.items():
            if isinstance(value, bool):
                value = str(value).lower()

            # auto-fill openqp path from $OPENQP_ROOT
            if sec == "openqp" and key == "openqp" and (value == "" or value is None):
                if openqp_root:
                    value = openqp_root
                else:
                    raise EnvironmentError(
                        "OpenQP path not set. Either:\n"
                        "  1. Set the $OPENQP_ROOT environment variable, or\n"
                        "  2. Pass it explicitly:\n"
                        '     manager.set_namd_config("openqp", "openqp", "/path/to/openqp")'
                    )

            # handle empty values (like 'restart' with no value)
            if value == "" or value is None:
                lines.append(f"{key}")
            else:
                lines.append(f"{key} {value}")

        lines.append("")

    return "\n".join(lines)


# ====================================================================
#  Built-in templates
# ====================================================================
BUILTIN_TEMPLATES = {
    "qp_default": {
        "description": "Pure QM single-point calculation",
        "mode": "qp",
        "sections": {
            "input": {
                "functional": "bhhlyp", "basis": "6-31g*",
                "method": "hf", "runtype": "energy", "charge": 0,
            },
            "scf": {"type": "rhf", "maxit": 100, "multiplicity": 1},
            "guess": {"type": "huckel"},
            "dftgrid": {"rad_type": "becke"},
        },
    },
    "qmmm_default": {
        "description": "QM/MM with OpenMM (MRSF-TDDFT)",
        "mode": "qmmm",
        "sections": {
            "input": {
                "functional": "bhhlyp", "basis": "6-31g*",
                "method": "tdhf", "qmmm_flag": True,
            },
            "scf": {"type": "rohf", "maxit": 100, "multiplicity": 1},
            "guess": {"type": "huckel"},
            "dftgrid": {"rad_type": "becke"},
            "tdhf": {"type": "mrsf", "nstate": 6},
            "qmmm": {
                "pdb_file": "", "forcefield_files": "amber14-all.xml",
                "qm_atoms": "", "cutoff": "PME", "embedding": "electrostatic",
                "n_steps": 1000, "timestep": 1.0, "temperature": 300.0,
            },
        },
    },
    "namd_default": {
        "description": "Non-adiabatic MD with PyRAI2MD (MRSF-TDDFT)",
        "mode": "namd",
        "sections": {
            "input": {
                "functional": "bhhlyp", "charge": 0,
                "method": "tdhf", "basis": "6-31G*",
            },
            "guess": {"type": "huckel"},
            "scf": {
                "type": "rohf", "maxit": 200, "multiplicity": 3,
                "converger_type": "soscf", "conv": "1e-8",
            },
            "tdhf": {
                "type": "mrsf", "nstate": 6, "maxit": 100,
                "maxit_zv": 200, "nvdav": 100, "zvconv": "1.0e-10",
            },
            "dftgrid": {"rad_npts": 96, "ang_npts": 302, "pruned": ""},
            "properties": {
                "export": True, "nac": "nacme", "back_door": True, "grad": 5,
            },
            "nac": {},
        },
        "namd": {
            "control": {
                "qc_ncpu": 16, "gl_seed": 1, "jobtype": "md",
                "qm": "openqp", "abinit": "openqp",
            },
            "molecule": {
                "ci": 5, "spin": 0,
                "coupling": "1 2, 1 3, 1 4, 1 5, 2 3, 2 4, 2 5, 3 4, 3 5, 4 5",
            },
            "openqp": {"openqp": "", "threads": 40, "use_hpc": -1, "align_mo": 1},
            "md": {
                "initcond": 0, "temp": 300, "randvelo": 0,
                "step": 80, "direct": 80, "checkpoint": 1, "buffer": 0,
                "size": 20.67, "root": 5, "activestate": 1,
                "sfhp": "fssh", "nactype": "dcm", "phasecheck": 0,
                "adjust": 0, "reflect": 1, "deco": "OFF",
                "substep": 20, "verbose": 1, "thermo": 0, "restart": "",
            },
        },
    },
    "namd_qmmm_default": {
        "description": "Non-adiabatic MD with QM/MM (MRSF-TDDFT + OpenMM)",
        "mode": "namd_qmmm",
        "sections": {
            "input": {
                "functional": "bhhlyp", "charge": 0,
                "method": "tdhf", "basis": "6-31G*", "qmmm_flag": True,
            },
            "guess": {"type": "huckel"},
            "scf": {
                "type": "rohf", "maxit": 200, "multiplicity": 3,
                "converger_type": "soscf", "conv": "1e-8",
            },
            "tdhf": {
                "type": "mrsf", "nstate": 6, "maxit": 100,
                "maxit_zv": 200, "nvdav": 100, "zvconv": "1.0e-10",
            },
            "dftgrid": {"rad_npts": 96, "ang_npts": 302, "pruned": ""},
            "properties": {
                "export": True, "nac": "nacme", "back_door": True, "grad": 5,
            },
            "nac": {},
            "qmmm": {
                "pdb_file": "", "forcefield_files": "amber14-all.xml",
                "qm_atoms": "", "cutoff": "PME", "embedding": "electrostatic",
                "n_steps": 1000, "timestep": 1.0, "temperature": 300.0,
            },
        },
        "namd": {
            "control": {
                "qc_ncpu": 16, "gl_seed": 1, "jobtype": "md",
                "qm": "openqp", "abinit": "openqp",
            },
            "molecule": {
                "ci": 5, "spin": 0,
                "coupling": "1 2, 1 3, 1 4, 1 5, 2 3, 2 4, 2 5, 3 4, 3 5, 4 5",
            },
            "openqp": {"openqp": "", "threads": 40, "use_hpc": -1, "align_mo": 1},
            "md": {
                "initcond": 0, "temp": 300, "randvelo": 0,
                "step": 80, "direct": 80, "checkpoint": 1, "buffer": 0,
                "size": 20.67, "root": 5, "activestate": 1,
                "sfhp": "fssh", "nactype": "dcm", "phasecheck": 0,
                "adjust": 0, "reflect": 1, "deco": "OFF",
                "substep": 20, "verbose": 1, "thermo": 0, "restart": "",
            },
        },
    },
}


# ====================================================================
#  OpenQPEngine
# ====================================================================
class OpenQPEngine:
    Runner = None

    def __init__(self, config=None):
        self.config = config or {}
        self.project = self.config.get("project", "default")
        self.workdir = self.config.get("workdir", ".")
        self.input_dict = self.config.get("input_dict", {})
        self.back_door_data = self.config.get("back_door_data", {})
        self.results = None

    def setup_input(self, folder, openqp_content, title):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        self.workdir = str(folder)
        self.project = title
        text, _ = _resolve_content(openqp_content)
        (folder / f"{title}.inp").write_text(text.strip() + "\n", encoding="utf-8")
        return folder

    def setup_input_with_xyz(self, folder, openqp_content, xyz_content, title):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        self.workdir = str(folder)
        self.project = title
        oqp_text, _ = _resolve_content(openqp_content)
        (folder / f"{title}.inp").write_text(oqp_text.strip() + "\n", encoding="utf-8")
        xyz_text, _ = _resolve_content(xyz_content)
        (folder / f"{title}.xyz").write_text(xyz_text.strip() + "\n", encoding="utf-8")
        return folder

    def prepare(self):
        if OpenQPEngine.Runner is None:
            from oqp.pyoqp import Runner
            OpenQPEngine.Runner = Runner
        Path(self.workdir).mkdir(parents=True, exist_ok=True)
        print(f"OpenQP engine prepared  (project={self.project}, workdir={self.workdir})")

    def run(self):
        if OpenQPEngine.Runner is None:
            self.prepare()
        pyoqp = OpenQPEngine.Runner(
            project=self.project,
            input_file=f"{self.workdir}/{self.project}.inp",
            input_dict=self.input_dict,
            log=f"{self.workdir}/{self.project}.log",
            silent=1, usempi=False,
        )
        original_dir = os.getcwd()
        try:
            os.chdir(self.workdir)
            pyoqp.back_door(self.back_door_data)
            pyoqp.run()
            self.results = pyoqp.results()
        finally:
            os.chdir(original_dir)
        print(f"OpenQP finished  (project={self.project})")
        return self.results

    def run_qmmm(self):
        from oqp.library.qmmm_md import QMMM_MD
        input_file = f"{self.project}.inp"
        original_dir = os.getcwd()
        try:
            os.chdir(self.workdir)
            print(f"Running QM/MM in: {self.workdir} (cfg={input_file})")
            md = QMMM_MD(oqp_cfg=input_file)
            md.run()
        finally:
            os.chdir(original_dir)
        print(f"QM/MM finished  (project={self.project})")


# ====================================================================
#  PyRAI2MDEngine
# ====================================================================
class PyRAI2MDEngine:
    def __init__(self, mode, config=None):
        self.mode = mode
        self.config = config or {}

    def prepare(self):
        print("Preparing PyRAI2MD engine")

    def _set_runtype(self, openqp_content):
        runtype = "prop" if self.mode == "namd" else "qmmm"
        lines, found = [], False
        for line in openqp_content.splitlines():
            if line.strip().startswith("runtype"):
                lines.append(f"runtype={runtype}")
                found = True
            else:
                lines.append(line)
        if not found:
            new_lines = []
            for line in lines:
                new_lines.append(line)
                if line.strip().lower() == "[input]":
                    new_lines.append(f"runtype={runtype}")
            lines = new_lines
        return "\n".join(lines)

    def setup_files(self, folder, input_content, openqp_content, xyz_content, title):
        """
        Write the three NAMD files into folder:
            title           — PyRAI2MD control file
            title.openqp    — OpenQP electronic structure input
            title.xyz       — molecular coordinates
        """
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)

        # PyRAI2MD control file (no extension)
        inp_text, _ = _resolve_content(input_content)
        (folder / title).write_text(inp_text.strip() + "\n", encoding="utf-8")

        # OpenQP input (.openqp) — inject runtype
        oqp_text, _ = _resolve_content(openqp_content)
        oqp_text = self._set_runtype(oqp_text)
        (folder / f"{title}.openqp").write_text(oqp_text.strip() + "\n", encoding="utf-8")

        # XYZ coordinates
        xyz_text, _ = _resolve_content(xyz_content)
        (folder / f"{title}.xyz").write_text(xyz_text.strip() + "\n", encoding="utf-8")

        return folder

    def run(self, folder, title):
        from PyRAI2MD.pyrai2md import PYRAI2MD
        folder = Path(folder).resolve()
        original_dir = os.getcwd()
        try:
            os.chdir(folder)
            print(f"Running PyRAI2MD in: {folder}")
            pmd = PYRAI2MD(title)
            pmd.run()
        finally:
            os.chdir(original_dir)


# ====================================================================
#  SimulationManager
# ====================================================================
class SimulationManager:
    """
    Unified interface for QM, QM/MM, NAMD, and NAMD-QM/MM simulations.

    For NAMD modes, three files are generated in the working folder:
        title           — PyRAI2MD control (&CONTROL, &MOLECULE, &openqp, &MD)
        title.openqp    — OpenQP input ([input], [scf], [tdhf], [properties], ...)
        title.xyz       — molecular coordinates

    Parameters
    ----------
    mode : str
        One of 'qp', 'qmmm', 'namd', 'namd_qmmm'.
    folder : str
        Working directory.
    title : str
        Project name (used for all filenames).
    xyz : str or path
        XYZ coordinates.
    pdb : str or path, optional
        PDB file for QM/MM modes.
    qm_atoms : str, optional
        QM atom selection for QM/MM (e.g. '0-2').
    sections : dict, optional
        OpenQP section overrides.
    namd : dict, optional
        PyRAI2MD section overrides (control, molecule, openqp, md).
    template : str, optional
        JSON template path or built-in name.
    openqp_content, xyz_content, input_content : str, optional
        Pre-built inputs (bypass auto-generation).
    configs : dict, optional
        Engine-level config overrides.
    """

    VALID_MODES = {"qp", "qmmm", "namd", "namd_qmmm"}

    def __init__(
        self,
        mode=None,
        folder=".",
        title="test1-1",
        xyz=None,
        pdb=None,
        qm_atoms=None,
        sections=None,
        namd=None,
        template=None,
        openqp_content=None,
        xyz_content=None,
        input_content=None,
        configs=None,
    ):
        # ── Load template ──
        tpl = {}
        if template is not None:
            tpl = self._load_template(template)

        self.mode = mode or tpl.get("mode")
        if self.mode not in self.VALID_MODES:
            raise ValueError(
                f"Unsupported mode '{self.mode}'. Valid: {sorted(self.VALID_MODES)}"
            )

        self.folder = folder
        self.title = title
        self.configs = configs or {}
        self._template_source = template

        self._xyz = xyz
        self._pdb = pdb
        self._qm_atoms = qm_atoms

        # ── Build OpenQP sections ──
        self._sections = deepcopy(tpl.get("sections", {}))
        if sections:
            for sec_name, sec_params in sections.items():
                if sec_name in self._sections:
                    self._sections[sec_name].update(sec_params)
                else:
                    self._sections[sec_name] = deepcopy(sec_params)
        if not self._sections:
            self._sections = self._default_sections()

        # inject molecule-specific values
        if self._pdb is not None and "qmmm" in self._sections:
            pdb_path = Path(self._pdb)
            self._sections["qmmm"]["pdb_file"] = (
                pdb_path.name if pdb_path.is_file() else f"{self.title}.pdb"
            )
        if self._qm_atoms is not None and "qmmm" in self._sections:
            self._sections["qmmm"]["qm_atoms"] = self._qm_atoms

        # ── Build NAMD params: defaults → template → user overrides ──
        if self.mode in {"namd", "namd_qmmm"}:
            self._namd_params = self._default_namd()
            # merge template on top of defaults
            for sec_name, sec_params in tpl.get("namd", {}).items():
                if sec_name in self._namd_params and isinstance(sec_params, dict):
                    self._namd_params[sec_name].update(deepcopy(sec_params))
                else:
                    self._namd_params[sec_name] = deepcopy(sec_params)
            # merge user overrides on top
            if namd:
                for sec_name, sec_params in namd.items():
                    if sec_name in self._namd_params and isinstance(sec_params, dict):
                        self._namd_params[sec_name].update(sec_params)
                    else:
                        self._namd_params[sec_name] = (
                            deepcopy(sec_params) if isinstance(sec_params, dict) else sec_params
                        )
        else:
            self._namd_params = {}

        # ── Check dependencies ──
        self._check_dependencies()

        # ── Generate inputs ──
        self._generated = {}
        self._generate_inputs(openqp_content, xyz_content, input_content)

        # ── Engines ──
        self.openqp = None
        self.pyrai2md = None
        self._setup_engines()
        self._setup_files()

    def _default_sections(self):
        """Default OpenQP sections per mode."""

        if self.mode == "qp":
            return {
                "input": {
                    "functional": "bhhlyp", "basis": "6-31g*",
                    "method": "hf", "runtype": "energy", "charge": 0,
                },
                "scf": {"type": "rhf", "maxit": 100, "multiplicity": 1},
                "guess": {"type": "huckel"},
                "dftgrid": {"rad_type": "becke"},
            }

        if self.mode == "qmmm":
            return {
                "input": {
                    "functional": "bhhlyp", "basis": "6-31g*",
                    "method": "tdhf", "qmmm_flag": True,
                },
                "scf": {
                    "type": "rhf", "maxit": 200,
                    "converger_type": "diis", "conv": "1e-8",
                },
                "guess": {"type": "huckel"},
                "dftgrid": {"rad_npts": 96, "ang_npts": 302, "pruned": ""},
                "qmmm": {
                    "pdb_file": "", "forcefield_files": "tip3p.xml",
                    "qm_atoms": "", "cutoff": "PME", "embedding": "electrostatic",
                    "n_steps": 1000, "timestep": 0.1, "temperature": 300.0,
                },
            }

        if self.mode == "namd":
            return {
                "input": {
                    "functional": "bhhlyp", "charge": 0,
                    "method": "tdhf", "basis": "6-31G*",
                },
                "guess": {"type": "huckel"},
                "scf": {
                    "type": "rohf", "maxit": 200, "multiplicity": 3,
                    "converger_type": "soscf", "conv": "1e-8",
                },
                "tdhf": {
                    "type": "mrsf", "nstate": 6, "maxit": 100,
                    "maxit_zv": 200, "nvdav": 100, "zvconv": "1.0e-10",
                },
                "dftgrid": {"rad_npts": 96, "ang_npts": 302, "pruned": ""},
                "properties": {
                    "export": True, "nac": "nacme", "back_door": True, "grad": 5,
                },
                "nac": {},
            }

        if self.mode == "namd_qmmm":
            return {
                "input": {
                    "functional": "bhhlyp", "charge": 0,
                    "method": "tdhf", "basis": "6-31G*", "qmmm_flag": True,
                },
                "guess": {"type": "huckel"},
                "scf": {
                    "type": "rohf", "maxit": 200, "multiplicity": 3,
                    "converger_type": "soscf", "conv": "1e-8",
                },
                "tdhf": {
                    "type": "mrsf", "nstate": 6, "maxit": 100,
                    "maxit_zv": 200, "nvdav": 100, "zvconv": "1.0e-10",
                },
                "dftgrid": {"rad_npts": 96, "ang_npts": 302, "pruned": ""},
                "properties": {
                    "export": True, "nac": "nacme", "back_door": True, "grad": 5,
                },
                "nac": {},
                "qmmm": {
                    "pdb_file": "", "forcefield_files": "amber14-all.xml",
                    "qm_atoms": "", "cutoff": "PME", "embedding": "electrostatic",
                    "n_steps": 1000, "timestep": 1.0, "temperature": 300.0,
                },
            }

        return {}

    def _default_namd(self):
        return {
            "control": {
                "qc_ncpu": 16, "gl_seed": 1, "jobtype": "md",
                "qm": "openqp", "abinit": "openqp",
            },
            "molecule": {"ci": 5, "spin": 0, "coupling": "1 2, 1 3, 1 4, 1 5, 2 3, 2 4, 2 5, 3 4, 3 5, 4 5"},
            "openqp": {"openqp": "", "threads": 40, "use_hpc": -1, "align_mo": 1},
            "md": {
                "initcond": 0, "temp": 300, "randvelo": 0,
                "step": 80, "direct": 80, "checkpoint": 1, "buffer": 0,
                "size": 20.67, "root": 5, "activestate": 1,
                "sfhp": "fssh", "nactype": "dcm", "phasecheck": 0,
                "adjust": 0, "reflect": 1, "deco": "OFF",
                "substep": 20, "verbose": 1, "thermo": 0, "restart": "",
            },
        }

    # ================================================================
    #  Dependency checking
    # ================================================================
    def _check_dependencies(self):
        """Check that required packages for the current mode are installed."""
        required = _MODE_DEPENDENCIES.get(self.mode, [])
        missing = []

        for pkg in required:
            available, info = _check_package(pkg)
            if not available:
                missing.append(pkg)

        if missing:
            pkg_names = {
                "oqp": "OpenQP",
                "openmm": "OpenMM",
                "PyRAI2MD": "PyRAI2MD",
            }
            names = [f"{pkg_names.get(p, p)} ({p})" for p in missing]
            raise ImportError(
                f"Mode '{self.mode}' requires the following packages that are not installed:\n"
                + "\n".join(f"  - {n}" for n in names)
                + "\n\nUse SimulationManager.check_dependencies() to see all package status."
            )

    @staticmethod
    def check_dependencies():
        """
        Print the installation status of all required packages.

        Example
        -------
            SimulationManager.check_dependencies()
            # OpenQP      (QM engine)    ✓ installed (0.3.1)
            # OpenMM      (MM engine)    ✓ installed (8.1.1)
            # PyRAI2MD    (NAMD driver)  ✗ not installed
        """
        results = _check_all_dependencies()
        print("Package dependencies:")
        print("-" * 60)
        for pkg, info in results.items():
            status = f"\u2713 {info['info']}" if info["available"] else "\u2717 not installed"
            print(f"  {info['label']}    {status}")

        print()
        print("Mode requirements:")
        print("-" * 60)
        for mode, pkgs in _MODE_DEPENDENCIES.items():
            all_ok = all(results[p]["available"] for p in pkgs)
            mark = "\u2713" if all_ok else "\u2717"
            pkg_list = ", ".join(pkgs)
            print(f"  {mark} {mode:12s}  needs: {pkg_list}")

    # ================================================================
    #  Template system
    # ================================================================
    @staticmethod
    def _load_template(source):
        if isinstance(source, str) and source in BUILTIN_TEMPLATES:
            return deepcopy(BUILTIN_TEMPLATES[source])
        path = Path(source)
        if not path.is_file():
            raise FileNotFoundError(
                f"Template not found: '{source}'. Built-in: {list(BUILTIN_TEMPLATES.keys())}"
            )
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)

    def save_template(self, filepath, description=None):
        """Save current settings as a reusable JSON template."""
        tpl = {
            "description": description or f"SimulationManager template ({self.mode})",
            "mode": self.mode,
            "sections": deepcopy(self._sections),
        }
        if self._namd_params:
            tpl["namd"] = deepcopy(self._namd_params)
        # clear molecule-specific values
        if "qmmm" in tpl["sections"]:
            tpl["sections"]["qmmm"]["pdb_file"] = ""
            tpl["sections"]["qmmm"]["qm_atoms"] = ""
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(tpl, f, indent=2, ensure_ascii=False)
        print(f"Template saved: {filepath}")

    @classmethod
    def from_template(cls, template, **overrides):
        """Create a SimulationManager from a template with per-job overrides."""
        return cls(template=template, **overrides)

    @staticmethod
    def list_templates(directory=None):
        """List built-in and user templates."""
        print("Built-in templates:")
        print("-" * 60)
        for name, tpl in BUILTIN_TEMPLATES.items():
            print(f"  {name:24s} {tpl.get('description', '')}")
        search_dir = Path(directory) if directory else Path(".")
        user_templates = []
        for f in sorted(search_dir.glob("*.json")):
            try:
                with open(f) as fh:
                    data = json.load(fh)
                if "mode" in data and "sections" in data:
                    user_templates.append((f.name, data.get("description", "")))
            except (json.JSONDecodeError, KeyError):
                continue
        if user_templates:
            print(f"\nUser templates in {search_dir}:")
            print("-" * 60)
            for name, desc in user_templates:
                print(f"  {name:24s} {desc}")

    @staticmethod
    def create_template(filepath, mode, description=None, sections=None, namd=None):
        """Create a template JSON file from scratch."""
        tpl = {
            "description": description or f"Custom template ({mode})",
            "mode": mode,
            "sections": sections or {},
        }
        if namd:
            tpl["namd"] = namd
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(tpl, f, indent=2, ensure_ascii=False)
        print(f"Template created: {filepath}")

    # ================================================================
    #  Input generation
    # ================================================================
    def _generate_inputs(self, openqp_content, xyz_content, input_content):
        # ── OpenQP .inp / .openqp ──
        if openqp_content is not None:
            text, _ = _resolve_content(openqp_content)
            self._generated["openqp"] = text
        else:
            xyz_for_system = self._xyz if self.mode == "qp" else None
            self._generated["openqp"] = _sections_to_inp(self._sections, xyz=xyz_for_system)

        # ── XYZ ──
        if xyz_content is not None:
            text, _ = _resolve_content(xyz_content)
            self._generated["xyz"] = text
        elif self._xyz is not None:
            text, _ = _resolve_content(self._xyz)
            self._generated["xyz"] = text

        # ── PDB ──
        if self._pdb is not None:
            text, _ = _resolve_content(self._pdb)
            self._generated["pdb"] = text

        # ── PyRAI2MD control file ──
        if input_content is not None:
            text, _ = _resolve_content(input_content)
            self._generated["namd"] = text
        elif self.mode in {"namd", "namd_qmmm"} and self._namd_params:
            self._generated["namd"] = _namd_dict_to_control(self._namd_params, self.title)

    # ================================================================
    #  Config inspection & modification
    # ================================================================
    def show_config(self, name=None):
        """
        Inspect generated inputs.

        name=None      → print all to console
        name='openqp'  → return OpenQP input string
        name='namd'    → return PyRAI2MD control string
        name='xyz'     → return XYZ string
        name='pdb'     → return PDB string
        name='sections'→ return OpenQP sections dict
        name='namd_sections' → return NAMD sections dict
        """
        if name == "sections":
            return deepcopy(self._sections)
        if name == "namd_sections":
            return deepcopy(self._namd_params)
        if name is not None:
            return self._generated.get(name)

        sep = "=" * 60
        for key, content in self._generated.items():
            print(f"\n{sep}")
            print(f"  {key.upper()} INPUT")
            print(sep)
            print(content)
        print(f"\n{sep}")
        print(f"  MODE: {self.mode}  |  FOLDER: {self.folder}  |  TITLE: {self.title}")
        if self._template_source:
            print(f"  TEMPLATE: {self._template_source}")
        print(sep)

    def set_config(self, section, key, value):
        """
        Modify OpenQP sections and regenerate .inp/.openqp.

        Example:
            manager.set_config("input", "basis", "cc-pvdz")
            manager.set_config("tdhf", "nstate", 8)
            manager.set_config("qmmm", "temperature", 350.0)
        """
        if section not in self._sections:
            self._sections[section] = {}
        self._sections[section][key] = value

        xyz_for_system = self._xyz if self.mode == "qp" else None
        self._generated["openqp"] = _sections_to_inp(self._sections, xyz=xyz_for_system)
        print(f"Config updated: [{section}] {key} = {value}")

    def set_namd_config(self, section, key, value):
        """
        Modify PyRAI2MD control sections and regenerate control file.

        Parameters
        ----------
        section : str
            One of 'control', 'molecule', 'openqp', 'md'.
        key : str
            Parameter name (e.g. 'step', 'temp', 'nstate').
        value : any
            New value.

        Example:
            manager.set_namd_config("md", "step", 200)
            manager.set_namd_config("md", "temp", 350)
            manager.set_namd_config("molecule", "ci", 8)
            manager.set_namd_config("openqp", "threads", 80)
        """
        if section not in self._namd_params:
            self._namd_params[section] = {}
        self._namd_params[section][key] = value

        self._generated["namd"] = _namd_dict_to_control(self._namd_params, self.title)
        print(f"NAMD config updated: &{section.upper()} {key} = {value}")

    def set_xyz(self, xyz):
        text, _ = _resolve_content(xyz)
        self._generated["xyz"] = text
        self._xyz = xyz
        if self.mode == "qp":
            self._generated["openqp"] = _sections_to_inp(self._sections, xyz=xyz)
        print("XYZ content updated.")

    def set_pdb(self, pdb):
        text, _ = _resolve_content(pdb)
        self._generated["pdb"] = text
        self._pdb = pdb
        print("PDB content updated.")

    # ================================================================
    #  Engine setup
    # ================================================================
    def _setup_engines(self):
        if self.mode in {"qp", "qmmm", "namd", "namd_qmmm"}:
            self.openqp = OpenQPEngine(self.configs.get("openqp"))
        if self.mode in {"namd", "namd_qmmm"}:
            self.pyrai2md = PyRAI2MDEngine(self.mode, self.configs.get("pyrai2md"))

    def _setup_files(self):
        """Write all generated inputs to disk."""
        openqp_content = self._generated.get("openqp", "")
        xyz_content = self._generated.get("xyz", "")
        pdb_content = self._generated.get("pdb", "")

        folder = Path(self.folder)
        folder.mkdir(parents=True, exist_ok=True)

        # PDB
        if pdb_content and self.mode in {"qmmm", "namd_qmmm"}:
            pdb_filename = self._sections.get("qmmm", {}).get("pdb_file", f"{self.title}.pdb")
            (folder / pdb_filename).write_text(pdb_content.strip() + "\n", encoding="utf-8")

        # QP mode
        if self.mode == "qp" and openqp_content:
            self.openqp.setup_input(
                folder=self.folder, openqp_content=openqp_content, title=self.title,
            )

        # QMMM mode
        elif self.mode == "qmmm" and openqp_content:
            if xyz_content:
                self.openqp.setup_input_with_xyz(
                    folder=self.folder, openqp_content=openqp_content,
                    xyz_content=xyz_content, title=self.title,
                )
            else:
                self.openqp.setup_input(
                    folder=self.folder, openqp_content=openqp_content, title=self.title,
                )

        # NAMD / NAMD_QMMM: write title, title.openqp, title.xyz
        if self.mode in {"namd", "namd_qmmm"}:
            namd_content = self._generated.get("namd", "")
            self.pyrai2md.setup_files(
                folder=self.folder, input_content=namd_content,
                openqp_content=openqp_content, xyz_content=xyz_content,
                title=self.title,
            )

    # ================================================================
    #  Prepare & Run
    # ================================================================
    def prepare(self):
        if self.openqp is not None:
            self.openqp.prepare()
        if self.pyrai2md is not None:
            self.pyrai2md.prepare()

    def run(self):
        """Rewrite files from current config, then execute."""
        self._setup_files()
        self.prepare()
        if self.mode == "qp":
            return self.openqp.run()
        elif self.mode == "qmmm":
            return self.openqp.run_qmmm()
        elif self.mode == "namd":
            self.pyrai2md.run(folder=self.folder, title=self.title)
        elif self.mode == "namd_qmmm":
            self.pyrai2md.run(folder=self.folder, title=self.title)
