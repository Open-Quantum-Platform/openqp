"""High-level Python helpers for OpenQP calculations."""

from copy import deepcopy
from collections.abc import Mapping

from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
from oqp.pyoqp import Runner
from oqp.utils.geometry import (
    GeometryLookupError,
    builtin_geometry,
    geometry_from_sdf,
    get_geometry,
    normalize_system,
    pubchem_geometry,
)
from oqp.utils.input_parser import OQPConfigParser
from oqp.utils.kword_map import resolve_param_key


def dump_strings_from_parser(parser):
    """Extract a pure string dict from parser."""
    out = {}
    for sec in parser.sections():
        out[sec] = {}
        for opt, val in parser[sec].items():
            out[sec][opt] = val
    return out


OPTIMIZER_RUNTYPES = {"optimize", "meci", "mecp", "tci", "mep", "ts", "irc", "neb"}
RUNTYPE_SECTIONS = {
    "grad": "properties",
    "hess": "hess",
    "nac": "nac",
    "nacme": "nac",
    "ekt": "ekt",
    "prop": "properties",
}


class _SectionProxy:
    """Attribute/call proxy for a section in the OpenQP input schema."""

    def __init__(self, owner, section):
        object.__setattr__(self, "_owner", owner)
        object.__setattr__(self, "_section", section)

    def __call__(self, **kwargs):
        if self._section == "optimize":
            self._owner._set_optimize_options(**kwargs)
            return self._owner
        self._owner.section(self._section, **kwargs)
        return self._owner

    def __getattr__(self, option):
        schema = OQP_CONFIG_SCHEMA.get(self._section, {})
        if option in schema:
            return self._owner.config_typed.get(self._section, {}).get(option)
        if self._section == "optimize":
            backend = self._owner._optimizer_backend_section()
            backend_schema = OQP_CONFIG_SCHEMA.get(backend, {})
            if option in backend_schema:
                return self._owner.config_typed.get(backend, {}).get(option)
        raise AttributeError(f"Unknown OpenQP option '{self._section}.{option}'.")

    def __setattr__(self, option, value):
        if option.startswith("_"):
            object.__setattr__(self, option, value)
            return
        if self._section == "optimize":
            self._owner._set_optimize_options(**{option: value})
            return
        self._owner.section(self._section, **{option: value})


class _WorkflowSectionProxy(_SectionProxy):
    """Compatibility section proxy that can also select a workflow runtype."""

    def __init__(self, owner, section, runtype=None):
        super().__init__(owner, section)
        object.__setattr__(self, "_runtype", runtype)

    def __call__(self, **kwargs):
        runtype = object.__getattribute__(self, "_runtype")
        if runtype is not None:
            self._owner._control(runtype=runtype)
        return super().__call__(**kwargs)


class _WorkflowOptimizeProxy(_WorkflowSectionProxy):
    """Optimization workflow proxy with runtype-aware backend routing."""

    def __init__(self, owner, default_runtype="optimize"):
        super().__init__(owner, "optimize")
        object.__setattr__(self, "_default_runtype", default_runtype)

    def __call__(self, runtype=None, **kwargs):
        if runtype is None:
            runtype = object.__getattribute__(self, "_default_runtype")
        if runtype is not None:
            self._owner._control(runtype=runtype)
        if kwargs:
            self._owner._set_optimize_options(**kwargs)
        return self._owner


class _WorkflowGradientProxy(_WorkflowSectionProxy):
    """Gradient workflow proxy that stores state selection in [properties]."""

    def __init__(self, owner):
        super().__init__(owner, "properties", runtype="grad")

    def __call__(self, state=None, grad=None, **kwargs):
        self._owner._control(runtype="grad")
        return self._owner._set_gradient_options(
            state=state,
            grad=grad,
            **kwargs,
        )


class _WorkflowPcmProxy(_WorkflowSectionProxy):
    """PCM workflow proxy for the current energy-only PCM/ddX path."""

    def __init__(self, owner):
        super().__init__(owner, "pcm", runtype="energy")

    def __call__(self, **kwargs):
        self._owner._require_reference_scf_theory_for("PCM/ddX")
        backend = str(kwargs.get("backend", "ddx")).lower()
        mode = str(kwargs.get("mode", "reference_scf")).lower()
        scf_type = str(self._owner.config_typed.get("scf", {}).get("type", "rhf")).lower()
        if backend != "ddx":
            raise ValueError("PCM/ddX currently supports backend='ddx' only.")
        if mode != "reference_scf":
            raise ValueError("PCM/ddX currently supports mode='reference_scf' only.")
        if scf_type not in {"rhf", "rohf"}:
            raise ValueError("PCM/ddX currently supports RHF/ROHF reference SCF only.")
        return super().__call__(**kwargs)


class _WorkflowNmrProxy:
    """NMR shielding workflow proxy for HF/DFT reference-SCF theory."""

    def __init__(self, owner):
        self._owner = owner

    def __call__(self, gauge=None, **kwargs):
        self._owner._require_reference_scf_theory_for("NMR")
        if gauge is not None:
            kwargs["nmr_gauge"] = gauge
        nmr_gauge = str(kwargs.get("nmr_gauge", "cgo")).lower()
        scf_type = str(self._owner.config_typed.get("scf", {}).get("type", "rhf")).lower()
        if nmr_gauge not in {"cgo", "giao"}:
            raise ValueError("NMR gauge must be 'cgo' or 'giao'.")
        if nmr_gauge == "cgo" and scf_type != "rhf":
            raise ValueError("CGO NMR shielding supports closed-shell RHF references only.")
        self._owner._reject_nmr_unsupported_functionals()
        kwargs["scf_prop"] = "nmr"
        return self._owner.section("properties", **kwargs)


class _WorkflowMrsfSectionProxy(_WorkflowSectionProxy):
    """Workflow section proxy for features currently limited to MRSF-TDDFT."""

    def __init__(self, owner, section, runtype, workflow_name):
        super().__init__(owner, section, runtype=runtype)
        object.__setattr__(self, "_workflow_name", workflow_name)

    def __call__(self, **kwargs):
        self._owner._require_mrsf_theory_for(object.__getattribute__(self, "_workflow_name"))
        return super().__call__(**kwargs)


class _WorkflowEktProxy(_WorkflowMrsfSectionProxy):
    """MRSF-EKT workflow proxy."""

    def __init__(self, owner):
        super().__init__(owner, "ekt", runtype="ekt", workflow_name="EKT")

    def __call__(self, **kwargs):
        if not kwargs.get("ip") and not kwargs.get("ea"):
            raise ValueError("EKT requires ip=True, ea=True, or both.")
        return super().__call__(**kwargs)


class _WorkflowSocProxy:
    """SOC workflow proxy grouped under job.workflow."""

    def __init__(self, owner):
        self._owner = owner

    def __call__(self, **kwargs):
        return self._owner._soc_control(**kwargs)


class _WorkflowEnergyProxy:
    """Plain single-point energy workflow selector."""

    def __init__(self, owner):
        self._owner = owner

    def __call__(self):
        return self._owner._control(runtype="energy")


class _TheoryProxy:
    """Quantum-theory namespace for OpenQP Python scripts."""

    def __init__(self, owner):
        self._owner = owner

    def __call__(self, method, functional=None, basis=None, runtype=None,
                 nstate=3, reference=None, **keywords):
        return self._owner._theory(
            method,
            functional=functional,
            basis=basis,
            runtype=runtype,
            nstate=nstate,
            reference=reference,
            **keywords,
        )

    def hf(self, **kwargs):
        return self._owner.hf(**kwargs)

    def dft(self, functional=None, **kwargs):
        if functional is None:
            raise ValueError("DFT theory requires functional=...")
        return self._owner.dft(functional, **kwargs)

    def mp2(self, **kwargs):
        return self._owner.mp2(**kwargs)

    def tdhf(self, **kwargs):
        return self._owner._theory("tdhf", **kwargs)

    def dftb(self, **kwargs):
        return self._owner.dftb(**kwargs)

    def tddft(self, functional=None, **kwargs):
        if functional is None:
            raise ValueError("TDDFT theory requires functional=...")
        return self._owner.tddft(functional, **kwargs)

    def sf_tddft(self, functional=None, **kwargs):
        if functional is None:
            raise ValueError("SF-TDDFT theory requires functional=...")
        return self._owner.sf_tddft(functional, **kwargs)

    def sf(self, functional=None, **kwargs):
        return self.sf_tddft(functional=functional, **kwargs)

    def mrsf(self, **kwargs):
        return self._owner.mrsf(**kwargs)

    def mrsf_tddft(self, **kwargs):
        return self.mrsf(**kwargs)


class _WorkflowProxy:
    """Scientific workflow namespace for OpenQP Python scripts."""

    def __init__(self, owner):
        object.__setattr__(self, "_owner", owner)
        object.__setattr__(self, "energy", _WorkflowEnergyProxy(owner))
        object.__setattr__(self, "gradient", _WorkflowGradientProxy(owner))
        object.__setattr__(self, "hess", _WorkflowSectionProxy(owner, "hess", runtype="hess"))
        object.__setattr__(self, "hessian", object.__getattribute__(self, "hess"))
        object.__setattr__(self, "optimize", _WorkflowOptimizeProxy(owner))
        for runtype in ("meci", "mecp", "tci", "mep", "ts", "irc", "neb"):
            object.__setattr__(self, runtype, _WorkflowOptimizeProxy(owner, runtype))
        object.__setattr__(self, "nac", _WorkflowMrsfSectionProxy(owner, "nac", "nac", "NAC"))
        object.__setattr__(self, "nacme", _WorkflowMrsfSectionProxy(owner, "nac", "nacme", "NACME"))
        object.__setattr__(self, "ekt", _WorkflowEktProxy(owner))
        object.__setattr__(self, "pcm", _WorkflowPcmProxy(owner))
        object.__setattr__(self, "nmr", _WorkflowNmrProxy(owner))
        object.__setattr__(self, "soc", _WorkflowSocProxy(owner))
        # Nonadiabatic MD (Tully surface hopping); MRSF-TDDFT only, [md] section.
        # Gas-phase by default; combine with job.qmmm(...) for QM/MM NAMD and
        # pass soc=True (with optional soc_basis) for SOC-NAMD.
        object.__setattr__(self, "namd", _WorkflowMrsfSectionProxy(owner, "md", "namd", "NAMD"))

    def __call__(self, runtype=None, **kwargs):
        return self._owner._control(
            runtype=runtype,
            **kwargs,
        )

    def __getattr__(self, name):
        if name in OQP_CONFIG_SCHEMA:
            return _WorkflowSectionProxy(self._owner, name)
        raise AttributeError(f"Unknown OpenQP workflow section '{name}'.")


class _SettingsProxy:
    """Raw OpenQP input-section namespace for advanced keyword setup."""

    def __init__(self, owner):
        object.__setattr__(self, "_owner", owner)

    def __call__(self, **sections):
        return self._owner.update(sections)

    def basis(self, basis=None, **tags):
        return self._owner._set_atom_basis(basis, **tags)

    def atom_basis(self, basis=None, **tags):
        return self.basis(basis, **tags)

    def __getattr__(self, name):
        if name in OQP_CONFIG_SCHEMA:
            return _SectionProxy(self._owner, name)
        raise AttributeError(f"Unknown OpenQP settings section '{name}'.")


class _ControlProxy:
    """Hardware/runtime control namespace for OpenQP Python scripts."""

    def __init__(self, owner):
        object.__setattr__(self, "_owner", owner)

    def __call__(self, runtype=None, omp_threads=None, usempi=None, **kwargs):
        return self._owner._control(
            runtype=runtype,
            omp_threads=omp_threads,
            usempi=usempi,
            **kwargs,
        )


class OpenQP:
    """
    OpenQP-native convenience layer over the OpenQP input schema.

    This class is additive: it builds the same sectioned input dictionary used
    by Runner and existing OpenQP input files, while making section-style
    keyword editing concise in Python scripts.
    """

    def __init__(
        self,
        project="oqp_project",
        log=None,
        silent=0,
        usempi=True,
        config=None,
        **sections,
    ):
        self.project = project or "oqp_project"
        self.log = log if log is not None else f"{self.project}.log"
        self.silent = silent
        self.usempi = usempi
        self.unit = "Angstrom"
        self.runner = None
        self.mol = None
        self.theory = _TheoryProxy(self)
        self.workflow = _WorkflowProxy(self)
        self.settings = _SettingsProxy(self)
        self.control = _ControlProxy(self)

        parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA)
        self.config_str = dump_strings_from_parser(parser)
        self.config_typed = parser.validate()

        if config:
            self.update(config)
        if sections:
            self.update(sections)

    def __getattr__(self, name):
        if name in OQP_CONFIG_SCHEMA:
            return _SectionProxy(self, name)
        raise AttributeError(f"{self.__class__.__name__!s} has no attribute {name!r}.")

    @classmethod
    def from_pyscf(cls, mol, **kwargs):
        """
        Build an OpenQP job from a PySCF-like Mole object without importing PySCF.

        PySCF's spin convention is 2S, so spin=2 maps to multiplicity=3.
        """
        runtime_keys = {"project", "log", "silent", "usempi", "config"}
        runtime = {key: kwargs.pop(key) for key in list(kwargs) if key in runtime_keys}
        job = cls(**runtime)

        system = getattr(mol, "atom", None)
        unit = getattr(mol, "unit", "Angstrom")
        input_updates = {}
        for attr in ("basis", "charge"):
            if hasattr(mol, attr):
                value = getattr(mol, attr)
                if value is not None:
                    input_updates[attr] = value

        if system is not None:
            job.molecule(system, unit=unit, **input_updates)
        elif input_updates:
            job.input(**input_updates)

        if hasattr(mol, "spin"):
            spin = getattr(mol, "spin")
            if spin is not None:
                job.scf(multiplicity=int(spin) + 1)

        if kwargs:
            job.update(kwargs)
        return job

    def molecule(self, system=None, system2=None, basis=None, charge=None,
                 multiplicity=None, unit="Angstrom", geometry=None,
                 geometry2=None, source="auto", timeout=10, **kwargs):
        """Set molecular system data using OpenQP input-section keywords."""
        if system is not None and geometry is not None:
            raise ValueError("Use either system or geometry, not both.")
        if system2 is not None and geometry2 is not None:
            raise ValueError("Use either system2 or geometry2, not both.")
        if system2 is not None and basis is None and self._looks_like_basis(system2):
            basis = system2
            system2 = None
        if geometry is not None:
            system = get_geometry(geometry, source=source, timeout=timeout)
        if geometry2 is not None:
            system2 = get_geometry(geometry2, source=source, timeout=timeout)
        if system is None:
            raise ValueError("molecule requires a system or geometry value.")

        self.unit = unit
        updates = {"input.system": system}
        if system2 is not None:
            updates["input.system2"] = system2
        if basis is not None:
            updates["input.basis"] = basis
        if charge is not None:
            updates["input.charge"] = charge
        if multiplicity is not None:
            updates["scf.multiplicity"] = multiplicity
        for option, value in kwargs.items():
            updates[f"input.{option}"] = value
        return self.set(**updates)

    def qmmm(self, pdb_file=None, forcefield=None, forcefield_files=None,
             qm_atoms=None, cutoff=None, embedding=None, rigidwater=None,
             frontier_scheme=None, **kwargs):
        """Enable ESPF QM/MM embedding and configure the ``[qmmm]`` section.

        Sets ``[input] qmmm_flag=true`` and populates ``[qmmm]``. The QM geometry
        and atom selection come from ``job.molecule("file.pdb <indices>")`` for a
        single-point QM/MM energy; QM/MM molecular dynamics (``runtype=namd``)
        additionally reads ``pdb_file``/``forcefield``/``qm_atoms`` here. Combine
        with ``job.workflow.namd(...)`` for (SOC-)NAMD-QMMM dynamics.

        ``forcefield`` is an alias for the ``[qmmm] forcefield_files`` list.
        ``qm_atoms`` accepts a string (``"0-2"`` / ``"0 1 2"``) or a list of
        indices; a ``forcefield`` list is joined into a comma-separated string.

        ``frontier_scheme`` sets the MM frontier-host (M1) charge treatment when
        the QM/MM partition cuts a covalent bond in the ESPF path:
        ``'none'`` (default, the validated full-field ESPF baseline), or the
        optional redistributions ``'rcd'`` / ``'rc'`` / ``'z1'``. It is a no-op
        for whole-molecule QM regions. Covalent QM/MM boundaries are handled by
        the ground-state QM/MM MD path (``QMMM_MD``); the nonadiabatic
        ``runtype=namd`` path does not yet append link atoms to its QM molecule
        and raises on a covalent cut.

        Any other ``[qmmm]`` keyword can be passed through as a keyword argument.
        """
        if forcefield is not None and forcefield_files is not None:
            raise ValueError("Use either forcefield or forcefield_files, not both.")
        ff = forcefield if forcefield is not None else forcefield_files
        if isinstance(ff, (list, tuple)):
            ff = ",".join(str(item) for item in ff)
        if isinstance(qm_atoms, (list, tuple)):
            qm_atoms = " ".join(str(index) for index in qm_atoms)

        self.set(**{"input.qmmm_flag": True})
        updates = {}
        if pdb_file is not None:
            updates["pdb_file"] = pdb_file
        if ff is not None:
            updates["forcefield_files"] = ff
        if qm_atoms is not None:
            updates["qm_atoms"] = qm_atoms
        if cutoff is not None:
            updates["cutoff"] = cutoff
        if embedding is not None:
            updates["embedding"] = embedding
        if rigidwater is not None:
            updates["rigidwater"] = rigidwater
        if frontier_scheme is not None:
            updates["frontier_scheme"] = frontier_scheme
        updates.update(kwargs)
        if updates:
            self.section("qmmm", **updates)
        return self

    @staticmethod
    def _looks_like_basis(value):
        if not isinstance(value, str):
            return False
        text = value.strip()
        if text.lower().endswith((".xyz", ".mol", ".sdf", ".pdb")):
            return False
        return bool(text and "\n" not in text and ";" not in text and " " not in text)

    def _control(self, runtype=None, omp_threads=None, usempi=None, **kwargs):
        """Set run-level controls such as runtype, OpenMP, and optimization options."""
        updates = {}
        section_runtype = runtype
        if runtype is not None and str(runtype).lower() == "pcm":
            section_runtype = "pcm"
            runtype = "energy"
        if runtype is not None:
            updates["input.runtype"] = runtype
        if omp_threads is not None:
            updates["input.omp_threads"] = omp_threads
        if usempi is not None:
            self.usempi = self._as_bool(usempi, option="usempi")
        if updates:
            self.set(**updates)

        if kwargs:
            active_runtype = str(
                section_runtype if section_runtype is not None
                else self.config_typed.get("input", {}).get("runtype", "energy")
            ).lower()
            if active_runtype in OPTIMIZER_RUNTYPES:
                return self._set_optimize_options(**kwargs)
            if active_runtype == "soc":
                return self._soc_control(**kwargs)
            if active_runtype == "pcm":
                return self.section("pcm", **kwargs)
            if active_runtype == "grad":
                return self._set_gradient_options(**kwargs)
            if active_runtype in RUNTYPE_SECTIONS:
                return self.section(RUNTYPE_SECTIONS[active_runtype], **kwargs)
            raise KeyError(
                "Extra job.control(...) options are supported for known "
                "workflow runtypes. Use job.workflow.<name>(...) for "
                "workflow-specific keywords."
            )
        return self

    @staticmethod
    def _as_bool(value, option):
        """Parse common Python/string booleans for runtime-only controls."""
        if isinstance(value, str):
            normalized = value.strip().lower()
            if normalized in {"1", "true", "yes", "on"}:
                return True
            if normalized in {"0", "false", "no", "off"}:
                return False
            raise ValueError(f"{option} must be a boolean value.")
        return bool(value)

    def _set_gradient_options(self, state=None, grad=None, **kwargs):
        """Set gradient state controls while keeping legacy grad= accepted."""
        if state is not None and grad is not None:
            raise ValueError("Use either state=... or legacy grad=..., not both.")
        if state is not None:
            kwargs["grad"] = state
        elif grad is not None:
            kwargs["grad"] = grad
        return self.section("properties", **kwargs)

    def theory(self, method, functional=None, basis=None, runtype=None,
               nstate=3, reference=None, **keywords):
        """Backward-compatible direct theory dispatcher."""
        return self._theory(
            method,
            functional=functional,
            basis=basis,
            runtype=runtype,
            nstate=nstate,
            reference=reference,
            **keywords,
        )

    def _theory(self, method, functional=None, basis=None, runtype=None,
                nstate=3, reference=None, **keywords):
        """Set a compact OpenQP theory model."""
        method_key = str(method).lower().replace("_", "-")
        if method_key in {"hf", "hartree-fock"}:
            return self.hf(
                reference=reference or "rhf",
                runtype=runtype,
                basis=basis,
                **keywords,
            )
        if method_key in {"dft", "ks", "kohn-sham"}:
            if functional is None:
                raise ValueError("DFT theory requires functional=...")
            return self.dft(
                functional,
                reference=reference or "rhf",
                runtype=runtype,
                basis=basis,
                **keywords,
            )
        if method_key in {"mp2", "moller-plesset", "moller-plesset-2"}:
            if functional not in (None, ""):
                raise ValueError("MP2 theory requires an HF reference; do not pass functional.")
            return self.mp2(
                reference=reference or "rhf",
                runtype=runtype,
                basis=basis,
                **keywords,
            )
        if method_key in {"tdhf", "td-hf"}:
            multiplicity = keywords.pop("multiplicity", 1)
            return self._response_theory(
                functional="",
                basis=basis,
                runtype=runtype,
                nstate=nstate,
                reference=reference or "rhf",
                multiplicity=multiplicity,
                **keywords,
            )
        if method_key in {"tddft", "td-dft"}:
            if functional is None:
                raise ValueError("TDDFT theory requires functional=...")
            multiplicity = keywords.pop("multiplicity", 1)
            return self._response_theory(
                functional=functional,
                basis=basis,
                runtype=runtype,
                nstate=nstate,
                reference=reference or "rhf",
                multiplicity=multiplicity,
                **keywords,
            )
        if method_key in {"dftb", "openqp-dftb", "tddftb", "td-dftb"}:
            return self.dftb(
                runtype=runtype,
                response_type=keywords.pop("response_type", "tddftb"),
                nstate=nstate,
                **keywords,
            )
        if method_key in {"sf-tddftb", "sf-td-dftb", "sftddftb"}:
            return self.dftb(
                runtype=runtype,
                response_type="sf",
                nstate=nstate,
                **keywords,
            )
        if method_key in {"mrsf-tddftb", "mrsf-td-dftb", "mrsftddftb"}:
            return self.dftb(
                runtype=runtype,
                response_type="mrsf",
                nstate=nstate,
                **keywords,
            )
        if method_key in {"sf-tddft", "sf-td-dft", "sftddft"}:
            if functional is None:
                raise ValueError("SF-TDDFT theory requires functional=...")
            multiplicity = keywords.pop("multiplicity", 3)
            return self._response_theory(
                functional=functional,
                basis=basis,
                runtype=runtype,
                nstate=nstate,
                reference=reference or "rohf",
                multiplicity=multiplicity,
                response_type="sf",
                **keywords,
            )
        if method_key in {"mrsf", "mrsf-tddft", "mrsf-td-dft"}:
            return self.mrsf(
                nstate=nstate,
                reference=reference or "rohf",
                runtype=runtype,
                functional=functional,
                basis=basis,
                **keywords,
            )
        raise ValueError(
            "Unknown theory method. Use hf, dft, mp2, tdhf, tddft, "
            "sf-tddft, mrsf-tddft, or dftb."
        )

    def hf(self, reference="rhf", runtype=None, multiplicity=None,
           basis=None, **scf_keywords):
        """Use a compact OpenQP HF setup for ordinary single-reference jobs."""
        # Clear any functional left from a prior DFT setup; OpenQP switches to
        # DFT whenever input.functional is non-empty, so HF must reset it.
        input_updates = {"method": "hf", "functional": ""}
        if runtype is not None:
            input_updates["runtype"] = runtype
        if basis is not None:
            input_updates["basis"] = basis
        self.input(**input_updates)
        updates = {}
        if reference is not None:
            updates["type"] = reference
        if multiplicity is not None:
            updates["multiplicity"] = multiplicity
        updates.update(scf_keywords)
        if updates:
            self.scf(**updates)
        return self

    def dft(self, functional, reference="rhf", runtype=None,
            multiplicity=None, basis=None, **scf_keywords):
        """Use a compact OpenQP DFT setup for ordinary Kohn-Sham jobs."""
        input_updates = {
            "method": "hf",
            "functional": functional,
        }
        if runtype is not None:
            input_updates["runtype"] = runtype
        if basis is not None:
            input_updates["basis"] = basis
        self.input(**input_updates)
        updates = {}
        if reference is not None:
            updates["type"] = reference
        if multiplicity is not None:
            updates["multiplicity"] = multiplicity
        updates.update(scf_keywords)
        if updates:
            self.scf(**updates)
        return self

    def mp2(self, reference="rhf", runtype=None, multiplicity=None,
            basis=None, variant=None, same_spin_scale=None,
            opposite_spin_scale=None, **scf_keywords):
        """Use a compact OpenQP MP2 setup for energy-only post-SCF jobs."""
        if runtype is None:
            runtype = "energy"
        elif str(runtype).lower() != "energy":
            raise ValueError("MP2 currently supports runtype='energy' only.")
        if "functional" in scf_keywords:
            functional = scf_keywords.pop("functional")
            if functional:
                raise ValueError("MP2 theory requires an HF reference; do not pass functional.")

        has_custom_scale = same_spin_scale is not None or opposite_spin_scale is not None
        if has_custom_scale:
            if variant is None:
                variant = "custom"
            elif str(variant).lower() != "custom":
                raise ValueError("Custom MP2 scale factors require variant='custom'.")

        input_updates = {"method": "mp2", "functional": "", "runtype": runtype}
        if basis is not None:
            input_updates["basis"] = basis
        self.input(**input_updates)

        scf_updates = {}
        if reference is not None:
            scf_updates["type"] = reference
        if multiplicity is not None:
            scf_updates["multiplicity"] = multiplicity
        scf_updates.update(scf_keywords)
        if scf_updates:
            self.scf(**scf_updates)

        mp2_updates = {}
        if variant is not None:
            mp2_updates["variant"] = variant
        if same_spin_scale is not None:
            mp2_updates["same_spin_scale"] = same_spin_scale
        if opposite_spin_scale is not None:
            mp2_updates["opposite_spin_scale"] = opposite_spin_scale
        if mp2_updates:
            self.section("mp2", **mp2_updates)
        return self

    def _response_theory(self, functional="", basis=None, runtype=None,
                         nstate=3, reference="rhf", multiplicity=1,
                         response_type=None, **tdhf_keywords):
        """Use a compact OpenQP TDHF/TDDFT response setup."""
        input_updates = {"method": "tdhf", "functional": functional or ""}
        if runtype is not None:
            input_updates["runtype"] = runtype
        if basis is not None:
            input_updates["basis"] = basis
        self.input(**input_updates)

        scf_updates = {}
        if reference is not None:
            scf_updates["type"] = reference
        if multiplicity is not None:
            scf_updates["multiplicity"] = multiplicity
        if scf_updates:
            self.scf(**scf_updates)

        updates = {"nstate": nstate}
        if response_type is not None:
            requested_type = tdhf_keywords.pop("type", response_type)
            if str(requested_type).lower() != str(response_type).lower():
                raise ValueError(
                    f"This theory helper requires [tdhf] type={response_type}."
                )
            updates["type"] = requested_type
        updates.update(tdhf_keywords)
        return self.section("tdhf", **updates)

    def tddft(self, functional, reference="rhf", runtype=None,
              multiplicity=1, basis=None, nstate=3, **tdhf_keywords):
        """Use a compact OpenQP TDDFT setup."""
        return self._response_theory(
            functional=functional,
            basis=basis,
            runtype=runtype,
            nstate=nstate,
            reference=reference,
            multiplicity=multiplicity,
            **tdhf_keywords,
        )

    def sf_tddft(self, functional, reference="rohf", runtype=None,
                 multiplicity=3, basis=None, nstate=3, **tdhf_keywords):
        """Use a compact OpenQP spin-flip TDDFT setup."""
        return self._response_theory(
            functional=functional,
            basis=basis,
            runtype=runtype,
            nstate=nstate,
            reference=reference,
            multiplicity=multiplicity,
            response_type="sf",
            **tdhf_keywords,
        )

    def mrsf(self, nstate=3, reference="rohf", multiplicity=3,
             runtype=None, functional=None, basis=None, **tdhf_keywords):
        """Use a compact OpenQP MRSF-TDDFT setup with an optional functional."""
        input_updates = {"method": "tdhf"}
        if runtype is not None:
            input_updates["runtype"] = runtype
        if functional is not None:
            input_updates["functional"] = functional
        if basis is not None:
            input_updates["basis"] = basis
        self.input(**input_updates)
        self.scf(type=reference, multiplicity=multiplicity)
        updates = {"type": "mrsf", "nstate": nstate}
        updates.update(tdhf_keywords)
        return self.tdhf(**updates)

    def dftb(self, runtype=None, response_type="mrsf", nstate=3,
             parameter_path=None, **keywords):
        """Use the optional OpenQP-DFTB backend through the normal OpenQP workflow."""
        input_updates = {"method": "dftb", "functional": ""}
        if runtype is not None:
            input_updates["runtype"] = runtype
        self.input(**input_updates)

        dftb_schema = OQP_CONFIG_SCHEMA.get("dftb", {})
        dftb_updates = {}
        if parameter_path is not None:
            dftb_updates["parameter_path"] = parameter_path
        for key in list(keywords.keys()):
            if key in dftb_schema:
                dftb_updates[key] = keywords.pop(key)

        requested_type = str(keywords.pop("type", response_type)).lower()
        dftb_type = requested_type
        tdhf_type = {
            "ground": "tda",
            "dftb": "tda",
            "dftb0": "tda",
            "tddftb": "tda",
            "td-dftb": "tda",
            "tda": "tda",
            "sf": "sf",
            "sftddftb": "sf",
            "sf-tddftb": "sf",
            "mrsf": "mrsf",
            "mrsftddftb": "mrsf",
            "mrsf-tddftb": "mrsf",
        }.get(requested_type)
        if tdhf_type is None:
            raise ValueError("DFTB response_type must be ground, tddftb, sf, or mrsf.")

        # The backend's response-method names are ground/tddftb/sf/mrsf; map the
        # tda/td-dftb aliases onto tddftb so the explicit request and the auto
        # path (tdhf.type=tda/rpa -> tddftb) reach the same backend method.
        _BACKEND_TYPE = {"tda": "tddftb", "td-dftb": "tddftb"}
        dftb_updates["type"] = _BACKEND_TYPE.get(dftb_type, dftb_type)
        self.section("dftb", **dftb_updates)
        return self.tdhf(type=tdhf_type, nstate=nstate, **keywords)

    def soc(self, nstate=3, functional=None, reference="rohf",
            reference_multiplicity=3, soc_2e=1, scal_rel=2,
            basis=None, **tdhf_keywords):
        """Use a compact MRSF-TDDFT SOC setup."""
        return self._soc(
            nstate=nstate,
            functional=functional,
            reference=reference,
            reference_multiplicity=reference_multiplicity,
            soc_2e=soc_2e,
            scal_rel=scal_rel,
            basis=basis,
            **tdhf_keywords,
        )

    def _soc_control(self, soc_2e=1, scal_rel=2, **tdhf_keywords):
        """Select the SOC workflow after an MRSF-TDDFT theory has been set."""
        self._require_mrsf_theory_for_soc()
        theory_keys = {"basis", "functional", "reference", "reference_multiplicity"}
        misplaced = sorted(theory_keys & set(tdhf_keywords))
        if misplaced:
            raise ValueError(
                "Set SOC theory with job.theory('mrsf-tddft', ...) before "
                "job.workflow.soc(...). Move these options to job.theory: "
                f"{', '.join(misplaced)}."
            )
        if "multiplicity" in tdhf_keywords:
            raise ValueError(
                "SOC computes singlet and triplet response internally; "
                "do not set tdhf.multiplicity in job.workflow.soc()."
            )

        self.set(**{
            "input.runtype": "soc",
            "input.soc_2e": soc_2e,
            "scf.scal_rel": scal_rel,
        })
        if tdhf_keywords:
            self.tdhf(**tdhf_keywords)
        return self

    def _require_mrsf_theory_for_soc(self):
        self._require_mrsf_theory_for("SOC")

    def _require_mrsf_theory_for(self, workflow_name):
        method = str(self.config_typed.get("input", {}).get("method", "")).lower()
        response = str(self.config_typed.get("tdhf", {}).get("type", "")).lower()
        dftb_type = str(self.config_typed.get("dftb", {}).get("type", "auto")).lower()
        dftb_mrsf = method == "dftb" and dftb_type in {
            "auto", "mrsf", "mrsftddftb", "mrsf-tddftb"} and response == "mrsf"
        if not ((method == "tdhf" and response == "mrsf") or dftb_mrsf):
            raise ValueError(
                f"{workflow_name} is currently supported only with MRSF-TDDFT "
                "or MRSF-TDDFTB. Call job.theory('mrsf-tddft', ...) or "
                "job.dftb(response_type='mrsf', ...) before selecting this workflow."
            )

    def _require_reference_scf_theory_for(self, workflow_name):
        method = str(self.config_typed.get("input", {}).get("method", "")).lower()
        if method != "hf":
            raise ValueError(
                f"{workflow_name} is currently supported only with HF/DFT "
                "reference-SCF theory. Call job.theory('hf', ...) or "
                "job.theory('dft', ...) before selecting this workflow."
            )

    def _reject_nmr_unsupported_functionals(self):
        functional = str(
            self.config_typed.get("input", {}).get("functional", "")
        ).lower()
        if functional.startswith(("cam-", "dtcam-", "lc-", "lrc-", "wb97", "hse")):
            raise ValueError(
                "NMR shielding with range-separated functionals is not implemented."
            )
        if functional.startswith((
            "m06", "m08", "m11", "mn12", "mn15", "tpss",
            "scan", "rscan", "r2scan", "b97m", "revm06",
        )):
            raise ValueError(
                "NMR shielding with meta-GGA functionals is not implemented."
            )

    def _soc(self, nstate=3, functional=None, reference="rohf",
             reference_multiplicity=3, soc_2e=1, scal_rel=2,
             basis=None, **tdhf_keywords):
        """Use a compact MRSF-TDDFT SOC setup."""
        if "multiplicity" in tdhf_keywords:
            raise ValueError(
                "SOC computes singlet and triplet response internally; "
                "do not set tdhf.multiplicity in job.soc()."
            )

        input_updates = {"method": "tdhf", "runtype": "soc", "soc_2e": soc_2e}
        if functional is not None:
            input_updates["functional"] = functional
        if basis is not None:
            input_updates["basis"] = basis
        self.input(**input_updates)
        self.scf(type=reference, multiplicity=reference_multiplicity, scal_rel=scal_rel)
        updates = {"type": "mrsf", "nstate": nstate}
        updates.update(tdhf_keywords)
        return self.tdhf(**updates)

    def _set_optimize_options(self, **kwargs):
        """Set optimizer options, routing backend-specific keys by lib."""
        requested = dict(kwargs)
        lib_value = requested.pop(
            "lib",
            self.config_typed.get("optimize", {}).get("lib", "oqp"),
        )
        lib_name = str(lib_value).lower()
        backend = self._optimizer_backend_section(lib_name)

        optimize_schema = OQP_CONFIG_SCHEMA.get("optimize", {})
        backend_schema = OQP_CONFIG_SCHEMA.get(backend, {})

        updates = {"optimize.lib": lib_name}
        for option, value in requested.items():
            if option in optimize_schema:
                updates[f"optimize.{option}"] = value
            elif option in backend_schema:
                updates[f"{backend}.{option}"] = value
            else:
                valid = sorted(set(optimize_schema) | set(backend_schema))
                raise KeyError(
                    f"Unknown optimizer option '{option}' for lib='{lib_name}'. "
                    f"Valid options are: {', '.join(valid)}"
                )
        return self.set(**updates)

    def _set_atom_basis(self, basis=None, **tags):
        """Set atom-wise basis assignment as a detailed settings operation."""
        if tags:
            if basis is not None:
                raise ValueError("Use either a basis mapping/list or tag keywords, not both.")
            basis = tags
        if basis is None:
            raise ValueError("settings.basis requires atom-wise basis data.")

        if isinstance(basis, Mapping):
            if not basis:
                raise ValueError("settings.basis mapping cannot be empty.")
            library = "\n".join(f"{tag} {name}" for tag, name in basis.items())
            return self.section("input", basis="library", library=library)

        if isinstance(basis, str):
            entries = [entry.strip() for entry in basis.split(";") if entry.strip()]
            if len(entries) == 1:
                raise ValueError(
                    "Use job.theory(..., basis=...) for a single global basis. "
                    "job.settings.basis(...) is for atom-wise basis assignments."
                )
        else:
            try:
                entries = [str(entry).strip() for entry in basis]
            except TypeError as exc:
                raise TypeError(
                    "settings.basis expects a mapping of atom tags to basis names "
                    "or an ordered iterable of per-atom basis names."
                ) from exc

        if not entries:
            raise ValueError("settings.basis entries cannot be empty.")
        if any(not entry for entry in entries):
            raise ValueError("settings.basis entries must be non-empty basis names.")
        return self.section("input", basis=";".join(entries))

    def _optimizer_backend_section(self, lib_name=None):
        """Return the backend-specific option section for the selected optimizer."""
        if lib_name is None:
            lib_name = self.config_typed.get("optimize", {}).get("lib", "oqp")
        lib_name = str(lib_name).lower()
        if lib_name in {"oqp", "geometric"}:
            return lib_name
        return None

    def section(self, section, **kwargs):
        """Update one OpenQP section using keyword arguments."""
        if section not in OQP_CONFIG_SCHEMA:
            raise KeyError(f"Unknown OpenQP section '{section}'.")
        if section == "input" and "unit" in kwargs:
            self.unit = kwargs.pop("unit")
        return self.set(**{f"{section}.{opt}": value for opt, value in kwargs.items()})

    def set(self, **kwargs):
        """
        Update configuration with dotted OpenQP keys or section dictionaries.
        """
        parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA)
        for sec, opts in self.config_str.items():
            if not parser.has_section(sec):
                parser.add_section(sec)
            for opt, sval in opts.items():
                parser[sec][opt] = sval

        flat_updates = {}
        for key, value in kwargs.items():
            if key == "unit":
                self.unit = value
                continue
            if key in OQP_CONFIG_SCHEMA and isinstance(value, Mapping):
                for opt, opt_value in value.items():
                    flat_updates[f"{key}.{opt}"] = opt_value
            else:
                flat_updates[key] = value

        for key, value in flat_updates.items():
            canonical_key, canonical_value = self._canonicalize_key_value(key, value)
            sec, opt = resolve_param_key(canonical_key)
            if not parser.has_section(sec):
                parser.add_section(sec)
            parser[sec][opt] = str(canonical_value)

        self.config_str = dump_strings_from_parser(parser)
        self.config_typed = parser.validate()

        if self.runner is not None:
            self.runner.mol.load_config(self.config_str)
        return self

    def update(self, config=None, **kwargs):
        """Update from a sectioned dictionary plus optional keyword overrides."""
        updates = {}
        if config:
            for section, options in config.items():
                if isinstance(options, Mapping):
                    for option, value in options.items():
                        updates[f"{section}.{option}"] = value
                else:
                    updates[section] = options
        updates.update(kwargs)
        return self.set(**updates)

    def to_input_dict(self):
        """Return a copy of the sectioned input dictionary passed to Runner."""
        return deepcopy(self.config_str)

    def run(self, run_type=None):
        """Run the calculation and return the native OpenQP Molecule object."""
        if run_type:
            self.set(**{"input.runtype": run_type})
        self.runner = Runner(
            project=self.project,
            input_file=None,
            log=self.log,
            input_dict=self.config_str,
            silent=self.silent,
            usempi=self.usempi,
        )
        self.runner.run()
        self.mol = self.runner.mol
        return self.mol

    def _canonicalize_key_value(self, key, value):
        key = str(key)
        if key in ("input.system", "input.system2"):
            value = normalize_system(value, self.unit)
        return key, value


OQP = OpenQP


class OPENQP:
    def __init__(self, cfg: dict, silent: bool = False):
        parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA)
        for k, v in cfg.items():
            if k == "input.system":
                v = self._normalize_system(v)
            sec, opt = resolve_param_key(k)
            if not parser.has_section(sec):
                parser.add_section(sec)
            parser[sec][opt] = str(v)

        self.config_str = self._dump_strings_from_parser(parser)
        self.config_typed = parser.validate()

        self.runner = Runner(
            project="oqp_project",
            input_file=None,
            log="oqp_project.log",
            input_dict=self.config_str,
            silent=1 if silent else 0,
            usempi=True
        )
        self.mol = self.runner.mol

    @property
    def sp(self):
        """SinglePoint view used by the ESPF QM/MM driver
        (oqp.library.qmmm_driver) for embedded SCF / excitation calls
        (``self.op.sp.*``). Built on first access and imported lazily so that
        importing oqp.openqp / constructing OPENQP does not require the compiled
        library."""
        sp = self.__dict__.get("_sp")
        if sp is None:
            from oqp.library.single_point import SinglePoint
            sp = SinglePoint(self.mol)
            self.__dict__["_sp"] = sp
        return sp

    def _normalize_system(self, system):
        """
        Accepts:
          - path/to/file.xyz   -> pass-through if file exists
          - "H 0 0 0; H 0 0 0.74" or ","
          - "H 0 0 0\nH 0 0 0.74"
          - [("H",0,0,0), ("H",0,0,0.74)]
        Produces: "\nH 0 0 0\nH 0 0 0.74"  (NOTE: leading newline, NO trailing newline)
        """
        return normalize_system(system)

    @staticmethod
    def _dump_strings_from_parser(parser: OQPConfigParser):
        """Extract a pure string dict from parser (preserves comma formats for arrays)."""
        return dump_strings_from_parser(parser)

    def set(self, **kwargs):
        """
        Update configuration. Write kwargs as strings into a fresh parser,
        then validate to refresh typed view. Always pass strings to Molecule.
        """
        parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA)

        # load current STRING config back first (so array formats stay correct)
        for sec, opts in self.config_str.items():
            if not parser.has_section(sec):
                parser.add_section(sec)
            for opt, sval in opts.items():
                parser[sec][opt] = sval

        # apply updates (as strings)
        for k, v in kwargs.items():
            sec, opt = resolve_param_key(k)
            if not parser.has_section(sec):
                parser.add_section(sec)
            parser[sec][opt] = str(v)

        # refresh both copies
        self.config_str = self._dump_strings_from_parser(parser)
        self.config_typed = parser.validate()

        # push strings into Molecule
        self.runner.mol.load_config(self.config_str)
        return self

    def run(self, run_type=None):
        if run_type:
            # update both string + typed configs consistently
            self.set(**{"input.runtype": run_type})
        self.runner.run()
        return self.mol
