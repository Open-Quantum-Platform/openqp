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
from oqp.utils.tb_backends import is_tb_method
from oqp.utils.state_labels import canonical_dftb_type


def dump_strings_from_parser(parser):
    """Extract a pure string dict from parser."""
    out = {}
    for sec in parser.sections():
        out[sec] = {}
        for opt, val in parser[sec].items():
            out[sec][opt] = val
    return out


OPTIMIZER_RUNTYPES = {"optimize", "meci", "mecp", "tci", "mep", "ts", "irc", "neb"}
# Spellings `theory()` accepts for the CASPT2 family -> the input `method`.
# Both the hyphenated and the run-together forms, since `[input] method`
# accepts both.
_CASPT2_VARIANTS = {
    "caspt2": "caspt2",
    "ms-caspt2": "ms-caspt2", "mscaspt2": "ms-caspt2",
    "xms-caspt2": "xms-caspt2", "xmscaspt2": "xms-caspt2",
}
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


class _DFTBSectionProxy(_SectionProxy):
    """Callable [dftb] section proxy.

    ``job.dftb(...)`` runs the DFTB workflow helper (method=dftb plus the
    response-type plumbing), while ``job.dftb.option`` reads and
    ``job.dftb.option = value`` writes the [dftb] section through the standard
    schema interface -- a plain method here would shadow the ``__getattr__``
    section proxy and break attribute access for this one section.
    """

    def __call__(self, *args, **kwargs):
        return self._owner._dftb(*args, **kwargs)


class _XTBSectionProxy(_SectionProxy):
    """Callable [xtb] section proxy.

    ``job.xtb(...)`` runs the OpenQP-xTB workflow helper (method=xtb plus the
    response-type plumbing), while ``job.xtb.option`` reads and
    ``job.xtb.option = value`` writes the [xtb] section through the standard
    schema interface -- a plain method here would shadow the ``__getattr__``
    section proxy and break attribute access for this one section.
    """

    def __call__(self, *args, **kwargs):
        return self._owner._xtb(*args, **kwargs)


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


class _WorkflowNAMDProxy(_WorkflowMrsfSectionProxy):
    """NAMD workflow proxy with early same-spin diagnostic validation."""

    def __init__(self, owner):
        super().__init__(owner, "md", runtype="namd", workflow_name="NAMD")

    def __call__(self, **kwargs):
        current = self._owner.config_typed.get("md", {})
        soc = kwargs.get("soc", current.get("soc", False))
        nacme_explicit = "nacme_check" in kwargs
        nacme_check = kwargs.get(
            "nacme_check", current.get("nacme_check", "baeck_an"))
        soc_requested = (soc is True) or (
            str(soc).strip().lower() in {"true", "1", "on", "yes"})
        if (soc_requested and nacme_explicit
                and str(nacme_check).strip().lower() != "off"):
            raise ValueError(
                "SOC NAMD does not support nacme_check; use nacme_check='off'."
            )
        if soc_requested and not nacme_explicit:
            # The global default is the same-spin Baeck-An diagnostic.  SOC
            # stores its complex spin-adiabatic overlap/TDC instead.
            kwargs["nacme_check"] = "off"
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

    def ccsd(self, **kwargs):
        return self._owner.ccsd(**kwargs)

    def ccsd_t(self, **kwargs):
        return self._owner.ccsd_t(**kwargs)

    def tdhf(self, **kwargs):
        return self._owner._theory("tdhf", **kwargs)

    def dftb(self, **kwargs):
        return self._owner._dftb(**kwargs)

    def xtb(self, **kwargs):
        return self._owner._xtb(**kwargs)

    def ground_dftb(self, **kwargs):
        return self._owner.ground_dftb(**kwargs)

    def tddftb(self, **kwargs):
        return self._owner.tddftb(**kwargs)

    def sf_tddftb(self, **kwargs):
        return self._owner.sf_tddftb(**kwargs)

    def mrsf_tddftb(self, **kwargs):
        return self._owner.mrsf_tddftb(**kwargs)


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
        object.__setattr__(self, "namd", _WorkflowNAMDProxy(owner))

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
        and atom selection come from ``job.molecule("file.pdb <indices>")``.
        When ``pdb_file`` or ``qm_atoms`` is omitted, it is inferred from that
        molecular-system value for both single-point and molecular-dynamics
        workflows. Combine with ``job.workflow.namd(...)`` for
        (SOC-)NAMD-QMMM dynamics.

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
        system = self.config_str.get("input", {}).get("system")
        inferred_pdb, inferred_qm_atoms = self._qmmm_selection_from_system(system)
        if pdb_file is None:
            pdb_file = inferred_pdb
        if qm_atoms is None:
            qm_atoms = inferred_qm_atoms

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
    def _qmmm_selection_from_system(system):
        """Return the PDB path and atom selector from ``input.system``."""
        if not isinstance(system, str):
            return None, None
        first_line = next(
            (line.strip() for line in system.splitlines() if line.strip()),
            "",
        )
        if not first_line:
            return None, None
        fields = first_line.split()
        if not fields[0].lower().endswith(".pdb"):
            return None, None
        selector = " ".join(fields[1:]) or None
        return fields[0], selector

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
        if method_key in {"ccsd", "cc"}:
            if functional not in (None, ""):
                raise ValueError(
                    "Coupled cluster requires an HF reference; do not pass functional.")
            return self.ccsd(
                reference=reference or "rhf",
                runtype=runtype,
                basis=basis,
                **keywords,
            )
        if method_key in {"ccsd(t)", "ccsd-t", "ccsdt"}:
            if functional not in (None, ""):
                raise ValueError(
                    "Coupled cluster requires an HF reference; do not pass functional.")
            return self.ccsd_t(
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
        # Wavefunction stack.  `nstate` is not forwarded: its default of 3 is
        # a response-method convention and would silently state-average or
        # over-solve here.  SA-CASSCF is the one method for which it is
        # meaningful, and it takes it as the number of averaged states.
        if method_key in {"fci", "full-ci"}:
            return self.fci(
                runtype=runtype, basis=basis,
                reference=reference or "rhf", **keywords)
        if method_key in {"casci", "cas-ci"}:
            return self.casci(
                runtype=runtype, basis=basis,
                reference=reference or "rhf", **keywords)
        if method_key == "casscf":
            return self.casscf(
                runtype=runtype, basis=basis,
                reference=reference or "rhf", **keywords)
        if method_key in {"sa-casscf", "sacasscf"}:
            return self.sa_casscf(
                runtype=runtype, basis=basis, nstate=keywords.pop("nstate", nstate),
                reference=reference or "rhf", **keywords)
        if method_key in _CASPT2_VARIANTS:
            return self.caspt2(
                variant=keywords.pop("variant", _CASPT2_VARIANTS[method_key]),
                runtype=runtype, basis=basis,
                reference=reference or "rhf", **keywords)
        if method_key in {"nevpt2", "sc-nevpt2", "scnevpt2"}:
            if method_key != "nevpt2":
                keywords.setdefault("contraction", "strong")
            return self.nevpt2(
                runtype=runtype, basis=basis,
                reference=reference or "rhf", **keywords)
        if method_key in {"mrmp2", "mcqdpt2", "xmcqdpt2", "qdpt2"}:
            return self.qdpt2(
                variant=keywords.pop("variant",
                                     "mrmp2" if method_key == "qdpt2" else method_key),
                runtype=runtype, basis=basis,
                reference=reference or "rhf", **keywords)

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
        if method_key in {"dftb", "openqp-dftb"}:
            return self._dftb(
                runtype=runtype,
                response_type=keywords.pop("response_type", "ground"),
                nstate=nstate,
                **keywords,
            )
        if method_key in {"tddftb", "td-dftb"}:
            return self._dftb(
                runtype=runtype,
                response_type=keywords.pop("response_type", "tddftb"),
                nstate=nstate,
                **keywords,
            )
        if method_key in {"sf-tddftb", "sf-td-dftb", "sftddftb"}:
            return self._dftb(
                runtype=runtype,
                response_type="sf",
                nstate=nstate,
                **keywords,
            )
        if method_key in {"mrsf-tddftb", "mrsf-td-dftb", "mrsftddftb"}:
            return self._dftb(
                runtype=runtype,
                response_type="mrsf",
                nstate=nstate,
                **keywords,
            )
        if method_key in {"xtb", "openqp-xtb"}:
            return self._xtb(
                runtype=runtype,
                response_type=keywords.pop("response_type", "ground"),
                nstate=nstate,
                **keywords,
            )
        if method_key in {"tdxtb", "td-xtb"}:
            return self._xtb(
                runtype=runtype,
                response_type=keywords.pop("response_type", "tddftb"),
                nstate=nstate,
                **keywords,
            )
        if method_key in {"sf-xtb", "sf-td-xtb", "sfxtb"}:
            return self._xtb(
                runtype=runtype,
                response_type="sf",
                nstate=nstate,
                **keywords,
            )
        if method_key in {"mrsf-xtb", "mrsf-td-xtb", "mrsfxtb"}:
            return self._xtb(
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
            "sf-tddft, mrsf-tddft, dftb, or xtb."
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

    def ccsd(self, reference="rhf", runtype=None, multiplicity=None,
             basis=None, nfzc=None, conv=None, maxit=None, ndiis=None,
             cholesky=None, cholesky_tol=None, cholesky_direct=None,
             triples=False, **scf_keywords):
        """Compact coupled-cluster setup for energy-only post-SCF jobs.

        reference may be rhf, uhf or rohf.  Open-shell references go through
        the spin-orbital solver, which stores sixteen times the integrals of
        the closed-shell path, so keep those systems small.
        """
        if runtype is None:
            runtype = "energy"
        elif str(runtype).lower() != "energy":
            raise ValueError("Coupled cluster currently supports runtype='energy' only.")
        if "functional" in scf_keywords:
            functional = scf_keywords.pop("functional")
            if functional:
                raise ValueError(
                    "Coupled cluster requires an HF reference; do not pass functional.")

        input_updates = {"method": "ccsd(t)" if triples else "ccsd",
                         "functional": "", "runtype": runtype}
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

        # The factorisation controls belong to [cc] like the rest.  Left out of
        # this list they fall through to **scf_keywords and are applied to
        # [scf], so job.ccsd_t(cholesky=False) failed with an unknown scf
        # keyword rather than configuring the route it names.
        cc_updates = {}
        for key, value in (("nfzc", nfzc), ("conv", conv),
                           ("maxit", maxit), ("ndiis", ndiis),
                           ("cholesky", cholesky),
                           ("cholesky_tol", cholesky_tol),
                           ("cholesky_direct", cholesky_direct)):
            if value is not None:
                cc_updates[key] = value
        if cc_updates:
            self.section("cc", **cc_updates)
        return self

    def ccsd_t(self, **kwargs):
        """CCSD with the perturbative triples correction."""
        kwargs["triples"] = True
        return self.ccsd(**kwargs)

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

    # ------------------------------------------------------------------ wavefunction stack
    # FCI / CASCI / CASSCF and the PT2 families built on them.  These share
    # one setup: an RHF reference in [scf], an active space in [cas], a CI
    # solver in [ci], optionally state averaging in [state_average], and for
    # the PT2 methods a perturbation in [pt2].  `_wf_setup` does that wiring
    # once; the public helpers below only choose the method name and the
    # section defaults that distinguish one from another.
    def _wf_setup(self, method, runtype=None, basis=None, reference="rhf",
                  multiplicity=None, active_electrons=None,
                  active_orbitals=None, frozen_core=None, nroot=None,
                  cas=None, ci=None, casscf=None, state_average=None,
                  pt2=None, **scf_keywords):
        input_updates = {"method": method, "functional": ""}
        if runtype is not None:
            input_updates["runtype"] = runtype
        if basis is not None:
            input_updates["basis"] = basis
        # A functional would silently switch OpenQP to DFT, and the whole
        # stack requires an HF reference; refuse rather than compute the
        # wrong thing.
        if scf_keywords.pop("functional", "") not in ("", None):
            raise ValueError(
                f"{method} requires an HF reference; do not pass functional."
            )
        self.input(**input_updates)

        scf_updates = {}
        if reference is not None:
            scf_updates["type"] = reference
        if multiplicity is not None:
            scf_updates["multiplicity"] = multiplicity
        scf_updates.update(scf_keywords)
        if scf_updates:
            self.scf(**scf_updates)

        cas_updates = dict(cas or {})
        if active_electrons is not None:
            cas_updates["active_electrons"] = active_electrons
        if active_orbitals is not None:
            cas_updates["active_orbitals"] = active_orbitals
        if frozen_core is not None:
            cas_updates["frozen_core"] = frozen_core
        if cas_updates:
            self.section("cas", **cas_updates)

        ci_updates = dict(ci or {})
        if nroot is not None:
            ci_updates["nroot"] = nroot
        if ci_updates:
            self.section("ci", **ci_updates)

        if casscf:
            self.section("casscf", **casscf)
        if state_average:
            self.section("state_average", **state_average)
        if pt2:
            self.section("pt2", **pt2)
        return self

    @staticmethod
    def _require_active_space(method, active_electrons, active_orbitals):
        if active_electrons is None or active_orbitals is None:
            raise ValueError(
                f"{method} requires active_electrons=... and active_orbitals=..."
            )

    def fci(self, nroot=1, frozen_core=None, runtype=None, basis=None,
            reference="rhf", **keywords):
        """Use a compact OpenQP full-CI setup on an RHF reference."""
        return self._wf_setup(
            "fci", runtype=runtype, basis=basis, reference=reference,
            frozen_core=frozen_core, nroot=nroot, **keywords)

    def casci(self, active_electrons=None, active_orbitals=None,
              frozen_core=None, nroot=1, runtype=None, basis=None,
              reference="rhf", **keywords):
        """Use a compact OpenQP CASCI setup (fixed reference orbitals)."""
        self._require_active_space("CASCI", active_electrons, active_orbitals)
        return self._wf_setup(
            "casci", runtype=runtype, basis=basis, reference=reference,
            active_electrons=active_electrons, active_orbitals=active_orbitals,
            frozen_core=frozen_core, nroot=nroot, **keywords)

    def casscf(self, active_electrons=None, active_orbitals=None,
               frozen_core=None, nroot=1, root=None, converger=None,
               hessian=None, max_macro_iterations=None, runtype=None,
               basis=None, reference="rhf", **keywords):
        """Use a compact OpenQP CASSCF setup (orbital + CI optimization)."""
        self._require_active_space("CASSCF", active_electrons, active_orbitals)
        opts = dict(keywords.pop("casscf", None) or {})
        for key, value in (("root", root), ("converger", converger),
                           ("hessian", hessian),
                           ("max_macro_iterations", max_macro_iterations)):
            if value is not None:
                opts[key] = value
        return self._wf_setup(
            "casscf", runtype=runtype, basis=basis, reference=reference,
            active_electrons=active_electrons, active_orbitals=active_orbitals,
            frozen_core=frozen_core, nroot=nroot, casscf=opts, **keywords)

    def sa_casscf(self, active_electrons=None, active_orbitals=None,
                  frozen_core=None, nstate=2, weights=None, target_roots=None,
                  runtype=None, basis=None, reference="rhf", **keywords):
        """Use a compact OpenQP state-averaged CASSCF setup.

        `nstate` is the number of averaged states; it also sets [ci] nroot,
        since every averaged root has to be solved for."""
        self._require_active_space("SA-CASSCF", active_electrons, active_orbitals)
        sa = dict(keywords.pop("state_average", None) or {})
        sa.setdefault("enabled", True)
        sa.setdefault("nstate", nstate)
        if weights is not None:
            sa["weights"] = weights
            sa.setdefault("equal_weights", False)
        if target_roots is not None:
            sa["target_roots"] = target_roots
        return self._wf_setup(
            "sa-casscf", runtype=runtype, basis=basis, reference=reference,
            active_electrons=active_electrons, active_orbitals=active_orbitals,
            frozen_core=frozen_core,
            # target_roots may be noncontiguous or not start at zero, e.g.
            # sa_casscf(nstate=2, target_roots=[1, 2]).  Solving only `nstate`
            # roots then makes validation reject root 2, so the CI has to be
            # sized to the highest root actually requested.
            nroot=keywords.pop("nroot", None)
                  or max([nstate] + [int(r) + 1 for r in (target_roots or ())]),
            state_average=sa, **keywords)

    def caspt2(self, active_electrons=None, active_orbitals=None,
               frozen_core=None, nroot=1, variant=None, h0=None,
               ipea_shift=None, imaginary_shift=None, level_shift=None,
               runtype=None, basis=None, reference="rhf", **keywords):
        """Use a compact OpenQP CASPT2 setup.

        `variant` selects `caspt2` (single state, the default), `ms-caspt2`
        or `xms-caspt2`; it is the input `method`, not a [pt2] key."""
        self._require_active_space("CASPT2", active_electrons, active_orbitals)
        method = str(variant or "caspt2").lower().replace("_", "-")
        if method in {"ms", "multistate"}:
            method = "ms-caspt2"
        elif method in {"xms", "extended-multistate"}:
            method = "xms-caspt2"
        if method not in {"caspt2", "ms-caspt2", "xms-caspt2"}:
            raise ValueError(
                "CASPT2 variant must be 'caspt2', 'ms-caspt2' or 'xms-caspt2'."
            )
        opts = dict(keywords.pop("pt2", None) or {})
        for key, value in (("h0", h0), ("ipea_shift", ipea_shift),
                           ("imaginary_shift", imaginary_shift),
                           ("level_shift", level_shift)):
            if value is not None:
                opts[key] = value
        return self._wf_setup(
            method, runtype=runtype, basis=basis, reference=reference,
            active_electrons=active_electrons, active_orbitals=active_orbitals,
            frozen_core=frozen_core, nroot=nroot, pt2=opts or None, **keywords)

    def nevpt2(self, active_electrons=None, active_orbitals=None,
               frozen_core=None, nroot=1, contraction=None, runtype=None,
               basis=None, reference="rhf", **keywords):
        """Use a compact OpenQP NEVPT2 setup.

        NEVPT2 is CASPT2's determinant machinery with the Dyall H0, so it is
        `method=caspt2` plus `[pt2] h0=dyall`.  `contraction='strong'` gives
        SC-NEVPT2; the default is the uncontracted form."""
        self._require_active_space("NEVPT2", active_electrons, active_orbitals)
        opts = dict(keywords.pop("pt2", None) or {})
        opts.setdefault("h0", "dyall")
        if contraction is not None:
            opts["contraction"] = contraction
        return self._wf_setup(
            "caspt2", runtype=runtype, basis=basis, reference=reference,
            active_electrons=active_electrons, active_orbitals=active_orbitals,
            frozen_core=frozen_core, nroot=nroot, pt2=opts, **keywords)

    def qdpt2(self, active_electrons=None, active_orbitals=None,
              frozen_core=None, nroot=1, variant=None, edshft=None,
              runtype=None, basis=None, reference="rhf", **keywords):
        """Use a compact OpenQP QDPT2 setup in the GAMESS convention.

        `variant` selects `mrmp2` (single state, the default), `mcqdpt2`
        (multistate) or `xmcqdpt2` (Granovsky-extended)."""
        self._require_active_space("QDPT2", active_electrons, active_orbitals)
        method = str(variant or "mrmp2").lower().replace("-", "").replace("_", "")
        if method not in {"mrmp2", "mcqdpt2", "xmcqdpt2"}:
            raise ValueError(
                "QDPT2 variant must be 'mrmp2', 'mcqdpt2' or 'xmcqdpt2'."
            )
        opts = dict(keywords.pop("pt2", None) or {})
        if edshft is not None:
            opts["edshft"] = edshft
        return self._wf_setup(
            method, runtype=runtype, basis=basis, reference=reference,
            active_electrons=active_electrons, active_orbitals=active_orbitals,
            frozen_core=frozen_core, nroot=nroot, pt2=opts or None, **keywords)

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

    @property
    def dftb(self):
        """Callable [dftb] section proxy.

        ``job.dftb(...)`` runs the DFTB workflow helper below, while
        ``job.dftb.option`` / ``job.dftb.option = value`` read and write the
        [dftb] section like every other schema section.
        """
        return _DFTBSectionProxy(self, "dftb")

    def tddftb(self, **kwargs):
        """Use the conventional singlet TD-DFTB (TDA) response helper."""
        return self._dftb(response_type="tddftb", **kwargs)

    def ground_dftb(self, **kwargs):
        """Use the ground-state SCC-DFTB helper explicitly."""
        return self._dftb(response_type="ground", **kwargs)

    def sf_tddftb(self, **kwargs):
        """Use the spin-flip TD-DFTB response helper."""
        return self._dftb(response_type="sf", **kwargs)

    def mrsf_tddftb(self, **kwargs):
        """Use the mixed-reference spin-flip TD-DFTB response helper."""
        return self._dftb(response_type="mrsf", **kwargs)

    def _dftb(self, runtype=None, response_type="mrsf", nstate=3,
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
        # Resolve the response type BEFORE draining schema keywords: `type` is a
        # [dftb] schema key, so the generic drain below would otherwise consume
        # an explicit job.dftb(type=...) and silently fall back to response_type.
        requested_type = canonical_dftb_type(keywords.pop("type", response_type))
        for key in list(keywords.keys()):
            if key in dftb_schema and key != "type":
                dftb_updates[key] = keywords.pop(key)

        tdhf_type = {
            "ground": "tda",
            "ground_noscc": "tda",
            "tddftb": "tda",
            "sf": "sf",
            "mrsf": "mrsf",
        }.get(requested_type)
        if tdhf_type is None:
            raise ValueError("DFTB response_type must be ground, tddftb, sf, or mrsf.")

        dftb_updates["type"] = requested_type
        self.section("dftb", **dftb_updates)
        if requested_type in {"sf", "mrsf"}:
            self.scf(type="rohf", multiplicity=3)
        return self.tdhf(type=tdhf_type, nstate=nstate, **keywords)

    @property
    def xtb(self):
        """Callable [xtb] section proxy.

        ``job.xtb(...)`` runs the xTB workflow helper below, while
        ``job.xtb.option`` / ``job.xtb.option = value`` read and write the
        [xtb] section like every other schema section.
        """
        return _XTBSectionProxy(self, "xtb")

    def _xtb(self, runtype=None, response_type="mrsf", nstate=3,
             parameter_path=None, **keywords):
        """Use the optional OpenQP-xTB backend through the normal OpenQP workflow."""
        input_updates = {"method": "xtb", "functional": ""}
        if runtype is not None:
            input_updates["runtype"] = runtype
        self.input(**input_updates)

        xtb_schema = OQP_CONFIG_SCHEMA.get("xtb", {})
        xtb_updates = {}
        if parameter_path is not None:
            xtb_updates["parameter_path"] = parameter_path
        # Resolve the response type BEFORE draining schema keywords: `type` is an
        # [xtb] schema key, so the generic drain below would otherwise consume
        # an explicit job.xtb(type=...) and silently fall back to response_type.
        requested_type = str(keywords.pop("type", response_type)).lower()
        for key in list(keywords.keys()):
            if key in xtb_schema and key != "type":
                xtb_updates[key] = keywords.pop(key)

        xtb_type = requested_type
        tdhf_type = {
            "ground": "tda",
            "dftb": "tda",
            "dftb0": "tda",
            "noscc": "tda",
            "ground_noscc": "tda",
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
            raise ValueError("xTB response_type must be ground, tddftb, sf, or mrsf.")

        # The backend's response-method names are ground/tddftb/sf/mrsf; map the
        # tda/td-dftb aliases onto tddftb so the explicit request and the auto
        # path (tdhf.type=tda/rpa -> tddftb) reach the same backend method.
        _BACKEND_TYPE = {"tda": "tddftb", "td-dftb": "tddftb"}
        xtb_updates["type"] = _BACKEND_TYPE.get(xtb_type, xtb_type)
        self.section("xtb", **xtb_updates)
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
        # TB backends ([dftb]/[xtb] section named after the method): the MRSF
        # response type may be explicit or derived (type=auto).
        tb_type = str(self.config_typed.get(method, {}).get("type", "auto")).lower()
        tb_mrsf = is_tb_method(method) and tb_type in {
            "auto", "mrsf", "mrsftddftb", "mrsf-tddftb"} and response == "mrsf"
        if not ((method == "tdhf" and response == "mrsf") or tb_mrsf):
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
        active_runtype = str(
            self.config_typed.get("input", {}).get("runtype", "energy")
        ).lower()
        public_aliases = {
            # the crossing driver decides which search key algorithm names
            "algorithm": (
                "mecp_search" if active_runtype == "mecp" else "meci_search"
            ),
            "sigma": "pen_sigma",
            "alpha": "pen_alpha",
            "delta_beta": "pen_delta",
            "beta_schedule": "pen_jump",
            "gap": "energy_gap",
        }
        for public, internal in public_aliases.items():
            if public not in requested:
                continue
            if internal in requested:
                raise ValueError(
                    f"Optimizer option '{internal}' is specified twice through "
                    f"'{public}' and '{internal}'."
                )
            requested[internal] = requested.pop(public)
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
            if option in {"states", "pen_jump"} and isinstance(value, (list, tuple)):
                value = ",".join(str(item) for item in value)
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
