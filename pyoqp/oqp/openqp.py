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

    @staticmethod
    def _looks_like_basis(value):
        if not isinstance(value, str):
            return False
        text = value.strip()
        if text.lower().endswith((".xyz", ".mol", ".sdf", ".pdb")):
            return False
        return bool(text and "\n" not in text and ";" not in text and " " not in text)

    def control(self, runtype=None, omp_threads=None, **kwargs):
        """Set run-level controls such as runtype, OpenMP, and optimization options."""
        updates = {}
        if runtype is not None:
            updates["input.runtype"] = runtype
        if omp_threads is not None:
            updates["input.omp_threads"] = omp_threads
        if updates:
            self.set(**updates)

        if kwargs:
            active_runtype = str(
                runtype if runtype is not None
                else self.config_typed.get("input", {}).get("runtype", "energy")
            ).lower()
            if active_runtype in OPTIMIZER_RUNTYPES:
                return self._set_optimize_options(**kwargs)
            raise KeyError(
                "Extra job.control(...) options are supported for optimization "
                "runtypes. Use section helpers for other workflow-specific keywords."
            )
        return self

    def theory(self, method, functional=None, basis=None, runtype=None,
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
            "Unknown theory method. Use hf, dft, or mrsf-tddft."
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

    def soc(self, nstate=3, functional=None, reference="rohf",
            reference_multiplicity=3, soc_2e=1, basis=None, **tdhf_keywords):
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
        self.scf(type=reference, multiplicity=reference_multiplicity)
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
    def __init__(self, cfg: dict):
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
            silent=0,
            usempi=True
        )
        self.mol = self.runner.mol

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
