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


class _SectionProxy:
    """Attribute/call proxy for a section in the OpenQP input schema."""

    def __init__(self, owner, section):
        object.__setattr__(self, "_owner", owner)
        object.__setattr__(self, "_section", section)

    def __call__(self, **kwargs):
        self._owner.section(self._section, **kwargs)
        return self._owner

    def __getattr__(self, option):
        schema = OQP_CONFIG_SCHEMA.get(self._section, {})
        if option not in schema:
            raise AttributeError(f"Unknown OpenQP option '{self._section}.{option}'.")
        return self._owner.config_typed.get(self._section, {}).get(option)

    def __setattr__(self, option, value):
        if option.startswith("_"):
            object.__setattr__(self, option, value)
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

    def molecule(self, system=None, basis=None, charge=None, unit="Angstrom",
                 geometry=None, source="auto", timeout=10, **kwargs):
        """Set molecular system data using OpenQP input-section keywords."""
        if system is not None and geometry is not None:
            raise ValueError("Use either system or geometry, not both.")
        if geometry is not None:
            system = get_geometry(geometry, source=source, timeout=timeout)
        if system is None:
            raise ValueError("molecule requires a system or geometry value.")

        self.unit = unit
        updates = {"input.system": system}
        if basis is not None:
            updates["input.basis"] = basis
        if charge is not None:
            updates["input.charge"] = charge
        for option, value in kwargs.items():
            updates[f"input.{option}"] = value
        return self.set(**updates)

    def hf(self, reference="rhf", runtype="energy", multiplicity=None, **scf_keywords):
        """Use a compact OpenQP HF setup for ordinary single-reference jobs."""
        self.input(method="hf", runtype=runtype)
        updates = {}
        if reference is not None:
            updates["type"] = reference
        if multiplicity is not None:
            updates["multiplicity"] = multiplicity
        updates.update(scf_keywords)
        if updates:
            self.scf(**updates)
        return self

    def dft(self, functional, reference="rhf", runtype="energy",
            multiplicity=None, **scf_keywords):
        """Use a compact OpenQP DFT setup for ordinary Kohn-Sham jobs."""
        self.input(method="hf", functional=functional, runtype=runtype)
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
             runtype="energy", **tdhf_keywords):
        """Use a compact OpenQP MRSF-TDDFT setup with an open-shell reference."""
        self.input(method="tdhf", runtype=runtype)
        self.scf(type=reference, multiplicity=multiplicity)
        updates = {"type": "mrsf", "nstate": nstate}
        updates.update(tdhf_keywords)
        return self.tdhf(**updates)

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
