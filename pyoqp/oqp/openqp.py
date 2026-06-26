"""High-level Python helpers for OpenQP calculations."""

from copy import deepcopy
import os
from collections.abc import Mapping

from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
from oqp.pyoqp import Runner
from oqp.utils.constants import ANGSTROM_TO_BOHR
from oqp.utils.input_parser import OQPConfigParser
from oqp.utils.kword_map import resolve_param_key


_BOHR_UNITS = {"bohr", "au", "a.u.", "atomic_unit", "atomic_units"}
_INPUT_ALIASES = {
    "atom": "input.system",
    "system": "input.system",
    "basis": "input.basis",
    "charge": "input.charge",
    "functional": "input.functional",
    "method": "input.method",
    "run_type": "input.runtype",
    "runtype": "input.runtype",
    "ispher": "input.ispher",
    "omp_threads": "input.omp_threads",
    "scf_type": "scf.type",
    "reference": "scf.type",
    "multiplicity": "scf.multiplicity",
    "td_type": "tdhf.type",
    "tdhf_type": "tdhf.type",
    "nstate": "tdhf.nstate",
}
_SECTION_ALIASES = {
    "input.atom": "input.system",
    "input.run_type": "input.runtype",
    "tdhf.states": "tdhf.nstate",
}


def _is_bohr(unit):
    return str(unit).strip().lower() in _BOHR_UNITS


def _format_coord(value, unit):
    if not _is_bohr(unit):
        return str(value)
    coord = float(str(value).replace("D", "E").replace("d", "e"))
    return f"{coord * ANGSTROM_TO_BOHR:.12g}"


def _is_coord_sequence(value):
    if isinstance(value, (str, bytes)):
        return False
    try:
        iter(value)
        len(value)
    except TypeError:
        return False
    return True


def _normalize_atom_row(item, unit):
    if isinstance(item, str):
        return _normalize_system_line(item, unit)

    if not _is_coord_sequence(item) or len(item) < 2:
        raise ValueError("Each atom must be 'El x y z', (symbol, x, y, z), or (symbol, (x, y, z)).")

    element = item[0]
    if len(item) >= 4:
        coords = item[1:4]
        extra = item[4:]
    elif _is_coord_sequence(item[1]) and len(item[1]) >= 3:
        coords = item[1][:3]
        extra = item[2:]
    else:
        raise ValueError("Each atom must be 'El x y z', (symbol, x, y, z), or (symbol, (x, y, z)).")

    fields = [str(element), *(_format_coord(coord, unit) for coord in coords), *(str(x) for x in extra)]
    return " ".join(fields)


def _normalize_system_line(line, unit):
    tokens = line.split()
    if len(tokens) < 4:
        return line.strip()
    fields = [tokens[0], *(_format_coord(coord, unit) for coord in tokens[1:4]), *tokens[4:]]
    return " ".join(fields)


def normalize_system(system, unit="Angstrom"):
    """
    Normalize inline molecular geometries for OpenQP input.system.

    Coordinates are written in Angstrom because the OpenQP input reader converts
    the text geometry into its internal Bohr representation. PySCF-style Bohr
    inputs can be passed with unit="Bohr".
    """
    if isinstance(system, (list, tuple)):
        lines = [_normalize_atom_row(item, unit) for item in system]
    elif isinstance(system, str):
        s = system.strip()

        if os.path.exists(s):
            return s

        if "\n" in s:
            lines = [_normalize_system_line(ln, unit) for ln in s.split("\n") if ln.strip()]
        elif ";" in s or "," in s:
            sep = ";" if ";" in s else ","
            parts = [p.strip() for p in s.split(sep) if p.strip()]
            if not all(len(p.split()) >= 4 for p in parts):
                raise ValueError("Inline atom must be 'El x y z' per entry.")
            lines = [_normalize_system_line(part, unit) for part in parts]
        else:
            toks = s.split()
            if len(toks) >= 4:
                lines = [_normalize_system_line(s, unit)]
            else:
                return s
    else:
        raise TypeError("input.system must be str or list/tuple of atoms.")

    return "\n" + "\n".join(lines)


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
    Pythonic, PySCF-compatible convenience layer over the OpenQP input schema.

    This class is additive: it builds the same sectioned input dictionary used
    by Runner and existing OpenQP input files, while allowing common molecular
    arguments such as atom, basis, charge, and spin.
    """

    def __init__(
        self,
        atom=None,
        system=None,
        basis=None,
        charge=None,
        spin=None,
        multiplicity=None,
        unit="Angstrom",
        project="oqp_project",
        log=None,
        silent=0,
        usempi=True,
        **kwargs,
    ):
        if atom is not None and system is not None:
            raise ValueError("Use either atom= or system=, not both.")

        self.project = project or "oqp_project"
        self.log = log if log is not None else f"{self.project}.log"
        self.silent = silent
        self.usempi = usempi
        self.unit = unit
        self.runner = None
        self.mol = None

        parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA)
        self.config_str = dump_strings_from_parser(parser)
        self.config_typed = parser.validate()

        updates = {}
        molecule = system if system is not None else atom
        if molecule is not None:
            updates["input.system"] = molecule
        if basis is not None:
            updates["input.basis"] = basis
        if charge is not None:
            updates["input.charge"] = charge
        if spin is not None and multiplicity is None:
            multiplicity = int(spin) + 1
        if multiplicity is not None:
            updates["scf.multiplicity"] = multiplicity

        updates.update(kwargs)
        if updates:
            self.set(**updates)

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
        params = {}
        for attr in ("atom", "basis", "charge", "spin", "unit"):
            if hasattr(mol, attr):
                value = getattr(mol, attr)
                if value is not None:
                    params[attr] = value
        params.update(kwargs)
        return cls(**params)

    def section(self, section, **kwargs):
        """Update one OpenQP section using keyword arguments."""
        if section not in OQP_CONFIG_SCHEMA:
            raise KeyError(f"Unknown OpenQP section '{section}'.")
        return self.set(**{f"{section}.{opt}": value for opt, value in kwargs.items()})

    def set(self, **kwargs):
        """
        Update configuration with dotted OpenQP keys, section dictionaries, or
        convenience aliases such as basis=, atom=, spin=, and nstate=.
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
        if key == "spin":
            return "scf.multiplicity", int(value) + 1
        if key in _INPUT_ALIASES:
            key = _INPUT_ALIASES[key]
        elif key in _SECTION_ALIASES:
            key = _SECTION_ALIASES[key]

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
