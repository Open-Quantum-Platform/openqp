from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
from oqp.utils.input_parser import OQPConfigParser  # your class
from oqp.pyoqp import Runner                        # your Runner
from oqp.utils.kword_map import resolve_param_key   # your pure dotted resolver
import os
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
        # list/tuple of atoms
        if isinstance(system, (list, tuple)):
            rows = []
            for item in system:
                if isinstance(item, (list, tuple)) and len(item) >= 4:
                    el, x, y, z = item[:4]
                    rows.append(f"{el} {x} {y} {z}")
                else:
                    raise ValueError("Each atom must be (symbol, x, y, z).")
            lines = rows

        elif isinstance(system, str):
            s = system.strip()

            if os.path.exists(s):
                return s

            if "\n" in s:
                lines = [ln.strip() for ln in s.split("\n") if ln.strip()]
            elif ";" in s or "," in s:
                sep = ";" if ";" in s else ","
                parts = [p.strip() for p in s.split(sep) if p.strip()]
                if not all(len(p.split()) >= 4 for p in parts):
                    raise ValueError("Inline atom must be 'El x y z' per entry.")
                lines = parts
            else:
                toks = s.split()
                if len(toks) >= 4:
                    lines = [s]
                else:
                    return s
        else:
            raise TypeError("input.system must be str or list/tuple of atoms.")

        return "\n" + "\n".join(lines)

    @staticmethod
    def _dump_strings_from_parser(parser: OQPConfigParser):
        """Extract a pure string dict from parser (preserves comma formats for arrays)."""
        out = {}
        for sec in parser.sections():
            out[sec] = {}
            for opt, val in parser[sec].items():
                out[sec][opt] = val  # already a string
        return out

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
