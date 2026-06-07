"""OpenQP input checker with actionable diagnostics."""

from __future__ import annotations

from dataclasses import dataclass, field
import multiprocessing
import os
from typing import Any

from oqp.utils.mpi_utils import MPIManager


SUPPORTED_RUNTYPES = {
    "energy", "grad", "hess", "nac", "nacme", "bp", "optimize",
    "meci", "mecp", "tci", "mep", "ts", "irc", "neb", "prop", "data", "ekt",
}
NOT_AVAILABLE_RUNTYPES = {"soc", "md"}
ALL_RUNTYPES = SUPPORTED_RUNTYPES | NOT_AVAILABLE_RUNTYPES
METHODS = {"hf", "tdhf"}
SCF_TYPES = {"rhf", "rohf", "uhf"}
TDHF_TYPES = {"rpa", "tda", "sf", "mrsf", "umrsf", "mrsf_ekt_ip", "mrsf_ekt_ea"}
GUESS_TYPES = {"huckel", "modhuckel", "hcore", "json", "auto", "sap", "minao"}
SCF_CONVERGERS = {"diis", "soscf", "trah"}
OPTIONAL_SCF_CONVERGERS = SCF_CONVERGERS | {"none", ""}
DIIS_TYPES = {"none", "cdiis", "ediis", "adiis", "vdiis"}
OPT_LIBS = {"scipy", "geometric", "oqp"}
SCIPY_OPTIMIZERS = {"bfgs", "cg", "l-bfgs-b", "newton-cg"}
MECI_SEARCH = {"penalty", "ubp", "hybrid"}
SCF_PROPS = {"el_mom", "mulliken", "nmr"}
NMR_GAUGES = {"cgo", "giao"}
INIT_SCF_TYPES = {"no", "rhf", "uhf", "rohf", "rks", "uks", "roks"}

WIKI_HELP = {
    "input.runtype": "Use energy, ekt, grad, hess, nac, nacme, optimize, meci, mecp, mep, ts, irc, neb, prop, or data. soc and md are recognized but not implemented yet.",
    "input.method": "Use method=hf for HF/DFT or method=tdhf for TDHF/TDDFT/SF/MRSF.",
    "input.system": "Set system to an XYZ file path or inline coordinates with one atom per indented line.",
    "input.basis": "Set basis to a basis name, a comma-separated per-atom list, or library with tagged atoms and [input] library mappings.",
    "scf.type": "RHF is for multiplicity 1 closed-shell references. SF/MRSF needs an open-shell reference, usually ROHF.",
    "tdhf.type": "Use rpa or tda for ordinary TDHF/TDDFT, sf or mrsf for spin-flip, umrsf only with UHF, and legacy mrsf_ekt_ip/mrsf_ekt_ea only with energy runtype. EKT analysis must use [input] runtype=ekt with [tdhf] type=mrsf and [ekt] IP, EA, or both.",
    "tdhf.nstate": "nstate must cover the highest excited-state index requested anywhere else in the input.",
    "guess.type": "Use huckel or modhuckel (weighted Wolfsberg-Helmholz) for native extended-Huckel guesses, hcore for the bare core Hamiltonian, sap for the native superposition-of-atomic-potentials guess, minao for projected minimal-basis densities, json with a JSON restart file, or auto for JSON-if-present otherwise Huckel.",
    "optimize.lib": "oqp is the default optimizer backend: the built-in NumPy/SciPy optimizer (redundant internals/DLC/TRIC + restricted-step RFO) supporting optimize, ts, meci, mecp, tci, neb, irc, and mep. geometric (the external geomeTRIC package) supports optimize, MECI, MECP, TS, IRC, and NEB. scipy supports optimize, meci, mecp, and mep.",
    "nac.states": "Use state pairs such as 1 2,2 3 for NAC calculations. Each index must be a TDHF excited state.",
}


_FALSE_BOOL = {"false", "0", "f", ".false.", "off", "no"}
_TRUE_BOOL = {"true", "1", "t", ".true.", "on", "yes"}

@dataclass
class Diagnostic:
    severity: str
    path: str
    message: str
    value: Any = None
    expected: str | None = None
    action: str | None = None
    wiki: str | None = None

    def to_text(self, index: int) -> str:
        lines = [f"{index}. [{self.severity}] {self.path}: {self.message}"]
        if self.value is not None:
            lines.append(f"   current: {self.value!r}")
        if self.expected:
            lines.append(f"   expected: {self.expected}")
        if self.action:
            lines.append(f"   fix: {self.action}")
        if self.wiki:
            lines.append(f"   note: {self.wiki}")
        return "\n".join(lines)


@dataclass
class CheckReport:
    diagnostics: list[Diagnostic] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not any(item.severity == "ERROR" for item in self.diagnostics)

    @property
    def errors(self) -> list[Diagnostic]:
        return [item for item in self.diagnostics if item.severity == "ERROR"]

    @property
    def warnings(self) -> list[Diagnostic]:
        return [item for item in self.diagnostics if item.severity == "WARNING"]

    def add(
        self,
        severity: str,
        path: str,
        message: str,
        *,
        value: Any = None,
        expected: str | None = None,
        action: str | None = None,
        wiki: str | None = None,
    ) -> None:
        self.diagnostics.append(
            Diagnostic(
                severity=severity,
                path=path,
                message=message,
                value=value,
                expected=expected,
                action=action,
                wiki=wiki,
            )
        )

    def to_text(self) -> str:
        title = "OpenQP input check: PASS" if self.ok else "OpenQP input check: FAIL"
        if not self.diagnostics:
            return f"{title}\nNo problems found."

        counts = (
            f"{len(self.errors)} error(s), {len(self.warnings)} warning(s), "
            f"{len(self.diagnostics) - len(self.errors) - len(self.warnings)} info message(s)"
        )
        body = "\n\n".join(item.to_text(index + 1) for index, item in enumerate(self.diagnostics))
        return f"{title}\n{counts}\n\n{body}"


def _get(config: dict[str, Any], section: str, option: str, default: Any = None) -> Any:
    return config.get(section, {}).get(option, default)


def _as_lower(value: Any) -> Any:
    return value.lower() if isinstance(value, str) else value


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, str):
        if not value.strip():
            return []
        return [item.strip() for item in value.split(",") if item.strip()]
    return [value]


def _norm_path(raw_path: str, system_text: str) -> str:
    if os.path.isabs(raw_path):
        return raw_path
    parent = os.getcwd()
    first_line = system_text.split("\n", 1)[0].strip()
    if first_line and os.path.splitext(first_line)[1].lower() == ".xyz":
        parent = os.path.dirname(os.path.abspath(first_line))
    return os.path.abspath(os.path.join(parent, raw_path))



def _parse_bool_like(value: Any, path: str, report: CheckReport, *, allow_true: bool = True, allow_auto: bool = False, strict_false_only: bool = False, experimental_warning=None) -> bool | str | None:
    """Validate boolean-like values used by Python-only symmetry gates."""
    if isinstance(value, bool):
        if strict_false_only and value:
            report.add(
                "ERROR",
                path,
                "Symmetry reduction flags are intentionally disabled by default.",
                value=value,
                expected="False",
                action="Keep symmetry reduction off until dedicated kernels are production-ready.",
            )
            return None
        return value

    if isinstance(value, (int, float)):
        if value in (0, 1):
            parsed = bool(int(value))
            if strict_false_only and parsed:
                report.add(
                    "ERROR",
                    path,
                    "Symmetry reduction flags are intentionally disabled by default.",
                    value=value,
                    expected="False",
                    action="Use 0 or false for disabled symmetry reductions.",
                )
                return None
            return parsed
        report.add(
            "ERROR",
            path,
            "Boolean option expects false-like or true-like values.",
            value=value,
            expected="False-like or True-like value",
            action="Use false/true, 0/1, or accepted string tokens.",
        )
        return None

    if not isinstance(value, str):
        report.add(
            "ERROR",
            path,
            "Boolean option expects a boolean-like value.",
            value=value,
            expected="False-like or True-like value",
            action="Use false/true, 0/1, or accepted string tokens.",
        )
        return None

    lowered = value.strip().lower()
    if lowered in _FALSE_BOOL:
        return False
    if lowered == "full" and experimental_warning:
        report.add(
            "WARNING", path,
            experimental_warning + " ('full' additionally enables the "
            "non-abelian full point group, ~1e-7 accuracy)",
            value=value, expected="False (default)",
            action="Validate results against a C1 reference run.",
        )
        return True
    if lowered in _TRUE_BOOL:
        if strict_false_only:
            report.add(
                "ERROR",
                path,
                "Symmetry reduction flags are intentionally disabled by default.",
                value=value,
                expected="False",
                action="Keep symmetry reductions off until validation and production kernels are in place.",
            )
            return None
        if experimental_warning:
            report.add(
                "WARNING",
                path,
                experimental_warning,
                value=value,
                expected="False (default)",
                action="Validate results against a C1 reference run.",
            )
        return True if allow_true else False
    if allow_auto and lowered == "auto":
        return "auto"

    report.add(
        "ERROR",
        path,
        "Boolean option expects false/true-like values.",
        value=value,
        expected="false/true, 0/1, or auto",
        action="Use recognized tokens like false, true, 0, 1, yes, no, or auto where supported.",
    )
    return None


def _check_symmetry(config: dict[str, Any], report: CheckReport) -> None:
    """Validate the new symmetry metadata block without enabling reductions by default."""
    section = config.get("symmetry")
    if not section:
        return
    if not isinstance(section, dict):
        report.add(
            "ERROR",
            "symmetry",
            "[symmetry] must be a mapping.",
            value=section,
            expected="[symmetry] block",
            action="Pass symmetry options as a mapping.",
        )
        return

    _parse_bool_like(section.get("enabled", "false"), "symmetry.enabled", report, allow_auto=True)
    if "point_group" in section and not isinstance(section.get("point_group"), str):
        report.add(
            "ERROR",
            "symmetry.point_group",
            "point_group must be a string label.",
            value=section.get("point_group"),
            expected="auto or point-group label",
            action="Use auto or a valid point-group abbreviation.",
        )
    if "subgroup" in section and not isinstance(section.get("subgroup"), str):
        report.add(
            "ERROR",
            "symmetry.subgroup",
            "subgroup must be a string label.",
            value=section.get("subgroup"),
            expected="auto or subgroup label",
            action="Use auto or a valid subgroup abbreviation.",
        )

    _parse_bool_like(section.get("label_mo", True), "symmetry.label_mo", report)
    _parse_bool_like(section.get("label_states", True), "symmetry.label_states", report)
    _parse_bool_like(section.get("label_modes", True), "symmetry.label_modes", report)
    _parse_bool_like(
        section.get("use_integral_symmetry", "False"),
        "symmetry.use_integral_symmetry",
        report,
        allow_true=True,
        experimental_warning=(
            "Experimental: petite-list/skeleton-Fock reduction. The molecule "
            "is reoriented to the symmetry standard orientation at load time."
        ),
    )
    _parse_bool_like(
        section.get("use_response_symmetry", "False"),
        "symmetry.use_response_symmetry",
        report,
        allow_true=True,
        experimental_warning=(
            "Experimental: irrep-blocked Davidson updates for the response "
            "solver. Validate excitation energies against an unblocked run."
        ),
    )
    _parse_bool_like(section.get("strict", False), "symmetry.strict", report)

    try:
        tolerance = float(section.get("tolerance", 1.0e-5))
    except (TypeError, ValueError):
        report.add(
            "ERROR",
            "symmetry.tolerance",
            "symmetry.tolerance must be numeric.",
            value=section.get("tolerance"),
            expected="positive float",
            action="Set tolerance to a positive float such as 1.0e-5.",
        )
    else:
        if tolerance <= 0.0:
            report.add(
                "ERROR",
                "symmetry.tolerance",
                "symmetry.tolerance must be positive.",
                value=tolerance,
                expected="> 0.0",
                action="Use positive tolerance with stricter or looser default (e.g., 1.0e-5).",
            )


def _iter_coordinate_lines(system: str) -> tuple[list[str], str | None]:
    lines = system.split("\n")
    if not lines:
        return [], None

    first_line = lines[0].strip()
    if first_line:
        return [], first_line
    return [line for line in lines[1:] if line.strip()], None


def _max_state(values: list[Any]) -> int:
    state_max = 0
    for value in values:
        if isinstance(value, (list, tuple)):
            ints = [int(x) for x in value]
            if ints:
                state_max = max(state_max, max(ints))
        elif value not in (None, ""):
            state_max = max(state_max, int(value))
    return state_max


def _omp_threads(report: CheckReport) -> int:
    raw = os.environ.get("OMP_NUM_THREADS", "1")
    try:
        value = int(raw)
        if value < 1:
            raise ValueError
        return value
    except (TypeError, ValueError):
        report.add(
            "WARNING",
            "env.OMP_NUM_THREADS",
            "OMP_NUM_THREADS is missing or invalid; assuming 1 thread.",
            value=raw,
            expected="positive integer",
            action="Set OMP_NUM_THREADS=1 or another positive integer before running OpenQP.",
        )
        return 1


def _check_system(config: dict[str, Any], report: CheckReport) -> None:
    system = _get(config, "input", "system", "")
    if not system:
        report.add(
            "ERROR",
            "input.system",
            "Molecular geometry is missing.",
            expected="XYZ filename or inline coordinates",
            action="Set [input] system to a .xyz path or multiline coordinates.",
            wiki=WIKI_HELP["input.system"],
        )
        return

    inline_lines, xyz_path = _iter_coordinate_lines(system)
    if xyz_path:
        resolved = os.path.abspath(xyz_path)
        if not os.path.exists(resolved):
            report.add(
                "ERROR",
                "input.system",
                "Referenced XYZ file does not exist.",
                value=xyz_path,
                action="Fix the path or place the XYZ file in the working directory.",
                wiki=WIKI_HELP["input.system"],
            )
        return

    if not inline_lines:
        report.add(
            "ERROR",
            "input.system",
            "Inline geometry block is empty.",
            action="Add one atom per line under system= with element/x/y/z columns.",
            wiki=WIKI_HELP["input.system"],
        )
        return

    for index, line in enumerate(inline_lines, start=1):
        parts = line.split()
        if len(parts) < 4:
            report.add(
                "ERROR",
                "input.system",
                f"Atom line {index} is incomplete.",
                value=line,
                expected="symbol x y z",
                action="Provide at least four columns for each atom.",
                wiki=WIKI_HELP["input.system"],
            )


def _check_basis(config: dict[str, Any], report: CheckReport) -> None:
    basis = _get(config, "input", "basis", "")
    system = _get(config, "input", "system", "")
    library = _get(config, "input", "library", "")

    if not basis:
        report.add(
            "ERROR",
            "input.basis",
            "Basis set is missing.",
            action="Set [input] basis to a supported basis name.",
            wiki=WIKI_HELP["input.basis"],
        )
        return

    if basis != "library":
        return

    inline_lines, xyz_path = _iter_coordinate_lines(system)
    lines = inline_lines
    if xyz_path and os.path.exists(os.path.abspath(xyz_path)):
        with open(os.path.abspath(xyz_path), "r", encoding="utf-8") as handle:
            xyz_lines = handle.read().splitlines()
        try:
            num_atoms = int(xyz_lines[0])
            lines = xyz_lines[2:2 + num_atoms]
        except (IndexError, ValueError):
            lines = xyz_lines

    if not library.strip():
        report.add(
            "ERROR",
            "input.library",
            "Tagged basis mode requires [input] library mappings.",
            action="Add lines such as tag basis-name to [input] library.",
            wiki=WIKI_HELP["input.basis"],
        )
        return

    missing_tags = []
    for index, line in enumerate(lines, start=1):
        if len(line.split()) < 5:
            missing_tags.append(index)

    if missing_tags:
        report.add(
            "ERROR",
            "input.system",
            "basis=library requires a tag column on every atom line.",
            value=missing_tags,
            expected="symbol x y z tag",
            action="Append a per-atom tag and define matching entries in [input] library.",
            wiki=WIKI_HELP["input.basis"],
        )


def _check_guess(config: dict[str, Any], report: CheckReport) -> None:
    guess_type = _as_lower(_get(config, "guess", "type", "sap"))
    guess_file = _get(config, "guess", "file", "")
    guess_file2 = _get(config, "guess", "file2", "")
    continue_geom = _get(config, "guess", "continue_geom", False)
    swapmo = _get(config, "guess", "swapmo", "")

    if guess_type not in GUESS_TYPES:
        report.add(
            "ERROR",
            "guess.type",
            "Unknown guess type.",
            value=guess_type,
            expected=", ".join(sorted(GUESS_TYPES)),
            action="Choose one of the supported guess modes.",
            wiki=WIKI_HELP["guess.type"],
        )

    if guess_type == "json":
        if not guess_file:
            report.add(
                "ERROR",
                "guess.file",
                "guess.type=json requires a JSON restart file.",
                action="Set [guess] file to an existing .json restart file.",
                wiki=WIKI_HELP["guess.type"],
            )
        else:
            resolved = _norm_path(guess_file, _get(config, "input", "system", ""))
            if not os.path.exists(resolved):
                report.add(
                    "ERROR",
                    "guess.file",
                    "JSON restart file does not exist.",
                    value=guess_file,
                    action="Fix the file path or generate the JSON file first.",
                    wiki=WIKI_HELP["guess.type"],
                )

    if continue_geom and guess_type != "json":
        report.add(
            "WARNING",
            "guess.continue_geom",
            "continue_geom only takes effect with guess.type=json.",
            value=continue_geom,
            action="Use guess.type=json or disable continue_geom.",
        )

    if guess_file2:
        resolved = _norm_path(guess_file2, _get(config, "input", "system", ""))
        if not os.path.exists(resolved):
            report.add(
                "WARNING",
                "guess.file2",
                "Second restart file is set but does not exist.",
                value=guess_file2,
                action="Fix guess.file2 or remove it if NACME should use system2 instead.",
            )

    if swapmo:
        try:
            values = [int(item.strip()) for item in swapmo.split(",") if item.strip()]
        except ValueError:
            report.add(
                "ERROR",
                "guess.swapmo",
                "swapmo must contain comma-separated orbital indices.",
                value=swapmo,
                action="Use pairs like 12,13,20,21.",
            )
        else:
            if len(values) % 2 != 0:
                report.add(
                    "ERROR",
                    "guess.swapmo",
                    "swapmo requires an even number of orbital indices.",
                    value=swapmo,
                    expected="index pairs",
                    action="Provide pairs like i,j,k,l.",
                )


def _check_scf(config: dict[str, Any], report: CheckReport) -> None:
    scf_type = _as_lower(_get(config, "scf", "type", "rhf"))
    multiplicity = _get(config, "scf", "multiplicity", 1)
    converger = _as_lower(_get(config, "scf", "converger_type", "diis"))
    init_converger = _as_lower(_get(config, "scf", "init_converger", "diis"))
    diis_type = _as_lower(_get(config, "scf", "diis_type", "cdiis"))
    alternative_scf = _as_lower(_get(config, "scf", "alternative_scf", "trah"))
    init_scf = _as_lower(_get(config, "scf", "init_scf", "no"))
    functional = _get(config, "input", "functional", "")

    if scf_type not in SCF_TYPES:
        report.add(
            "ERROR",
            "scf.type",
            "Unknown SCF reference type.",
            value=scf_type,
            expected=", ".join(sorted(SCF_TYPES)),
            action="Choose rhf, rohf, or uhf.",
            wiki=WIKI_HELP["scf.type"],
        )

    if multiplicity < 1:
        report.add(
            "ERROR",
            "scf.multiplicity",
            "Multiplicity must be at least 1.",
            value=multiplicity,
            action="Use a positive spin multiplicity.",
            wiki=WIKI_HELP["scf.type"],
        )

    if multiplicity > 1 and scf_type == "rhf":
        report.add(
            "ERROR",
            "scf.type",
            "RHF is inconsistent with multiplicity > 1.",
            value=f"{scf_type}/{multiplicity}",
            expected="rohf or uhf for open-shell references",
            action="Switch scf.type to rohf or uhf.",
            wiki=WIKI_HELP["scf.type"],
        )

    if converger not in SCF_CONVERGERS:
        report.add(
            "ERROR",
            "scf.converger_type",
            "Unknown SCF converger.",
            value=converger,
            expected=", ".join(sorted(SCF_CONVERGERS)),
            action="Choose diis, soscf, or trah.",
        )

    if init_converger not in OPTIONAL_SCF_CONVERGERS:
        report.add(
            "ERROR",
            "scf.init_converger",
            "Unknown initial SCF converger.",
            value=init_converger,
            expected="none, empty, or one of " + ", ".join(sorted(SCF_CONVERGERS)),
            action="Use none to leave the initial converger unset, or choose diis, soscf, or trah.",
        )

    if diis_type not in DIIS_TYPES:
        report.add(
            "ERROR",
            "scf.diis_type",
            "Unknown DIIS mode.",
            value=diis_type,
            expected=", ".join(sorted(DIIS_TYPES)),
            action="Choose one of the implemented DIIS types.",
        )

    if alternative_scf not in SCF_CONVERGERS:
        report.add(
            "ERROR",
            "scf.alternative_scf",
            "Unknown fallback SCF converger.",
            value=alternative_scf,
            expected=", ".join(sorted(SCF_CONVERGERS)),
            action="Choose diis, soscf, or trah.",
        )

    if init_scf not in INIT_SCF_TYPES:
        report.add(
            "ERROR",
            "scf.init_scf",
            "Unknown initial SCF mode.",
            value=init_scf,
            expected=", ".join(sorted(INIT_SCF_TYPES)),
            action="Use no, rhf, uhf, rohf, rks, uks, or roks.",
        )

    if init_scf in {"rks", "uks", "roks"} and not functional:
        report.add(
            "ERROR",
            "scf.init_scf",
            "KS-style initial SCF requires a DFT functional.",
            value=init_scf,
            action="Set [input] functional or use rhf/uhf/rohf/no for init_scf.",
        )


def _check_tdhf(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    if method != "tdhf":
        return

    runtype = _as_lower(_get(config, "input", "runtype", "energy"))
    scf_type = _as_lower(_get(config, "scf", "type", "rhf"))
    scf_mult = _get(config, "scf", "multiplicity", 1)
    td_type = _as_lower(_get(config, "tdhf", "type", "rpa"))
    td_mult = _get(config, "tdhf", "multiplicity", 1)
    nstate = _get(config, "tdhf", "nstate", 1)
    nvdav = _get(config, "tdhf", "nvdav", 50)

    if td_type not in TDHF_TYPES:
        report.add(
            "ERROR",
            "tdhf.type",
            "Unknown TDHF response type.",
            value=td_type,
            expected=", ".join(sorted(TDHF_TYPES)),
            action="Choose rpa, tda, sf, mrsf, umrsf, mrsf_ekt_ip, or mrsf_ekt_ea.",
            wiki=WIKI_HELP["tdhf.type"],
        )
        return

    if td_type in {"mrsf_ekt_ip", "mrsf_ekt_ea"} and runtype != "energy":
        report.add(
            "ERROR",
            "tdhf.type",
            "Legacy tdhf.type=mrsf_ekt_ip/mrsf_ekt_ea is energy-only. EKT analysis must use [input] runtype=ekt with [tdhf] type=mrsf.",
            value=f"{runtype}/{td_type}",
            expected="input.runtype=ekt with tdhf.type=mrsf, or input.runtype=energy for legacy direct EKT calls",
            action="Set [input] runtype=ekt and [tdhf] type=mrsf for MRSF-EKT IP/EA analysis.",
            wiki=WIKI_HELP["tdhf.type"],
        )

    if td_type in {"rpa", "tda"} and scf_mult != td_mult:
        report.add(
            "INFO",
            "tdhf.multiplicity",
            "Response multiplicity differs from the SCF reference multiplicity.",
            value=td_mult,
            action="This is valid for state-specific singlet/triplet targets; keep it if intentional.",
        )

    if td_type in {"sf", "mrsf"} and scf_mult == td_mult:
        report.add(
            "INFO",
            "tdhf.multiplicity",
            "Spin-flip response multiplicity matches the SCF reference multiplicity.",
            value=td_mult,
            action="This can be intentional; verify the target state labeling if results look unexpected.",
        )

    if td_type in {"sf", "mrsf"} and scf_type != "rohf":
        report.add(
            "ERROR",
            "scf.type",
            "SF/MRSF requires an ROHF reference in the current code path.",
            value=scf_type,
            expected="rohf",
            action="Set [scf] type=rohf.",
            wiki=WIKI_HELP["tdhf.type"],
        )

    if td_type == "umrsf" and scf_type != "uhf":
        report.add(
            "ERROR",
            "scf.type",
            "UMRSF-TDDFT requires a UHF reference.",
            value=scf_type,
            expected="uhf",
            action="Set [scf] type=uhf.",
            wiki=WIKI_HELP["tdhf.type"],
        )

    if nstate < 1:
        report.add(
            "ERROR",
            "tdhf.nstate",
            "nstate must be at least 1.",
            value=nstate,
            action="Set nstate to the number of excited states to compute.",
            wiki=WIKI_HELP["tdhf.nstate"],
        )

    if nvdav < nstate:
        report.add(
            "WARNING",
            "tdhf.nvdav",
            "Davidson subspace is smaller than the number of requested states.",
            value=nvdav,
            expected=f">= {nstate}",
            action="Increase nvdav to at least nstate.",
        )


def _check_properties(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    runtype = _as_lower(_get(config, "input", "runtype", "energy"))
    grad_states = _as_list(_get(config, "properties", "grad", []))
    td_prop = _get(config, "properties", "td_prop", False)
    scf_prop = [_as_lower(item) for item in _as_list(_get(config, "properties", "scf_prop", []))]
    nmr_gauge = _as_lower(_get(config, "properties", "nmr_gauge", "cgo"))
    scf_type = _as_lower(_get(config, "scf", "type", "rhf"))

    for prop in scf_prop:
        if prop not in SCF_PROPS:
            report.add(
                "WARNING",
                "properties.scf_prop",
                "Unknown SCF property keyword.",
                value=prop,
                expected=", ".join(sorted(SCF_PROPS)),
                action="Remove the keyword or confirm that downstream code can ignore it.",
            )

    if "nmr" in scf_prop:
        if nmr_gauge not in NMR_GAUGES:
            report.add(
                "ERROR",
                "properties.nmr_gauge",
                "Unknown NMR shielding gauge formulation.",
                value=nmr_gauge,
                expected=", ".join(sorted(NMR_GAUGES)),
                action="Use nmr_gauge=giao (gauge-origin independent; RHF/UHF/ROHF) "
                       "or nmr_gauge=cgo (common gauge origin; closed-shell RHF).",
            )
        elif nmr_gauge == "cgo" and scf_type in ("uhf", "rohf"):
            report.add(
                "ERROR",
                "properties.nmr_gauge",
                "CGO NMR shielding supports closed-shell RHF references only.",
                value=f"{nmr_gauge} with scf.type={scf_type}",
                action="Use nmr_gauge=giao for open-shell (UHF/ROHF) NMR shielding.",
            )
        functional = _as_lower(_get(config, "input", "functional", ""))
        # Pre-flight mirror of the Fortran guards (the runtime also aborts):
        # the NMR response implements global hybrids only. Name-based, so it
        # cannot be exhaustive; unknown functionals are caught at runtime.
        if functional.startswith(("cam-", "dtcam-", "lc-", "lrc-", "wb97", "hse")):
            report.add(
                "ERROR",
                "input.functional",
                "NMR shielding with range-separated (CAM/LC) functionals is not implemented.",
                value=functional,
                action="Use HF, a pure functional, or a global hybrid (e.g. pbe0, b3lyp, bhhlyp).",
            )
        elif functional.startswith(("m06", "m08", "m11", "mn12", "mn15", "tpss",
                                    "scan", "rscan", "r2scan", "b97m", "revm06")):
            report.add(
                "ERROR",
                "input.functional",
                "NMR shielding with meta-GGA (tau-dependent) functionals is not implemented.",
                value=functional,
                action="Use HF, an LDA/GGA functional, or a global hybrid GGA (e.g. pbe0, b3lyp).",
            )
    if td_prop:
        report.add(
            "WARNING",
            "properties.td_prop",
            "TD properties are marked in documentation as not available yet.",
            value=td_prop,
            action="Disable td_prop unless you have a confirmed downstream implementation.",
        )

    if runtype != "grad":
        return

    if any(int(state) < 0 for state in grad_states if state not in (None, "")):
        report.add(
            "ERROR",
            "properties.grad",
            "Gradient state indices must be non-negative.",
            value=grad_states,
            action="Use state 0 for HF/DFT reference gradients or positive excited-state indices for TDHF.",
            wiki=WIKI_HELP["nac.states"],
        )

    if method == "hf" and any(int(state) > 0 for state in grad_states if state not in (None, "")):
        report.add(
            "ERROR",
            "properties.grad",
            "HF/DFT gradients are only available for state 0.",
            value=grad_states,
            expected="0",
            action="Set grad=0 or switch to method=tdhf for excited-state gradients.",
        )

    if method == "tdhf" and 0 in [int(state) for state in grad_states if state not in (None, "")]:
        report.add(
            "ERROR",
            "properties.grad",
            "TDHF gradients do not support state 0 in the current implementation.",
            value=grad_states,
            expected="positive excited-state indices",
            action="Use grad=1,2,... for TDHF/MRSF gradients.",
        )



def _check_requested_states(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    runtype = _as_lower(_get(config, "input", "runtype", "energy"))
    if method != "tdhf":
        return

    requested = []
    if runtype == "grad":
        requested.extend(_as_list(_get(config, "properties", "grad", [])))
    if runtype in {"optimize", "mep", "ts"}:
        requested.append(_get(config, "optimize", "istate", 0))
    if runtype in {"meci", "mecp"}:
        requested.append(_get(config, "optimize", "istate", 0))
        requested.append(_get(config, "optimize", "jstate", 0))
    if runtype == "tci":
        requested.append(_get(config, "optimize", "istate", 0))
        requested.append(_get(config, "optimize", "jstate", 0))
        requested.append(_get(config, "optimize", "kstate", 0))
    if runtype == "hess":
        requested.append(_get(config, "hess", "state", 0))
    if runtype in {"nac", "bp", "nacme"}:
        requested.extend(_as_list(_get(config, "nac", "states", [])))

    max_state = _max_state(requested)
    nstate = _get(config, "tdhf", "nstate", 1)

    if max_state > nstate:
        report.add(
            "ERROR",
            "tdhf.nstate",
            "Requested state index exceeds the number of computed excited states.",
            value=nstate,
            expected=f">= {max_state}",
            action="Increase [tdhf] nstate or lower the requested state indices.",
            wiki=WIKI_HELP["tdhf.nstate"],
        )


def _check_runtype(config: dict[str, Any], report: CheckReport,
                   input_dir: str | None = None) -> None:
    runtype = _as_lower(_get(config, "input", "runtype", "energy"))
    method = _as_lower(_get(config, "input", "method", "hf"))

    if runtype not in ALL_RUNTYPES:
        report.add(
            "ERROR",
            "input.runtype",
            "Unknown runtype.",
            value=runtype,
            expected=", ".join(sorted(ALL_RUNTYPES)),
            action="Choose one of the implemented run types.",
            wiki=WIKI_HELP["input.runtype"],
        )
        return

    if runtype in NOT_AVAILABLE_RUNTYPES:
        report.add(
            "ERROR",
            "input.runtype",
            "This runtype is recognized but not implemented.",
            value=runtype,
            action="Choose another runtype or implement the missing workflow first.",
            wiki=WIKI_HELP["input.runtype"],
        )
        return

    if runtype == "ekt":
        td_type = _as_lower(_get(config, "tdhf", "type", "rpa"))
        ekt_ip = bool(_get(config, "ekt", "ip", False))
        ekt_ea = bool(_get(config, "ekt", "ea", False))
        if method != "tdhf" or td_type != "mrsf":
            report.add(
                "ERROR",
                "input.runtype",
                "EKT runtype only supports MRSF-TDDFT.",
                value=f"{method}/{td_type}",
                expected="input.method=tdhf and tdhf.type=mrsf",
                action="Set [input] method=tdhf and [tdhf] type=mrsf for MRSF-EKT analysis.",
                wiki=WIKI_HELP["tdhf.type"],
            )
        if not ekt_ip and not ekt_ea:
            report.add(
                "ERROR",
                "ekt",
                "EKT runtype must request IP, EA, or both.",
                value={"ip": ekt_ip, "ea": ekt_ea},
                expected="[ekt] ip=True and/or ea=True",
                action="Enable [ekt] ip, ea, or both for EKT analysis.",
                wiki=WIKI_HELP["tdhf.type"],
            )
        return

    # UMRSF-TDDFT only implements the energy path. Every other runtype
    # eventually drives a gradient, Hessian, or Z-vector (grad/prop/data,
    # hess/thermo, nac/nacme, optimize/meci/mecp/mep/ts/irc/neb), none of
    # which exist for UMRSF yet. Reject them here at the single choke point
    # so validation fails early instead of dying at runtime.
    td_type = _as_lower(_get(config, "tdhf", "type", "rpa"))
    if method == "tdhf" and td_type == "umrsf" and runtype != "energy":
        report.add(
            "ERROR",
            "tdhf.type",
            "UMRSF-TDDFT only supports runtype=energy; "
            "gradients, Hessians, and Z-vectors are not implemented.",
            value=f"{td_type}/{runtype}",
            expected="energy",
            action="Use runtype=energy for UMRSF-TDDFT until UMRSF-TDDFT gradients/Z-vectors are implemented.",
            wiki=WIKI_HELP["tdhf.type"],
        )
        return

    if runtype == "grad":
        if method == "hf":
            if _max_state(_as_list(_get(config, "properties", "grad", []))) > 0:
                report.add(
                    "ERROR",
                    "properties.grad",
                    "HF/DFT grad runtype can only target state 0.",
                    value=_get(config, "properties", "grad", []),
                    action="Use grad=0 or switch to method=tdhf.",
                )
        elif method == "tdhf" and 0 in [int(state) for state in _as_list(_get(config, "properties", "grad", [])) if state not in (None, "")]:
            report.add(
                "ERROR",
                "properties.grad",
                "TDHF grad runtype cannot target state 0.",
                value=_get(config, "properties", "grad", []),
                action="Use grad=1,2,... for excited-state gradients.",
            )

    if runtype in {"optimize", "meci", "mecp", "mep", "ts", "irc", "neb"}:
        _check_optimize(config, report)

    if runtype == "neb":
        _check_neb(config, report, input_dir)

    if runtype in {"nac", "bp"}:
        _check_nac(config, report)

    if runtype == "nacme":
        _check_nacme(config, report)

    if runtype == "hess":
        _check_hess(config, report)


def _check_optimize(config: dict[str, Any], report: CheckReport) -> None:
    runtype = _as_lower(_get(config, "input", "runtype", "optimize"))
    method = _as_lower(_get(config, "input", "method", "hf"))
    lib = _as_lower(_get(config, "optimize", "lib", "oqp"))
    optimizer = _as_lower(_get(config, "optimize", "optimizer", "bfgs"))
    istate = _get(config, "optimize", "istate", 0)
    jstate = _get(config, "optimize", "jstate", 0)
    kstate = _get(config, "optimize", "kstate", 0)
    imult = _get(config, "optimize", "imult", 1)
    jmult = _get(config, "optimize", "jmult", 1)
    meci_search = _as_lower(_get(config, "optimize", "meci_search", "penalty"))

    if lib not in OPT_LIBS:
        report.add(
            "ERROR",
            "optimize.lib",
            "Unknown optimization library.",
            value=lib,
            expected=", ".join(sorted(OPT_LIBS)),
            action="Use scipy or geometric.",
            wiki=WIKI_HELP["optimize.lib"],
        )
        return

    if lib == "scipy" and optimizer not in SCIPY_OPTIMIZERS:
        report.add(
            "ERROR",
            "optimize.optimizer",
            "Unknown SciPy optimizer.",
            value=optimizer,
            expected=", ".join(sorted(SCIPY_OPTIMIZERS)),
            action="Use a supported SciPy optimizer.",
        )

    if method == "hf" and istate > 0:
        report.add(
            "ERROR",
            "optimize.istate",
            "HF/DFT optimization cannot target excited states.",
            value=istate,
            expected="0",
            action="Set istate=0 or switch to method=tdhf.",
        )

    if method == "tdhf" and runtype in {"optimize", "mep", "ts"} and istate == 0:
        report.add(
            "ERROR",
            "optimize.istate",
            "TDHF optimization cannot target state 0 in the current implementation.",
            value=istate,
            expected=">= 1",
            action="Set istate to a TDHF excited state index.",
        )

    if runtype == "meci":
        if method != "tdhf":
            report.add(
                "ERROR",
                "input.method",
                "MECI optimization requires method=tdhf.",
                value=method,
                expected="tdhf",
                action="Set [input] method=tdhf and configure [tdhf].",
            )
        if jstate <= istate:
            report.add(
                "ERROR",
                "optimize.jstate",
                "MECI requires jstate > istate.",
                value=f"{istate}/{jstate}",
                expected="jstate >= istate + 1",
                action="Set jstate to the next higher state.",
            )
        if meci_search not in MECI_SEARCH:
            report.add(
                "ERROR",
                "optimize.meci_search",
                "Unknown MECI search algorithm.",
                value=meci_search,
                expected=", ".join(sorted(MECI_SEARCH)),
                action="Use penalty, ubp, or hybrid.",
            )

    if runtype == "mecp" and imult == jmult:
        report.add(
            "ERROR",
            "optimize.jmult",
            "MECP requires different state multiplicities.",
            value=f"{imult}/{jmult}",
            expected="imult != jmult",
            action="Choose different multiplicities for the two crossing states.",
        )

    if lib == "scipy" and runtype in {"ts", "irc"}:
        report.add(
            "ERROR",
            "optimize.lib",
            "This runtype is not wired to the SciPy optimizer map.",
            value=f"{lib}/{runtype}",
            expected="geometric",
            action="Use [optimize] lib=geometric for runtype=ts/irc.",
        )

    if runtype == "neb" and lib not in {"geometric", "oqp"}:
        report.add(
            "ERROR",
            "optimize.lib",
            "NEB is wired through geomeTRIC and the oqp optimizer.",
            value=lib,
            expected="geometric or oqp",
            action="Set [optimize] lib=geometric or lib=oqp for runtype=neb.",
        )

    if lib == "geometric" and runtype not in {"optimize", "meci", "mecp", "ts", "irc", "neb"}:
        report.add(
            "ERROR",
            "optimize.lib",
            "geomeTRIC is currently connected only to state-specific geometry optimization, MECI, MECP, TS, IRC, and NEB.",
            value=f"{lib}/{runtype}",
            expected="optimize, meci, mecp, ts, irc, or neb",
            action="Use [input] runtype=optimize/meci/mecp/ts/irc/neb or choose scipy for this runtype.",
        )

    if lib == "oqp" and runtype not in {"optimize", "ts", "meci", "mecp", "tci", "neb", "irc", "mep"}:
        report.add(
            "ERROR",
            "optimize.lib",
            "The oqp optimizer currently supports optimize, ts, meci, mecp, tci, neb, irc, and mep.",
            value=f"{lib}/{runtype}",
            expected="optimize, ts, meci, mecp, tci, neb, irc, or mep",
            action="Choose a supported oqp runtype.",
        )

    if runtype == "tci":
        if method != "tdhf":
            report.add(
                "ERROR", "input.method",
                "Three-state CI (tci) optimization requires method=tdhf.",
                value=method, expected="tdhf",
                action="Set [input] method=tdhf and configure [tdhf].",
            )
        if not (istate < jstate < kstate):
            report.add(
                "ERROR", "optimize.kstate",
                "TCI requires istate < jstate < kstate.",
                value=f"{istate}/{jstate}/{kstate}",
                expected="istate < jstate < kstate",
                action="Set three increasing state indices (e.g. 1/2/3).",
            )
        if lib != "oqp":
            report.add(
                "ERROR", "optimize.lib",
                "Three-state CI (tci) is currently wired only through the oqp optimizer.",
                value=lib, expected="oqp",
                action="Set [optimize] lib=oqp for runtype=tci.",
            )


def _check_neb(config: dict[str, Any], report: CheckReport,
               input_dir: str | None = None) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    istate = _get(config, "optimize", "istate", 0)
    product = _get(config, "neb", "product", "")
    nimage = _get(config, "neb", "nimage", 5)

    if method == "hf" and istate != 0:
        report.add(
            "ERROR",
            "optimize.istate",
            "HF/DFT NEB currently supports only the ground state.",
            value=istate,
            expected="0",
            action="Set [optimize] istate=0, or use method=tdhf/[tdhf] type=mrsf for excited-state NEB.",
        )

    if method not in {"hf", "tdhf"}:
        report.add(
            "ERROR",
            "input.method",
            "NEB currently supports HF/DFT and TDHF/MRSF state-specific surfaces.",
            value=method,
            expected="hf or tdhf",
            action="Use method=hf with istate=0 or method=tdhf with a valid target state.",
        )

    if not product:
        report.add(
            "ERROR",
            "neb.product",
            "NEB product endpoint is missing.",
            expected="XYZ filename",
            action="Set [neb] product to a product-endpoint XYZ file.",
        )
    else:
        # The product endpoint may be given relative to the input file (where it
        # is normally stored beside the .inp), not just the current directory.
        candidates = [os.path.abspath(str(product))]
        if input_dir and not os.path.isabs(str(product)):
            candidates.append(os.path.abspath(os.path.join(input_dir, str(product))))
        if not any(os.path.exists(c) for c in candidates):
            report.add(
                "ERROR",
                "neb.product",
                "NEB product endpoint file does not exist.",
                value=product,
                action="Fix the product path, or place the XYZ file beside the input "
                       "file or in the working directory.",
            )

    if nimage < 3:
        report.add(
            "ERROR",
            "neb.nimage",
            "NEB requires at least reactant, one intermediate, and product images.",
            value=nimage,
            expected=">= 3",
            action="Set [neb] nimage=3 or larger.",
        )


def analytic_hessian_capability(config: dict[str, Any]) -> tuple[str, str]:
    """Return analytic-Hessian capability status and a precise reason.

    This checker is intentionally conservative: enabling a dispatch scaffold must
    not imply broad scientific support. Unsupported analytic Hessians fail at
    input-check time instead of silently falling back to numerical Hessians.
    """

    method = _as_lower(_get(config, "input", "method", "hf"))
    scf_type = _as_lower(_get(config, "scf", "type", "rhf"))
    td_type = _as_lower(_get(config, "tdhf", "type", "rpa"))
    functional = _as_lower(_get(config, "input", "functional", ""))
    state = _get(config, "hess", "state", 0)

    # The native analytic-Hessian derivative-integral machinery now covers the
    # features that were previously gated to the numerical Hessian:
    #   * ECP second derivatives -- libecpint deriv order 2, contracted in
    #     hf_hessian via add_ecphess + the ECP core-derivative in the CPHF response;
    #   * range-separated (CAM/LC) functionals -- the erfc-attenuated two-pass split
    #     in grd2_hess_driver (skeleton), grd2_driver (fock_deriv_contract response)
    #     and fock_jk (cphf), keyed off infos%dft%cam_flag.
    # Both are finite-difference validated for RHF/RKS, UHF/UKS and ROHF/ROKS.
    if method == "hf":
        if state != 0:
            return "unsupported_feature", "HF/DFT analytic Hessian supports only hess.state=0 in this scaffold."
        if scf_type == "rhf":
            return "supported", "Native OpenQP HF/DFT ground-state analytic Hessian dispatch is enabled."
        if scf_type in ("uhf", "rohf"):
            return "supported", f"Native OpenQP open-shell ({scf_type.upper()}) HF/DFT analytic Hessian dispatch is enabled."
        return "unsupported_scf_type", "Native analytic Hessian supports RHF/RKS, UHF/UKS and ROHF/ROKS references for this scftype. Use [hess] type=numerical."

    if method == "tdhf":
        if td_type == "mrsf":
            return "unsupported_tdhf_type", "MRSF-TDDFT analytic Hessian is not implemented; use type=numerical until the MRSF gradient/Z-vector finite-difference baseline is validated."
        if td_type == "umrsf":
            return "unsupported_tdhf_type", "UMRSF-TDDFT analytic Hessian is not implemented; use type=numerical until UMRSF-TDDFT gradients/Z-vectors are implemented and finite-difference validated."
        if td_type == "sf":
            return "unsupported_tdhf_type", "SF-TDDFT analytic Hessian is not implemented; use type=numerical until the SF gradient/Z-vector finite-difference baseline is validated."
        if td_type in {"tda", "rpa"}:
            return "unsupported_tdhf_type", f"TDDFT analytic Hessian is not implemented yet for tdhf.type={td_type}."
        return "unsupported_tdhf_type", f"Analytic Hessian does not support tdhf.type={td_type}."

    return "unsupported_method", f"Analytic Hessian does not support input.method={method}."


def _basis_max_angular_momentum(config: dict[str, Any]) -> int | None:
    """Return max L in the configured basis, or None if it cannot be inspected."""
    try:
        import basis_set_exchange as bse
    except Exception:
        return None

    basis = _get(config, "input", "basis", "")
    system = _get(config, "input", "system", "")
    library = _get(config, "input", "library", "")
    inline_lines, xyz_path = _iter_coordinate_lines(system)
    lines = inline_lines
    if xyz_path and os.path.exists(os.path.abspath(xyz_path)):
        with open(os.path.abspath(xyz_path), "r", encoding="utf-8") as handle:
            xyz_lines = handle.read().splitlines()
        try:
            num_atoms = int(xyz_lines[0])
            lines = xyz_lines[2:2 + num_atoms]
        except (IndexError, ValueError):
            lines = xyz_lines

    if not lines:
        return None

    per_atom_basis: list[str] = []
    if basis == "library":
        mapping: dict[str, str] = {}
        for raw in library.splitlines():
            parts = raw.split()
            if len(parts) >= 2:
                mapping[parts[0]] = " ".join(parts[1:])
        for line in lines:
            parts = line.split()
            if len(parts) >= 5 and parts[4] in mapping:
                per_atom_basis.append(mapping[parts[4]])
    else:
        names = [item.strip() for item in str(basis).split(";") if item.strip()]
        if len(names) == 1:
            per_atom_basis = names * len(lines)
        else:
            per_atom_basis = names

    if len(per_atom_basis) != len(lines):
        return None

    max_l = 0
    for line, basis_name in zip(lines, per_atom_basis):
        parts = line.split()
        if not parts:
            continue
        element = parts[0]
        data = bse.get_basis(basis_name, elements=[element])
        for item in data.get("elements", {}).values():
            for shell in item.get("electron_shells", []):
                max_l = max(max_l, max(int(l) for l in shell.get("angular_momentum", [])))
    return max_l


def _check_hess(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    state = _get(config, "hess", "state", 0)
    hess_type = _as_lower(_get(config, "hess", "type", "numerical"))
    read = _get(config, "hess", "read", False)
    restart = _get(config, "hess", "restart", False)
    nproc = _get(config, "hess", "nproc", 1)
    temperatures = _as_list(_get(config, "hess", "temperature", []))

    if hess_type == "analytical":
        capability, reason = analytic_hessian_capability(config)
        if capability != "supported":
            report.add(
                "ERROR",
                "hess.type",
                reason,
                value=hess_type,
                expected="supported analytical Hessian capability",
                action="Set [hess] type=numerical or use a supported analytic-Hessian method/state.",
            )
        max_l = _basis_max_angular_momentum(config)
        if max_l is not None and max_l >= 4:
            report.add(
                "ERROR",
                "input.basis",
                "Analytical Hessian native Rys nuclear-attraction second derivatives support basis angular momentum only up to L=3.",
                value=f"max L={max_l}",
                expected="max L <= 3",
                action="Use a basis without g/higher functions for analytical Hessian, or set [hess] type=numerical.",
            )

    if method == "hf" and state > 0:
        report.add(
            "ERROR",
            "hess.state",
            "HF/DFT Hessian currently supports only state 0.",
            value=state,
            expected="0",
            action="Set hess.state=0 or switch to method=tdhf.",
        )

    if method == "tdhf" and state == 0:
        report.add(
            "ERROR",
            "hess.state",
            "TDHF Hessian cannot target state 0 in the current implementation.",
            value=state,
            expected=">= 1",
            action="Set hess.state to an excited-state index.",
        )

    if nproc < 1:
        report.add(
            "ERROR",
            "hess.nproc",
            "Number of Hessian workers must be positive.",
            value=nproc,
            expected=">= 1",
            action="Set hess.nproc to 1 or more.",
        )

    for temp in temperatures:
        if float(temp) <= 0:
            report.add(
                "ERROR",
                "hess.temperature",
                "Thermochemistry temperature must be positive.",
                value=temp,
                action="Use positive Kelvin values.",
            )

    if not read:
        _add_cpu_info(report, "hess.nproc", nproc, restart)


def _check_nac(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    td_type = _as_lower(_get(config, "tdhf", "type", "rpa"))
    nproc = _get(config, "nac", "nproc", 1)
    states = _as_list(_get(config, "nac", "states", []))

    if method != "tdhf":
        report.add(
            "ERROR",
            "input.method",
            "NAC workflows require method=tdhf.",
            value=method,
            expected="tdhf",
            action="Set [input] method=tdhf.",
        )

    if td_type != "mrsf":
        report.add(
            "ERROR",
            "tdhf.type",
            "NAC workflows currently require tdhf.type=mrsf.",
            value=td_type,
            expected="mrsf",
            action="Set [tdhf] type=mrsf.",
        )

    if nproc < 1:
        report.add(
            "ERROR",
            "nac.nproc",
            "Number of NAC workers must be positive.",
            value=nproc,
            expected=">= 1",
            action="Set nac.nproc to 1 or more.",
        )

    if not states:
        report.add(
            "WARNING",
            "nac.states",
            "No NAC state pairs were provided.",
            action="Set [nac] states with pairs like 1 2,2 3.",
            wiki=WIKI_HELP["nac.states"],
        )
    else:
        for pair in states:
            if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                report.add(
                    "ERROR",
                    "nac.states",
                    "Each NAC entry must be a pair of state indices.",
                    value=pair,
                    expected="i j",
                    action="Use syntax like states=1 2,2 3.",
                    wiki=WIKI_HELP["nac.states"],
                )
                continue
            if min(int(pair[0]), int(pair[1])) < 1:
                report.add(
                    "ERROR",
                    "nac.states",
                    "NAC state indices must be excited-state indices >= 1.",
                    value=pair,
                    action="Use positive TDHF state numbers.",
                    wiki=WIKI_HELP["nac.states"],
                )

    _add_cpu_info(report, "nac.nproc", nproc, False)


def _check_nacme(config: dict[str, Any], report: CheckReport) -> None:
    _check_nac(config, report)

    guess_file2 = _get(config, "guess", "file2", "")
    system2 = _get(config, "input", "system2", "")
    align = _as_lower(_get(config, "nac", "align", "reorder"))

    if not guess_file2 and not system2:
        report.add(
            "ERROR",
            "input.system2",
            "NACME overlap needs previous-step information from guess.file2 or input.system2.",
            action="Set [guess] file2 to a restart JSON or provide [input] system2.",
        )

    if system2:
        inline_lines, xyz_path = _iter_coordinate_lines(system2)
        if xyz_path and not os.path.exists(os.path.abspath(xyz_path)):
            report.add(
                "ERROR",
                "input.system2",
                "Referenced system2 XYZ file does not exist.",
                value=xyz_path,
                action="Fix the system2 file path.",
            )
        if not xyz_path and not inline_lines:
            report.add(
                "ERROR",
                "input.system2",
                "system2 is set but empty.",
                action="Provide a previous geometry or remove system2.",
            )

    if align not in {"no", "reorder"}:
        report.add(
            "WARNING",
            "nac.align",
            "Only align=no and align=reorder are clearly handled by the Python layer.",
            value=align,
            action="Use align=reorder unless you have implemented another mode downstream.",
        )


def _add_cpu_info(report: CheckReport, path: str, nproc: int, restart: bool) -> None:
    if restart:
        return

    omp = _omp_threads(report)
    ncpu = multiprocessing.cpu_count()
    mpi = MPIManager()
    requested = mpi.size * omp if mpi.use_mpi else nproc * omp
    mode = "single node via MPI" if mpi.use_mpi else "local workers"
    severity = "WARNING" if requested > ncpu else "INFO"
    report.add(
        severity,
        path,
        f"Requested approximately {requested} CPU threads through {mode}; local machine reports {ncpu} CPUs.",
        value=requested,
        action="Lower nproc or OMP_NUM_THREADS if this run oversubscribes the node.",
    )


def check_input_values(
    config: dict[str, Any],
    *,
    raise_error: bool = True,
    emit: bool = True,
    input_dir: str | None = None,
) -> CheckReport:
    """Validate an already parsed OpenQP config and return a diagnostic report.

    ``input_dir`` is the directory of the input file, used to resolve paths
    (e.g. the NEB product endpoint) that are stored relative to the input file
    rather than the current working directory.
    """

    report = CheckReport()
    method = _as_lower(_get(config, "input", "method", "hf"))

    if method not in METHODS:
        report.add(
            "ERROR",
            "input.method",
            "Unknown electronic structure method.",
            value=method,
            expected=", ".join(sorted(METHODS)),
            action="Choose hf or tdhf.",
            wiki=WIKI_HELP["input.method"],
        )

    _check_system(config, report)
    _check_basis(config, report)
    _check_guess(config, report)
    _check_scf(config, report)
    _check_symmetry(config, report)
    _check_tdhf(config, report)
    _check_properties(config, report)
    _check_requested_states(config, report)
    _check_runtype(config, report, input_dir)

    if _get(config, "input", "d4", False) and not _get(config, "input", "functional", ""):
        report.add(
            "WARNING",
            "input.d4",
            "D4 dispersion is enabled but no DFT functional is set.",
            value=True,
            action="Set [input] functional for DFT or disable d4 for pure HF.",
        )

    if emit and report.diagnostics:
        print(report.to_text())

    if raise_error and not report.ok:
        raise SystemExit(1)

    return report
