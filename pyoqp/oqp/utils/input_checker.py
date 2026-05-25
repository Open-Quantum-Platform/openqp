"""OpenQP input checker with actionable diagnostics."""

from __future__ import annotations

from dataclasses import dataclass, field
import multiprocessing
import os
from typing import Any

from oqp.utils.mpi_utils import MPIManager


SUPPORTED_RUNTYPES = {
    "energy", "grad", "hess", "nac", "nacme", "bp", "optimize",
    "meci", "mecp", "mep", "ts", "prop", "data",
}
NOT_AVAILABLE_RUNTYPES = {"soc", "neb"}
ALL_RUNTYPES = SUPPORTED_RUNTYPES | NOT_AVAILABLE_RUNTYPES

METHODS = {"hf", "tdhf"}
SCF_TYPES = {"rhf", "rohf", "uhf"}
TDHF_TYPES = {"rpa", "tda", "sf", "mrsf", "umrsf"}
GUESS_TYPES = {"huckel", "hcore", "json", "auto", "pyscf", "sad", "sap"}
SCF_CONVERGERS = {"diis", "soscf", "trah"}
OPTIONAL_SCF_CONVERGERS = SCF_CONVERGERS | {"none", ""}
DIIS_TYPES = {"none", "cdiis", "ediis", "adiis", "vdiis"}
PCM_BACKENDS = {"ddx", "pcmsolver"}
PCM_MODES = {"reference_scf", "reference_scf_plus_post_state", "post_state_correction"}
PCM_MODELS = {"ddcosmo", "ddpcm", "ddlpb", "iefpcm", "cpcm"}
PCM_BACKEND_MODELS = {
    "ddx": {"ddcosmo", "ddpcm", "ddlpb"},
    "pcmsolver": {"iefpcm", "cpcm"},
}
OPT_LIBS = {"scipy", "dlfind"}
SCIPY_OPTIMIZERS = {"bfgs", "cg", "l-bfgs-b", "newton-cg"}
MECI_SEARCH = {"penalty", "ubp", "hybrid"}
SCF_PROPS = {"el_mom", "mulliken"}
DLFIND_SINGLE_ICOORD = {0, 1, 2, 3, 4}
DLFIND_LN_ICOORD = {10, 11, 12, 13, 14}
DLFIND_MIN_IOPT = {0, 1, 2, 3}
DLFIND_TS_IOPT = {9}
DLFIND_MECI_IMS = {1, 2, 3}
INIT_SCF_TYPES = {"no", "rhf", "uhf", "rohf", "rks", "uks", "roks"}

WIKI_HELP = {
    "input.runtype": "Use energy, grad, hess, nac, nacme, optimize, meci, mecp, mep, ts, prop, or data. soc and neb are recognized but not implemented yet.",
    "input.method": "Use method=hf for HF/DFT and method=tdhf for TDHF/TDDFT/SF/MRSF runs.",
    "input.system": "Set system to an XYZ file path or inline coordinates with one atom per indented line.",
    "input.basis": "Set basis to a basis name, a comma-separated per-atom list, or library with tagged atoms and [input] library mappings.",
    "scf.type": "RHF is for multiplicity 1 closed-shell references. SF/MRSF needs an open-shell reference, usually ROHF.",
    "tdhf.type": "Use rpa or tda for ordinary TDHF/TDDFT, sf or mrsf for spin-flip, and umrsf only with UHF.",
    "tdhf.nstate": "nstate must cover the highest excited-state index requested anywhere else in the input.",
    "guess.type": "Use json with a JSON restart file, auto for JSON-if-present otherwise Huckel, sad/sap for PySCF atomic-density/potential guesses, or pyscf to build a converged external guess.",
    "pcm.enabled": "PCM input is reserved for the planned energy-only solvent backend. Initial scope is RHF/ROHF reference_scf single-point energy; gradients and state-specific MRSF PCM are out of scope.",
    "pcm.backend": "Use backend=ddx for the preferred active ddCOSMO/ddPCM library candidate, or backend=pcmsolver for the classic PCM API candidate.",
    "pcm.mode": "Use mode=reference_scf for MRSF-compatible PCM on the RHF/ROHF reference. post_state_correction and reference_scf_plus_post_state are planned perturbative extensions.",
    "optimize.lib": "scipy supports optimize, meci, mecp, and mep. dlfind supports optimize, meci, and ts.",
    "dlfind.ims": "ims=0 is single-state, ims=1/2/3 are MECI modes and belong to runtype=meci.",
    "nac.states": "Use state pairs such as 1 2,2 3 for NAC calculations. Each index must be a TDHF excited state.",
}


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
    guess_type = _as_lower(_get(config, "guess", "type", "huckel"))
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


def _check_pcm(config: dict[str, Any], report: CheckReport) -> None:
    pcm = config.get("pcm", {})
    if not pcm:
        return

    enabled = bool(_get(config, "pcm", "enabled", False))
    backend = _as_lower(_get(config, "pcm", "backend", "ddx"))
    mode = _as_lower(_get(config, "pcm", "mode", "reference_scf"))
    model = _as_lower(_get(config, "pcm", "model", "ddpcm"))
    epsilon = _get(config, "pcm", "epsilon", 78.3553)
    runtype = _as_lower(_get(config, "input", "runtype", "energy"))

    if backend not in PCM_BACKENDS:
        report.add(
            "ERROR",
            "pcm.backend",
            "Unknown PCM backend.",
            value=backend,
            expected=", ".join(sorted(PCM_BACKENDS)),
            action="Choose ddx or pcmsolver.",
            wiki=WIKI_HELP["pcm.backend"],
        )

    if mode not in PCM_MODES:
        report.add(
            "ERROR",
            "pcm.mode",
            "Unknown PCM coupling mode.",
            value=mode,
            expected=", ".join(sorted(PCM_MODES)),
            action="Use reference_scf for the first MRSF-compatible solvent mode.",
            wiki=WIKI_HELP["pcm.mode"],
        )

    if model not in PCM_MODELS:
        report.add(
            "ERROR",
            "pcm.model",
            "Unknown PCM model.",
            value=model,
            expected=", ".join(sorted(PCM_MODELS)),
            action="Use ddpcm/ddcosmo with backend=ddx or iefpcm/cpcm with backend=pcmsolver.",
        )

    if backend in PCM_BACKEND_MODELS and model in PCM_MODELS and model not in PCM_BACKEND_MODELS[backend]:
        report.add(
            "ERROR",
            "pcm.model",
            f"PCM model {model} is not supported by backend {backend}.",
            value=model,
            expected=", ".join(sorted(PCM_BACKEND_MODELS[backend])),
            action="Use a ddCOSMO/ddPCM/ddLPB model with backend=ddx, or an IEFPCM/CPCM model with backend=pcmsolver.",
            wiki=WIKI_HELP["pcm.backend"],
        )

    if float(epsilon) <= 1.0:
        report.add(
            "ERROR",
            "pcm.epsilon",
            "PCM dielectric constant must be greater than 1.",
            value=epsilon,
            action="Use a physical solvent dielectric, e.g. 78.3553 for water.",
        )

    if enabled:
        report.add(
            "ERROR",
            "pcm.enabled",
            "PCM input parsing is scaffolded but runtime solvent coupling is not implemented yet.",
            value=True,
            expected="False until the backend Fock/energy coupling is implemented",
            action="Leave pcm.enabled=false or continue development by adding the RHF/ROHF reference_scf backend.",
            wiki=WIKI_HELP["pcm.enabled"],
        )

        if runtype != "energy":
            report.add(
                "ERROR",
                "input.runtype",
                "The first PCM implementation is scoped to single-point energies only.",
                value=runtype,
                expected="energy",
                action="Use runtype=energy; gradients/optimizations require a separate analytic-gradient implementation.",
                wiki=WIKI_HELP["pcm.enabled"],
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
            action="Choose rpa, tda, sf, mrsf, or umrsf.",
            wiki=WIKI_HELP["tdhf.type"],
        )
        return

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
            "UMRSF requires a UHF reference.",
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

    if method == "tdhf" and runtype == "grad" and grad_states:
        highest_grad = max(int(state) for state in grad_states if state not in (None, ""))
        nstate = _get(config, "tdhf", "nstate", 1)
        if highest_grad == nstate:
            report.add(
                "WARNING",
                "tdhf.nstate",
                "The requested gradient root is exactly the highest computed TD state.",
                value=nstate,
                expected=f">= {highest_grad + 1}",
                action="Consider increasing tdhf.nstate by 1 to avoid missing degenerate states.",
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


def _check_runtype(config: dict[str, Any], report: CheckReport) -> None:
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
        elif 0 in [int(state) for state in _as_list(_get(config, "properties", "grad", [])) if state not in (None, "")]:
            report.add(
                "ERROR",
                "properties.grad",
                "TDHF grad runtype cannot target state 0.",
                value=_get(config, "properties", "grad", []),
                action="Use grad=1,2,... for excited-state gradients.",
            )

    if runtype in {"optimize", "meci", "mecp", "mep", "ts"}:
        _check_optimize(config, report)

    if runtype in {"nac", "bp"}:
        _check_nac(config, report)

    if runtype == "nacme":
        _check_nacme(config, report)

    if runtype == "hess":
        _check_hess(config, report)


def _check_optimize(config: dict[str, Any], report: CheckReport) -> None:
    runtype = _as_lower(_get(config, "input", "runtype", "optimize"))
    method = _as_lower(_get(config, "input", "method", "hf"))
    lib = _as_lower(_get(config, "optimize", "lib", "scipy"))
    optimizer = _as_lower(_get(config, "optimize", "optimizer", "bfgs"))
    istate = _get(config, "optimize", "istate", 0)
    jstate = _get(config, "optimize", "jstate", 0)
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
            action="Use scipy or dlfind.",
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

    if lib == "dlfind":
        _check_dlfind(config, report)

    if lib == "scipy" and runtype == "ts":
        report.add(
            "ERROR",
            "optimize.lib",
            "Transition-state optimization is not wired to the SciPy optimizer map.",
            value=lib,
            expected="dlfind",
            action="Use [optimize] lib=dlfind for runtype=ts.",
        )

    if lib == "dlfind" and runtype not in {"optimize", "meci", "ts"}:
        report.add(
            "ERROR",
            "optimize.lib",
            "DL-FIND is not connected to this runtype.",
            value=f"{lib}/{runtype}",
            expected="optimize, meci, or ts",
            action="Switch to lib=scipy or choose a DL-FIND-supported runtype.",
        )


def _check_dlfind(config: dict[str, Any], report: CheckReport) -> None:
    runtype = _as_lower(_get(config, "input", "runtype", "optimize"))
    icoord = _get(config, "dlfind", "icoord", 3)
    iopt = _get(config, "dlfind", "iopt", 3)
    ims = _get(config, "dlfind", "ims", 0)

    if runtype == "optimize":
        if icoord not in DLFIND_SINGLE_ICOORD:
            report.add(
                "ERROR",
                "dlfind.icoord",
                "Single-state DL-FIND optimization only supports icoord 0-4.",
                value=icoord,
                action="Use icoord in 0,1,2,3,4.",
            )
        if iopt not in DLFIND_MIN_IOPT:
            report.add(
                "ERROR",
                "dlfind.iopt",
                "Single-state DL-FIND optimization only supports iopt 0-3.",
                value=iopt,
                action="Use iopt in 0,1,2,3.",
            )
        if ims != 0:
            report.add(
                "ERROR",
                "dlfind.ims",
                "Single-state optimization requires ims=0.",
                value=ims,
                action="Set ims=0.",
                wiki=WIKI_HELP["dlfind.ims"],
            )

    if runtype == "meci":
        if iopt not in DLFIND_MIN_IOPT:
            report.add(
                "ERROR",
                "dlfind.iopt",
                "DL-FIND MECI only supports iopt 0-3.",
                value=iopt,
                action="Use iopt in 0,1,2,3.",
            )
        if ims not in DLFIND_MECI_IMS:
            report.add(
                "ERROR",
                "dlfind.ims",
                "DL-FIND MECI requires ims=1, 2, or 3.",
                value=ims,
                action="Set ims to a MECI mode.",
                wiki=WIKI_HELP["dlfind.ims"],
            )
        if ims == 3 and icoord not in DLFIND_LN_ICOORD:
            report.add(
                "ERROR",
                "dlfind.icoord",
                "Lagrange-Newton MECI requires icoord 10-14.",
                value=icoord,
                action="Use icoord 10-14 with ims=3.",
            )
        if ims in {1, 2} and icoord not in DLFIND_SINGLE_ICOORD:
            report.add(
                "ERROR",
                "dlfind.icoord",
                "Penalty/gradient-projection MECI requires icoord 0-4.",
                value=icoord,
                action="Use icoord 0-4 with ims=1 or ims=2.",
            )

    if runtype == "ts":
        if iopt not in DLFIND_TS_IOPT:
            report.add(
                "ERROR",
                "dlfind.iopt",
                "Transition-state DL-FIND uses P-RFO (iopt=9).",
                value=iopt,
                action="Set [dlfind] iopt=9 for runtype=ts.",
            )
        if ims != 0:
            report.add(
                "ERROR",
                "dlfind.ims",
                "Transition-state search is not a MECI mode and requires ims=0.",
                value=ims,
                action="Set [dlfind] ims=0 for runtype=ts.",
            )
        if icoord not in DLFIND_SINGLE_ICOORD:
            report.add(
                "ERROR",
                "dlfind.icoord",
                "Transition-state DL-FIND search expects icoord 0-4.",
                value=icoord,
                action="Use icoord 0-4 for TS optimization.",
            )


def _check_hess(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    state = _get(config, "hess", "state", 0)
    hess_type = _as_lower(_get(config, "hess", "type", "numerical"))
    read = _get(config, "hess", "read", False)
    restart = _get(config, "hess", "restart", False)
    nproc = _get(config, "hess", "nproc", 1)
    temperatures = _as_list(_get(config, "hess", "temperature", []))

    if hess_type == "analytical":
        report.add(
            "ERROR",
            "hess.type",
            "Analytical Hessian is not implemented; the runtime exits for this case.",
            value=hess_type,
            expected="numerical",
            action="Set [hess] type=numerical.",
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
) -> CheckReport:
    """Validate an already parsed OpenQP config and return a diagnostic report."""

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
    _check_pcm(config, report)
    _check_scf(config, report)
    _check_tdhf(config, report)
    _check_properties(config, report)
    _check_requested_states(config, report)
    _check_runtype(config, report)

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
