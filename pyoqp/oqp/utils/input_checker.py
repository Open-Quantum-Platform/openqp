"""OpenQP input checker with actionable diagnostics."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import multiprocessing
import os
import re
from typing import Any

from oqp.utils.mpi_utils import MPIManager


SUPPORTED_RUNTYPES = {
    "energy", "grad", "hess", "nac", "nacme", "bp", "optimize",
    "meci", "mecp", "tci", "mep", "ts", "irc", "neb", "prop", "data", "ekt", "soc",
    "namd",
}
NOT_AVAILABLE_RUNTYPES = {"md"}
ALL_RUNTYPES = SUPPORTED_RUNTYPES | NOT_AVAILABLE_RUNTYPES
METHODS = {
    "hf", "tdhf", "mp2", "dftb", "xtb",
    "fci", "casci", "casscf", "sa-casscf", "sacasscf",
    "caspt2", "ms-caspt2", "mscaspt2", "xms-caspt2", "xmscaspt2",
    "mrmp2", "mcqdpt2", "xmcqdpt2",
}
SCF_TYPES = {"rhf", "rohf", "uhf"}
TDHF_TYPES = {"rpa", "tda", "sf", "mrsf", "umrsf", "mrsf_ekt_ip", "mrsf_ekt_ea"}
MP2_VARIANTS = {
    "mp2", "conventional",
    "scs", "scs-mp2",
    "sos", "sos-mp2",
    "os", "os-mp2", "opposite-spin",
    "ss", "ss-mp2", "same-spin", "sss", "sss-mp2",
    "scs-mi", "scs-mi-mp2",
    "custom",
}
DFTB_BACKENDS = {"auto", "native", "probe"}
DFTB_TYPES = {
    "auto", "ground", "dftb", "dftb0", "ground_noscc", "noscc",
    "tddftb", "tda", "td-dftb", "sf", "sftddftb", "sf-tddftb",
    "mrsf", "mrsftddftb", "mrsf-tddftb",
}
DFTB_SCC_MIXERS = {"linear", "anderson", "pulay", "broyden", "auto", "diis", "trust", "trah"}

# Canonical model keywords are "dtcam", "dtcam2"/"dtcam-erf" and "ob2".  The
# historical spellings stay ACCEPTED ALIASES so committed inputs keep working;
# "dftb+" in particular was renamed because DFTB+ is a different program and
# the preset is really the conventional OB2 / LC-DFTB2 protocol.  This set must
# stay in sync with openqp_dftb_preset_by_name in the native library.
DFTB_MODELS = {
    # canonical
    "dtcam",
    "dtcam2", "dtcam-erf", "dtcam_erf", "dtcamerf",
    "ob2",
    # legacy aliases
    "dtcam-tb", "dtcam_tb", "dtcamtb",
    "dtcam-tb2", "dtcam_tb2", "dtcamtb2",
    "dtcam-tb-erf", "dtcam_tb_erf", "dtcamtberf",
    "dftb+", "dftbplus", "dftb_plus",
}
# Production defaults when no model= is named and no preset-locked key is
# tuned: MRSF-TDDFTB runs the tuned DTCAM-TB operator by default (it restores
# the ionic-bright / covalent-dark crossing the plain protocol misses); the
# other OPEN-SHELL SCC routes (SF, open-shell ground) run the OB2
# compatibility protocol (conventional LC-DFTB2, SPC=0.5) by default.  Either
# operator can be requested explicitly with [dftb] model=dtcam or
# model=ob2.  Closed-shell routes (singlet ground, TD-DFTB) and DFTB0 stay
# preset-free: LC-DFTB2 is not implemented for the restricted reference, and
# the native library rejects lc_ground_state there.
DFTB_MODEL_DEFAULT_MRSF = "dtcam"
DFTB_MODEL_DEFAULT_OPEN_SHELL = "ob2"
# Keys a [dftb] model preset fixes; the checker refuses to combine them with
# model= (the preset overrides them inside openqp-dftb, so a user-tuned value
# would be silently discarded).  ``scc_mixer`` is intentionally excluded:
# selecting the first-try convergence algorithm (Broyden versus charge/spin
# TRAH) is a numerical recovery choice, and the native C API applies that
# explicit override after loading the operator preset.
DFTB_MODEL_LOCKED_KEYS = (
    "omega", "cam_alpha", "cam_beta", "lc_gamma", "lc_ground_state",
    "w_scale", "response_w_scale", "response_omega", "response_cam_alpha",
    "response_cam_beta", "c_mrsf", "c_mrsf_oo", "response_global_hybrid",
    "onsite_exchange_scale", "spc", "spc_coco", "spc_ovov", "spc_coov",
    "onsite_ss", "onsite_sp", "onsite_pp",
    "mrsf_shift_oo", "mrsf_shift_co", "mrsf_shift_ov", "mrsf_shift_cv",
    # The preset also fixes the numerical protocol (openqp-dftb
    # openqp_dftb_apply_dtcam_preset): Broyden mixing, SCC budget and
    # tolerance, the Davidson response solver, and the Z-vector gradient
    # path. Keys the preset does NOT touch (response_tolerance,
    # response_max_iterations, response_max_subspace, nstate, ...) stay
    # user-tunable and are deliberately absent here.
    "scc_mixing", "scc_history", "scc_max_step",
    "scc_tolerance", "max_scc_iterations", "response_solver", "zvector",
)

# Top-level [input] method shorthands for the DFTB response families.  Writing
# method=mrsf-tddftb / sf-tddftb / tddftb is equivalent to method=dftb with the
# matching [dftb]/[tdhf] type and the reference each family needs (ROKS triplet
# for the spin-flip families, closed-shell RHF for plain TD-DFTB).  Expanded
# before validation so the rest of the checker sees the canonical method=dftb.
DFTB_METHOD_ALIASES = {
    "mrsf-tddftb": "mrsf", "mrsftddftb": "mrsf", "mrsf-dftb": "mrsf",
    "sf-tddftb": "sf", "sftddftb": "sf", "sf-dftb": "sf",
    "tddftb": "tddftb", "td-dftb": "tddftb", "tda-dftb": "tddftb",
    "tddftb-tda": "tddftb", "tda-tddftb": "tddftb",
}


def expand_dftb_method_alias(config: dict[str, Any]) -> str:
    """Expand a top-level [input] method DFTB-response alias in place.

    method=mrsf-tddftb / sf-tddftb / tddftb becomes method=dftb plus the
    canonical [dftb] type and the family's reference ([scf] type + reference
    multiplicity, [tdhf] type).  Because the parsed config always carries schema
    defaults, the implied keys are set authoritatively (the alias defines them);
    use the explicit method=dftb route for a non-standard reference.  Returns the
    canonical dftb type when an alias was expanded, else "".
    """
    if not isinstance(config, dict):
        return ""
    method = _as_lower(_get(config, "input", "method", ""))
    dftb_type = DFTB_METHOD_ALIASES.get(method)
    if dftb_type is None:
        return ""
    config.setdefault("input", {})["method"] = "dftb"
    dftb = config.setdefault("dftb", {})
    dftb["type"] = dftb_type
    dftb["backend"] = _get(config, "dftb", "backend", "native") or "native"
    tdhf = config.setdefault("tdhf", {})
    scf = config.setdefault("scf", {})
    if dftb_type in ("sf", "mrsf"):
        tdhf["type"] = dftb_type
        scf["type"] = "rohf"
        scf["multiplicity"] = 3
    else:  # plain closed-shell TD-DFTB (TDA)
        tdhf["type"] = "tda"
        scf["type"] = "rhf"
        scf["multiplicity"] = 1
    return dftb_type


# New-generation operator keys the probe CLI never forwards (native only).
DFTB_NATIVE_ONLY_KEYS = (
    "c_mrsf", "c_mrsf_oo", "response_global_hybrid", "onsite_exchange_scale",
    "w_scale", "response_w_scale", "response_omega", "response_cam_alpha",
    "response_cam_beta", "onsite_ss", "onsite_sp", "onsite_pp",
)

PLANNED_WAVEFUNCTION_METHODS = set()
FCI_INTEGRAL_BACKENDS = {"native"}
FCI_SOLVERS = {"auto", "dense", "davidson"}
TARGET_SPIN_LABELS = {
    "singlet",
    "doublet",
    "triplet",
    "quartet",
    "quintet",
    "sextet",
    "septet",
    "octet",
    "nonet",
}
CAS_ORBITAL_SOURCES = {"rhf", "json"}
CAS_LOCALIZE_OPTIONS = {"none"}
CAS_SORT_OPTIONS = {"energy", "none"}
CI_ROOT_TRACKING_OPTIONS = {"energy"}
PT2_METHOD_ALIASES = {
    "caspt2": "caspt2",
    "ms-caspt2": "ms-caspt2",
    "mscaspt2": "ms-caspt2",
    "xms-caspt2": "xms-caspt2",
    "xmscaspt2": "xms-caspt2",
    "mrmp2": "mrmp2",
    "mcqdpt2": "mcqdpt2",
    "xmcqdpt2": "xmcqdpt2",
}
PT2_VARIANTS = {"auto", "caspt2", "ms-caspt2", "xms-caspt2",
                "mrmp2", "mcqdpt2", "xmcqdpt2"}
# QDPT (GAMESS-convention) family groupings used by the consistency checks
PT2_SINGLE_STATE_METHODS = {"caspt2", "mrmp2"}
PT2_MS_METHODS = {"ms-caspt2", "mcqdpt2"}
PT2_XMS_METHODS = {"xms-caspt2", "xmcqdpt2"}
PT2_QDPT_METHODS = {"mrmp2", "mcqdpt2", "xmcqdpt2"}
PT2_REFERENCES = {"casci", "casscf"}
PT2_CONTRACTIONS = {"uncontracted", "none", "full", "strong", "sc", "sc-nevpt2",
                    "ic", "internally-contracted"}
PT2_STRONG_CONTRACTIONS = {"strong", "sc", "sc-nevpt2", "ic", "internally-contracted"}
PT2_MULTISTATE_MODES = {"auto", "none", "ms", "xms"}
STATE_AVERAGE_SPIN_BLOCKS = {"diagnostic"}
STATE_AVERAGE_ROOT_TRACKING = {"overlap"}
GUESS_TYPES = {"huckel", "modhuckel", "hcore", "json", "auto", "sap", "minao"}
SCF_CONVERGERS = {"diis", "soscf", "trah", "auto", "ml"}
OPTIONAL_SCF_CONVERGERS = SCF_CONVERGERS | {"none", ""}
DIIS_TYPES = {"none", "cdiis", "ediis", "adiis", "vdiis"}
PCM_BACKENDS = {"ddx", "pcmsolver"}
PCM_MODES = {"reference_scf", "reference_scf_plus_post_state", "post_state_correction"}
PCM_MODELS = {"ddcosmo", "ddpcm", "ddlpb", "iefpcm", "cpcm"}
PCM_BACKEND_MODELS = {
    "ddx": {"ddcosmo", "ddpcm", "ddlpb"},
    "pcmsolver": {"iefpcm", "cpcm"},
}
OPT_LIBS = {"scipy", "geometric", "oqp"}
SCIPY_OPTIMIZERS = {"bfgs", "cg", "l-bfgs-b", "newton-cg"}
MECI_SEARCH = {"auto", "penalty", "ubp", "hybrid", "baeka"}
SCF_PROPS = {"el_mom", "mulliken", "lowdin", "resp", "nmr"}
NMR_GAUGES = {"cgo", "giao"}
INIT_SCF_TYPES = {"no", "rhf", "uhf", "rohf", "rks", "uks", "roks"}

WIKI_HELP = {
    "input.runtype": "Use energy, ekt, grad, hess, nac, nacme, optimize, meci, mecp, mep, ts, irc, neb, soc, prop, or data. md is recognized but not yet implemented.",
    "input.method": "Use method=hf for HF/DFT, method=tdhf for TDHF/TDDFT/SF/MRSF, method=mp2 for ground-state MP2, method=dftb for the optional OpenQP-DFTB backend, method=xtb for the optional OpenQP-xTB (LC-GFN1-xTB) backend, method=fci for legacy FCI-style CI, method=casci for fixed-orbital active-space CI, method=casscf for native CASSCF macroiterations, method=sa-casscf with [state_average] enabled=true for native SA-CASSCF, method=caspt2 for native determinant-space state-specific PT2, method=ms-caspt2 for native determinant-space multistate PT2, method=xms-caspt2 for native determinant-space extended multistate PT2, or the GAMESS-convention QDPT family: method=mrmp2 (single-state), method=mcqdpt2 (multistate, single-set), method=xmcqdpt2 (Granovsky extended H0).",
    "mp2.variant": "Use mp2, scs-mp2, sos-mp2, os-mp2, ss-mp2, scs-mi-mp2, or custom with explicit OS/SS scales.",
    "input.system": "Set system to an XYZ file path or inline coordinates with one atom per indented line.",
    "input.basis": "Set basis to a basis name, a comma-separated per-atom list, or library with tagged atoms and [input] library mappings.",
    "scf.type": "RHF is for multiplicity 1 closed-shell references. SF/MRSF needs an open-shell reference, usually ROHF.",
    "tdhf.type": "Use rpa or tda for ordinary TDHF/TDDFT, sf or mrsf for spin-flip, umrsf only with UHF, and legacy mrsf_ekt_ip/mrsf_ekt_ea only with energy runtype. EKT analysis must use [input] runtype=ekt with [tdhf] type=mrsf and [ekt] IP, EA, or both.",
    "tdhf.nstate": "nstate must cover the highest excited-state index requested anywhere else in the input.",
    "guess.type": "Use huckel or modhuckel (weighted Wolfsberg-Helmholz) for native extended-Huckel guesses, hcore for the bare core Hamiltonian, sap for the native superposition-of-atomic-potentials guess, minao for projected minimal-basis densities, json with a JSON restart file, or auto for JSON-if-present otherwise Huckel.",
    "pcm.enabled": "PCM input is reserved for the planned energy-only solvent backend. Initial scope is RHF/ROHF reference_scf single-point energy; gradients and state-specific MRSF PCM are out of scope.",
    "pcm.backend": "Use backend=ddx for the preferred active ddCOSMO/ddPCM library candidate, or backend=pcmsolver for the classic PCM API candidate.",
    "pcm.mode": "Use mode=reference_scf for MRSF-compatible PCM on the RHF/ROHF reference. post_state_correction and reference_scf_plus_post_state are planned perturbative extensions.",
    "optimize.lib": "oqp is the default optimizer backend: the built-in NumPy/SciPy optimizer (redundant internals/DLC/TRIC + restricted-step RFO) supporting optimize, ts, meci, mecp, tci, neb, irc, and mep. geometric (the external geomeTRIC package) supports optimize, MECI, MECP, TS, IRC, and NEB. scipy supports optimize, meci, mecp, and mep.",
    "nac.states": "Use state pairs such as 1 2,2 3 for NAC calculations. Each index must be a TDHF excited state.",
    "dftb.parameter_path": "Set [dftb] parameter_path to an .opdftb file or an SKF parameter directory; OPENQP_DFTB_PARAMETER_PATH may also provide it. When neither is set, the bundled openqp-dftb parameter set (OB2W0PT3 with the official spinw.txt) is used automatically if installed.",
    "dftb.namd": "OpenQP-DFTB can provide NAMD surface data through runtype=data. NAC/NACME requires an implemented DFTB state-overlap backend.",
    "xtb.parameter_path": "Set [xtb] parameter_path to an .opxtb parameter file (converter-generated); OPENQP_XTB_PARAMETER_PATH may also provide it.",
    "xtb.namd": "OpenQP-xTB can provide NAMD surface data through runtype=data. NAC/NACME requires an implemented xTB state-overlap backend.",
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


def _check_choice_literal(
    value: Any,
    path: str,
    choices: set[str],
    report: CheckReport,
    *,
    message: str,
    action: str,
    fallback: str,
    default_if_none: bool = False,
    default_if_blank: bool = False,
) -> str:
    if value is None and default_if_none:
        return fallback
    if isinstance(value, str):
        if not value.strip() and default_if_blank:
            return fallback
        normalized = value.strip().lower()
        if normalized in choices:
            return normalized
        report_value = normalized
    else:
        report_value = value

    report.add(
        "ERROR",
        path,
        message,
        value=report_value,
        expected=", ".join(sorted(choices)),
        action=action,
    )
    return fallback


def _parse_bool_literal(value: Any) -> tuple[bool, bool]:
    if isinstance(value, bool):
        return True, value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "t", "yes", "y", "on", ".true."}:
            return True, True
        if normalized in {"0", "false", "f", "no", "n", "off", ".false.", ""}:
            return True, False
        return False, False
    if isinstance(value, int) and value in (0, 1):
        return True, bool(value)
    return False, False


def _check_bool_literal(
    value: Any,
    path: str,
    report: CheckReport,
    *,
    action: str = "Use true or false.",
) -> bool:
    return _validate_bool_literal(value, path, report, action=action)[1]


def _validate_bool_literal(
    value: Any,
    path: str,
    report: CheckReport,
    *,
    action: str = "Use true or false.",
) -> tuple[bool, bool]:
    valid, parsed = _parse_bool_literal(value)
    if valid:
        return True, parsed
    report.add(
        "ERROR",
        path,
        "Boolean option must be true or false.",
        value=value,
        expected="true or false",
        action=action,
    )
    return False, False


def _parse_int_literal(value: Any) -> tuple[bool, int]:
    if isinstance(value, bool):
        return False, 0
    if isinstance(value, int):
        return True, int(value)
    if isinstance(value, float):
        if math.isfinite(value) and value.is_integer():
            return True, int(value)
        return False, 0
    if isinstance(value, str):
        raw = value.strip()
        if not raw:
            return False, 0
        try:
            return True, int(raw, 10)
        except ValueError:
            return False, 0
    return False, 0


def _check_int_literal(
    value: Any,
    path: str,
    report: CheckReport,
    *,
    minimum: int | None = None,
    expected: str = "integer",
    action: str = "Use an integer value.",
) -> tuple[bool, int]:
    valid, parsed = _parse_int_literal(value)
    if not valid:
        report.add(
            "ERROR",
            path,
            "Integer option must be an exact integer.",
            value=value,
            expected=expected,
            action=action,
        )
        return False, 0
    if minimum is not None and parsed < minimum:
        report.add(
            "ERROR",
            path,
            "Integer option is below the allowed minimum.",
            value=value,
            expected=f">= {minimum}",
            action=action,
        )
        return False, parsed
    return True, parsed


def _check_active_electron_literal(
    value: Any,
    path: str,
    report: CheckReport,
) -> tuple[bool, int]:
    expected = "non-negative integer or nalpha,nbeta pair"
    action = "Use 0 for default inference, a total count, or an alpha,beta pair such as 1,1."
    if isinstance(value, str) and "," in value:
        values = _as_list(value)
    elif isinstance(value, (list, tuple)):
        values = list(value)
    else:
        return _check_int_literal(
            value,
            path,
            report,
            minimum=0,
            expected=expected,
            action=action,
        )

    if len(values) != 2:
        report.add(
            "ERROR",
            path,
            "Active electron pair must contain exactly alpha and beta counts.",
            value=value,
            expected=expected,
            action=action,
        )
        return False, 0

    total = 0
    for item in values:
        valid, parsed = _check_int_literal(
            item,
            path,
            report,
            minimum=0,
            expected=expected,
            action=action,
        )
        if not valid:
            return False, 0
        total += parsed
    return True, total


def _check_float_literal(
    value: Any,
    path: str,
    report: CheckReport,
    *,
    expected: str = "finite number",
    action: str = "Use a finite numeric value.",
) -> tuple[bool, float]:
    if isinstance(value, bool):
        report.add(
            "ERROR",
            path,
            "Numeric option must be a finite real number, not a boolean.",
            value=value,
            expected=expected,
            action=action,
        )
        return False, 0.0
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        report.add(
            "ERROR",
            path,
            "Numeric option must be a finite real number.",
            value=value,
            expected=expected,
            action=action,
        )
        return False, 0.0
    if not math.isfinite(numeric):
        report.add(
            "ERROR",
            path,
            "Numeric option must be finite.",
            value=value,
            expected=expected,
            action=action,
        )
        return False, numeric
    return True, numeric


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
        # A QM/MM reference such as ``mol.pdb 0 1 2`` also has three finite
        # numeric fields after its first token.  Recognize supported geometry
        # filenames before applying the inline-coordinate heuristic.
        path_token, _ = _split_geometry_reference(first_line)
        if path_token.lower().endswith((".xyz", ".pdb")):
            return [], first_line

        fields = first_line.split()
        try:
            is_coordinate = len(fields) >= 4 and all(
                math.isfinite(float(value)) for value in fields[1:4]
            )
        except ValueError:
            is_coordinate = False
        if is_coordinate:
            return [line for line in lines if line.strip()], None
        return [], first_line
    return [line for line in lines[1:] if line.strip()], None


def _split_geometry_reference(reference: str) -> tuple[str, list[str]]:
    """Split an XYZ/PDB reference without breaking paths that contain spaces."""
    stripped = reference.strip()
    lower = stripped.lower()
    for suffix in (".xyz", ".pdb"):
        end = lower.find(suffix)
        if end >= 0:
            end += len(suffix)
            path = stripped[:end].strip()
            trailing = stripped[end:].strip().split()
            return path, trailing
    tokens = stripped.split()
    return (tokens[0], tokens[1:]) if tokens else ("", [])


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


def _check_qm_atom_indices(tokens: list[str], report: CheckReport) -> None:
    """Validate the QM-region atom indices that follow a .pdb path.

    Mirrors the parsing in ``oqpdata.read_system``: each token is either a
    single index or a ``start-end`` range. Indices must be non-negative and
    unique.
    """
    if not tokens:
        report.add(
            "ERROR",
            "input.system",
            "QM/MM PDB geometry has no QM atom indices.",
            expected="one or more atom indices, e.g. '9 10 17-19'",
            action="List the QM-region atom indices after the .pdb path.",
            wiki=WIKI_HELP["input.system"],
        )
        return

    indices: list[int] = []
    for tok in tokens:
        if "-" in tok:
            parts = tok.split("-")
            if len(parts) != 2 or not all(p.strip().isdigit() for p in parts):
                report.add(
                    "ERROR",
                    "input.system",
                    "Invalid QM atom index range.",
                    value=tok,
                    expected="start-end with non-negative integers, e.g. '17-19'",
                    action="Use an ascending integer range for the QM atom selection.",
                    wiki=WIKI_HELP["input.system"],
                )
                return
            start, end = int(parts[0]), int(parts[1])
            indices.extend(range(start, end + 1))
        elif tok.isdigit():
            indices.append(int(tok))
        else:
            report.add(
                "ERROR",
                "input.system",
                "Invalid QM atom index.",
                value=tok,
                expected="a non-negative integer or 'start-end' range",
                action="List QM atom indices as integers or ranges after the .pdb path.",
                wiki=WIKI_HELP["input.system"],
            )
            return

    if len(indices) != len(set(indices)):
        report.add(
            "ERROR",
            "input.system",
            "Repeated entries in the QM atom list.",
            action="Remove duplicate QM atom indices from the system selection.",
            wiki=WIKI_HELP["input.system"],
        )


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
        # Only the first whitespace-delimited token is the file path. QM/MM
        # runs append the QM-region atom indices after a .pdb path
        # (e.g. "mol.pdb 9 10 17-19"), so validate the path token alone rather
        # than the whole line (which never exists as a file).
        path_token, trailing_tokens = _split_geometry_reference(xyz_path)
        resolved = os.path.abspath(path_token)
        if path_token.lower().endswith(".pdb"):
            if not os.path.exists(resolved):
                report.add(
                    "ERROR",
                    "input.system",
                    "Referenced PDB file does not exist.",
                    value=path_token,
                    action="Fix the path or place the PDB file in the working directory.",
                    wiki=WIKI_HELP["input.system"],
                )
            else:
                _check_qm_atom_indices(trailing_tokens, report)
        elif not os.path.exists(resolved):
            report.add(
                "ERROR",
                "input.system",
                "Referenced XYZ file does not exist.",
                value=path_token,
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
    elif mode != "reference_scf":
        report.add(
            "ERROR",
            "pcm.mode",
            "Only reference_scf PCM is in the first implementation scope.",
            value=mode,
            expected="reference_scf",
            action="Do not request post-state or nonequilibrium PCM modes until those separate runtime paths are implemented and validated.",
            wiki=WIKI_HELP["pcm.mode"],
        )

    scf_type = _as_lower(_get(config, "scf", "type", "rhf"))

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

    try:
        epsilon_value = float(epsilon)
    except (TypeError, ValueError):
        report.add(
            "ERROR",
            "pcm.epsilon",
            "PCM dielectric constant must be numeric and greater than 1.",
            value=epsilon,
            action="Use a physical solvent dielectric, e.g. 78.3553 for water.",
        )
    else:
        if epsilon_value <= 1.0:
            report.add(
                "ERROR",
                "pcm.epsilon",
                "PCM dielectric constant must be greater than 1.",
                value=epsilon,
                action="Use a physical solvent dielectric, e.g. 78.3553 for water.",
            )

    if enabled:
        if scf_type not in {"rhf", "rohf"}:
            report.add(
                "ERROR",
                "scf.type",
                "PCM first scope supports RHF/ROHF reference SCF only.",
                value=scf_type,
                expected="rhf or rohf",
                action="Use RHF for closed-shell singlets or ROHF for the high-spin MRSF reference; UHF PCM is a separate future validation target.",
                wiki=WIKI_HELP["pcm.enabled"],
            )

        # The reference_scf RHF/ROHF energy path is implemented for the ddX
        # backend only. The pcmsolver backend remains a planned candidate.
        if backend != "ddx":
            report.add(
                "ERROR",
                "pcm.backend",
                "Only the ddx PCM backend has an implemented runtime energy path.",
                value=backend,
                expected="ddx",
                action="Use backend=ddx; the pcmsolver runtime coupling is not implemented yet.",
                wiki=WIKI_HELP["pcm.backend"],
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


def _dftb_key_customized(config: dict[str, Any], key: str) -> bool:
    """True when a [dftb] key differs from its schema default (i.e. the user
    tuned it). The parsed config always carries defaults, so presence alone
    cannot distinguish user intent."""
    from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
    spec = OQP_CONFIG_SCHEMA.get("dftb", {}).get(key)
    if spec is None:
        return False
    value = _get(config, "dftb", key, None)
    if value is None:
        return False
    default_raw = spec.get("default", "")
    try:
        return abs(float(value) - float(default_raw)) > 1e-300
    except (TypeError, ValueError):
        return str(value).strip().lower() != str(default_raw).strip().lower()


def apply_dftb_model_default(config: dict[str, Any]) -> str:
    """Materialize the production default [dftb] model preset.

    MRSF-TDDFTB defaults to the tuned DTCAM-TB operator (``model=dtcam``); the
    other open-shell SCC routes (SF-TDDFTB, open-shell ground) default to the
    conventional OB2 / LC-DFTB2 protocol (``model=ob2``, SPC=0.5).  The OB2
    protocol can be requested for MRSF with [dftb] model=ob2.  Closed-shell
    routes (singlet ground, TD-DFTB)
    and DFTB0 stay preset-free because the restricted reference has no
    long-range exchange yet.  ``model=none`` keeps the explicit-keys route,
    and any tuned preset-locked key implies manual operator control (legacy
    inputs keep meaning what they said).  Returns the effective model string.
    """
    if _as_lower(_get(config, "input", "method", "hf")) != "dftb":
        return ""
    dftb = config.get("dftb")
    if not isinstance(dftb, dict):
        return ""
    model = _as_lower(_get(config, "dftb", "model", ""))
    if model == "none":
        dftb["model"] = ""
        return ""
    if model:
        return model
    # Presets are resolved inside the native library; the probe backend
    # cannot carry them, so it keeps the explicit-keys route.
    if _as_lower(_get(config, "dftb", "backend", "native")) not in {
            "native", "auto"}:
        return ""
    from oqp.utils.state_labels import resolved_dftb_type  # noqa: PLC0415
    route = resolved_dftb_type(config)
    if any(_dftb_key_customized(config, key)
           for key in DFTB_MODEL_LOCKED_KEYS):
        return ""
    if route == "mrsf":
        dftb["model"] = DFTB_MODEL_DEFAULT_MRSF
    elif route == "sf" or (
            route == "ground"
            and _get(config, "dftb", "reference_multiplicity", 0)
            and int(_get(config, "dftb", "reference_multiplicity", 0)) > 1):
        dftb["model"] = DFTB_MODEL_DEFAULT_OPEN_SHELL
    else:
        # Closed-shell routes (singlet ground, TD-DFTB) and DFTB0: the
        # restricted reference has no LC yet, so no preset default.
        return ""
    return dftb["model"]


def _check_dftb(config: dict[str, Any], report: CheckReport) -> None:
    _check_tb(config, report, section="dftb")


def _check_xtb(config: dict[str, Any], report: CheckReport) -> None:
    _check_tb(config, report, section="xtb")


def _check_tb(config: dict[str, Any], report: CheckReport, *, section: str) -> None:
    """Validate the [dftb]/[xtb] TB backend section (shared option family)."""
    method = _as_lower(_get(config, "input", "method", "hf"))
    if method != section:
        return
    is_xtb = section == "xtb"
    disp = "OpenQP-xTB" if is_xtb else "OpenQP-DFTB"     # display name
    short = "xTB" if is_xtb else "DFTB"                  # terse family name
    lib_stem = "libopenqp_xtb_c" if is_xtb else "libopenqp_dftb_c"

    backend = _as_lower(_get(config, section, "backend", "native"))
    dftb_type = _as_lower(_get(config, section, "type", "auto"))
    # Canonicalize response-type aliases so the runtype/NAMD/SOC gates below see
    # a single spelling per family (e.g. mrsftddftb/mrsf-tddftb -> mrsf).
    _DFTB_TYPE_CANON = {
        "dftb": "ground", "dftb0": "ground_noscc", "noscc": "ground_noscc",
        "td-dftb": "tddftb", "tda": "tddftb",
        "sftddftb": "sf", "sf-tddftb": "sf",
        "mrsftddftb": "mrsf", "mrsf-tddftb": "mrsf",
    }
    dftb_type_canon = _DFTB_TYPE_CANON.get(dftb_type, dftb_type)
    parameter_path = _get(config, section, "parameter_path", "")
    executable = _get(config, section, "executable", "")
    runtype = _as_lower(_get(config, "input", "runtype", "energy"))
    td_type = _as_lower(_get(config, "tdhf", "type", "rpa"))
    nstate = _get(config, "tdhf", "nstate", 1)
    grad_states = _as_list(_get(config, "properties", "grad", []))
    istate = _get(config, "optimize", "istate", 0)
    jstate = _get(config, "optimize", "jstate", 0)

    allowed_backends = {"auto", "native"} if is_xtb else DFTB_BACKENDS
    if backend not in allowed_backends:
        message = (f"{disp} does not support the probe backend (DFTB-only "
                   "developer fallback)." if is_xtb and backend == "probe"
                   else f"Unknown {disp} backend.")
        action_tail = ("." if is_xtb
                       else "; probe is only an explicit developer fallback.")
        report.add(
            "ERROR",
            f"{section}.backend",
            message,
            value=backend,
            expected=", ".join(sorted(allowed_backends)),
            action=f"Use native (loads the standalone {lib_stem} shared library){action_tail}",
        )

    if backend == "probe" and not is_xtb and runtype in {"nacme", "nac", "namd", "soc"}:
        report.add(
            "ERROR",
            f"{section}.backend",
            f"{disp} nacme/nac/namd/soc require the native backend.",
            value=backend,
            expected="native",
            action=f"Use [{section}] backend=native: the probe executable only returns "
                   "energies/gradients and cannot publish the MO/response-vector "
                   "tags or the state-overlap/SOC entry points these workflows need.",
        )

    if backend == "probe" and not is_xtb:
        # The probe executable receives only a fixed subset of scalar controls;
        # options it cannot forward would be silently ignored, so reject them
        # rather than run with probe defaults.
        _probe_unforwarded = []
        if int(_get(config, section, "target_multiplicity", 1)) != 1:
            _probe_unforwarded.append("target_multiplicity")
        for _ch in ("spc_coco", "spc_ovov", "spc_coov"):
            if float(_get(config, section, _ch, -999.0)) != -999.0:
                _probe_unforwarded.append(_ch)
        for _sh in ("mrsf_shift_oo", "mrsf_shift_co", "mrsf_shift_ov", "mrsf_shift_cv"):
            if float(_get(config, section, _sh, 0.0)) != 0.0:
                _probe_unforwarded.append(_sh)
        _lc = _get(config, section, "lc_ground_state", False)
        if (_lc is True) or (str(_lc).lower() in ("true", "1", "on", "yes")):
            _probe_unforwarded.append("lc_ground_state")
        _zv = _get(config, section, "zvector", True)
        if (_zv is False) or (str(_zv).lower() in ("false", "0", "off", "no")):
            _probe_unforwarded.append("zvector")
        # Response-solver / reference controls the probe CLI also never forwards.
        if float(_get(config, section, "response_tolerance", 1.0e-6)) != 1.0e-6:
            _probe_unforwarded.append("response_tolerance")
        if int(_get(config, section, "response_max_iterations", 50)) != 50:
            _probe_unforwarded.append("response_max_iterations")
        if int(_get(config, section, "response_max_subspace", 100)) != 100:
            _probe_unforwarded.append("response_max_subspace")
        if str(_get(config, section, "response_solver", "auto")).lower() != "auto":
            _probe_unforwarded.append("response_solver")
        if int(_get(config, section, "reference_multiplicity", 0)) != 0:
            _probe_unforwarded.append("reference_multiplicity")
        _sc = _get(config, section, "spin_complete", True)
        if (_sc is False) or (str(_sc).lower() in ("false", "0", "off", "no")):
            _probe_unforwarded.append("spin_complete")
        for _op in DFTB_NATIVE_ONLY_KEYS:
            if _dftb_key_customized(config, _op):
                _probe_unforwarded.append(_op)
        if _probe_unforwarded:
            report.add(
                "ERROR",
                f"{section}.backend",
                f"{disp} backend=probe cannot forward these options to the "
                "probe executable; they would be silently ignored.",
                value=", ".join(_probe_unforwarded),
                expected="native",
                action=f"Use [{section}] backend=native for target_multiplicity, "
                       "per-channel spc_*, mrsf_shift_*, lc_ground_state, "
                       "zvector=false, or the DTCAM operator keys.",
            )

    model = _as_lower(_get(config, "dftb", "model", ""))
    if model == "none":
        # Explicit opt-out of the SF/MRSF dtcam default: explicit-keys route.
        model = ""
    if model:
        if model not in DFTB_MODELS:
            report.add(
                "ERROR",
                "dftb.model",
                "Unknown OpenQP-DFTB operator model preset.",
                value=model,
                expected=", ".join(sorted(DFTB_MODELS)),
                action="Use model=dtcam (DTCAM-TB paper vector), "
                       "model=ob2 (conventional OB2/LC-DFTB2 protocol), or "
                       "omit model and set the operator keys individually.",
            )
        if backend == "probe":
            report.add(
                "ERROR",
                "dftb.model",
                "Operator model presets require the native backend.",
                value=backend,
                expected="native",
                action="Use backend=native (presets are resolved inside openqp-dftb).",
            )
        # LC-DFTB2 exists only for the spin-polarized (ROKS) reference; the
        # native library rejects lc_ground_state on the closed-shell
        # restricted path, so surface that at input time for LC presets.
        if model in {"ob2", "dftb+", "dftbplus", "dftb_plus"} and \
                dftb_type_canon in {"ground", "ground_noscc", "tddftb"}:
            reference_multiplicity = _get(
                config, "dftb", "reference_multiplicity", 0)
            try:
                reference_multiplicity = int(reference_multiplicity)
            except (TypeError, ValueError):
                reference_multiplicity = 0
            if reference_multiplicity <= 1:
                report.add(
                    "ERROR",
                    "dftb.model",
                    "model=ob2 sets an LC-DFTB2 (lc_ground_state) "
                    "reference, which is only implemented for the "
                    "spin-polarized (ROKS) path.",
                    value=f"{model} with closed-shell type={dftb_type_canon}",
                    expected="an SF/MRSF route or reference_multiplicity > 1",
                    action="Use type=sf/mrsf, set reference_multiplicity > 1, "
                           "or drop model= (closed-shell runs stay preset-free).",
                )
        conflicting = [
            key for key in DFTB_MODEL_LOCKED_KEYS
            if _dftb_key_customized(config, key)
        ]
        if conflicting:
            report.add(
                "ERROR",
                "dftb.model",
                "A model preset fixes the operator and numerical protocol; "
                "only the explicit SCC mixer override is allowed.",
                value=", ".join(conflicting),
                expected=f"model={model} with no tuned operator keys",
                action="Remove the listed [dftb] keys or drop model= to tune manually.",
            )

    scf_prop = _as_list(_get(config, "properties", "scf_prop", []))
    if scf_prop:
        report.add(
            "ERROR",
            "properties.scf_prop",
            f"{disp} does not implement SCF property analyses "
            "(mulliken/lowdin/nmr/...); they dispatch to Gaussian-basis "
            f"routines the {short} adapter never initializes.",
            value=", ".join(str(v) for v in scf_prop),
            expected=f"omit properties.scf_prop for method={section}",
            action="Remove properties.scf_prop, or use an all-electron method for SCF properties.",
        )

    if bool(_get(config, "pcm", "enabled", False)):
        report.add(
            "ERROR",
            "pcm.enabled",
            f"{disp} does not consume a PCM solvent reaction field; the {short} "
            "energy/gradient path returns directly from the adapter and never "
            "calls ddX/PCM, so enabling PCM would silently produce gas-phase "
            f"{short} results instead of solvated ones.",
            value=True,
            expected=f"pcm.enabled=false for method={section}",
            action=f"Disable [pcm] enabled for method={section}, or use an all-electron "
                   "method for PCM solvation.",
            wiki=WIKI_HELP["pcm.enabled"],
        )

    if dftb_type not in DFTB_TYPES:
        report.add(
            "ERROR",
            f"{section}.type",
            f"Unknown {disp} calculation type.",
            value=dftb_type,
            expected="ground, tddftb, sf, mrsf, or auto",
            action="Use auto to derive the DFTB response type from [tdhf] type, or set it explicitly.",
        )

    explicit_response_types = {
        "tddftb": {"rpa", "tda"},
        "sf": {"sf"},
        "mrsf": {"mrsf"},
    }
    compatible_tdhf_types = explicit_response_types.get(dftb_type_canon)
    if (
        dftb_type != "auto"
        and compatible_tdhf_types is not None
        and td_type not in compatible_tdhf_types
    ):
        report.add(
            "ERROR",
            "dftb.type",
            "Explicit OpenQP-DFTB response type conflicts with [tdhf] type.",
            value={"dftb.type": dftb_type, "tdhf.type": td_type},
            expected="tdhf.type=" + "/".join(sorted(compatible_tdhf_types)),
            action="Align [dftb] type with [tdhf] type, or use [dftb] type=auto.",
        )

    auto_ground = (
        dftb_type == "auto"
        and (
            (runtype in {"energy", "grad"} and not any(int(state) > 0 for state in grad_states))
            or (runtype in {"optimize", "mep"} and int(istate) == 0)
        )
    )
    conventional_tddftb = dftb_type_canon == "tddftb" or (
        dftb_type == "auto" and td_type in {"rpa", "tda"} and not auto_ground
    )
    if conventional_tddftb:
        target_mult = int(_get(config, "dftb", "target_multiplicity", 1))
        reference_mult = int(_get(config, "dftb", "reference_multiplicity", 0))
        td_mult = int(_get(config, "tdhf", "multiplicity", 1))
        if target_mult != 1 or td_mult != 1 or reference_mult not in {0, 1}:
            report.add(
                "ERROR",
                "dftb.target_multiplicity",
                "Conventional TD-DFTB currently implements a closed-shell singlet "
                "reference and singlet TDA response only.",
                value={
                    "reference": reference_mult or 1,
                    "dftb_target": target_mult,
                    "tdhf_target": td_mult,
                },
                expected="reference and target multiplicity 1",
                action="Use SF-TDDFTB or MRSF-TDDFTB for triplet-state calculations.",
            )

    lc_gamma = _as_lower(_get(config, section, "lc_gamma", "ok" if is_xtb else "yukawa"))
    lc_gamma_allowed = {"yukawa", "erf", "ok"} if is_xtb else {"yukawa", "erf"}
    if lc_gamma not in lc_gamma_allowed:
        report.add(
            "ERROR",
            f"{section}.lc_gamma",
            f"Unknown {disp} long-range gamma kernel.",
            value=lc_gamma,
            expected="yukawa, erf, or ok" if is_xtb else "yukawa or erf",
            action=("Use ok (Ohno-Klopman, the GFN1 default), yukawa, or erf (erf(omega R)/R)."
                    if is_xtb else
                    "Use yukawa (DFTB+ LC-DFTB2 Yukawa-Slater gamma) or erf (erf(omega R)/R)."),
        )

    response_solver = _as_lower(_get(config, section, "response_solver", "auto"))
    if response_solver not in {"auto", "dense", "davidson"}:
        report.add(
            "ERROR",
            f"{section}.response_solver",
            f"Unknown {disp} response solver.",
            value=response_solver,
            expected="auto, dense, or davidson",
            action="Use auto unless you need to force a specific eigensolver.",
        )

    print_level = int(_get(config, "dftb", "print_level", 1))
    if print_level not in {0, 1, 2}:
        report.add(
            "ERROR",
            "dftb.print_level",
            "Unknown OpenQP-DFTB native trace level.",
            value=print_level,
            expected="0, 1, or 2",
            action="Use 0 for quiet, 1 for stage summaries, or 2 for iteration detail.",
        )

    target_multiplicity = int(_get(config, section, "target_multiplicity", 1))
    if target_multiplicity not in (1, 3):
        report.add(
            "ERROR",
            f"{section}.target_multiplicity",
            f"{disp} MRSF/SF response states must be singlet or triplet.",
            value=target_multiplicity,
            expected="1 or 3",
            action="Use 1 for singlet response states or 3 for triplet response states.",
        )

    param_env = "OPENQP_XTB_PARAMETER_PATH" if is_xtb else "OPENQP_DFTB_PARAMETER_PATH"
    if not parameter_path and not os.environ.get(param_env):
        if is_xtb:
            # openqp-xtb ships no bundled default parameter set; an explicit
            # converter-generated .opxtb source is always required.
            report.add(
                "ERROR",
                f"{section}.parameter_path",
                f"{disp} requires an explicit parameter source.",
                expected=".opxtb parameter file (converter-generated)",
                action=f"Set [{section}] parameter_path, or export {param_env}.",
                wiki=WIKI_HELP[f"{section}.parameter_path"],
            )
        else:
            # The openqp-dftb wheel ships a bundled default parameter set
            # (OB2W0PT3 with the official spinw.txt); when it is resolvable the
            # runtime uses it automatically and no explicit source is required.
            bundled = None
            try:
                import openqp_dftb  # noqa: PLC0415

                bundled = openqp_dftb.default_parameter_path()
            except (ImportError, AttributeError, FileNotFoundError):
                bundled = None
            if not bundled:
                report.add(
                    "ERROR",
                    "dftb.parameter_path",
                    "OpenQP-DFTB could not resolve a parameter source.",
                    expected="bundled parameters, an .opdftb file, or an SKF directory",
                    action=(
                        "Install an openqp-dftb wheel with bundled parameters, set "
                        "[dftb] parameter_path, or export OPENQP_DFTB_PARAMETER_PATH."
                    ),
                    wiki=WIKI_HELP["dftb.parameter_path"],
                )

    if backend == "probe" and not is_xtb and not executable and not os.environ.get("OPENQP_DFTB_STATE_GRADIENT_PROBE"):
        report.add(
            "WARNING",
            f"{section}.executable",
            "Probe backend selected without an explicit probe executable.",
            action="Set [dftb] executable or put openqp_dftb_state_gradient_probe on PATH.",
        )

    for path, value in (
        (f"{section}.scc_tolerance", _get(config, section, "scc_tolerance", 1.0e-8)),
        (f"{section}.response_tolerance", _get(config, section, "response_tolerance", 1.0e-6)),
        (f"{section}.scc_mixing", _get(config, section, "scc_mixing", 0.35)),
    ):
        if value <= 0.0:
            report.add(
                "ERROR",
                path,
                f"{short} numeric option must be positive.",
                value=value,
                action="Use a positive value.",
            )

    for path, value in (
        (f"{section}.omega", _get(config, section, "omega", 0.3)),
        (f"{section}.cam_alpha", _get(config, section, "cam_alpha", 0.0)),
        (f"{section}.cam_beta", _get(config, section, "cam_beta", 1.0)),
    ):
        if value < 0.0:
            report.add(
                "ERROR",
                path,
                f"{short} range-separation/CAM option must be non-negative.",
                value=value,
                action="Use zero or a positive value.",
            )

    if _get(config, section, "scc_mixing", 0.35) > 1.0:
        report.add(
            "ERROR",
            f"{section}.scc_mixing",
            "SCC mixing must be in the interval (0, 1].",
            value=_get(config, section, "scc_mixing", 0.35),
            action="Reduce the population-vector mixing fraction.",
        )

    for path in ("scc_history", "max_scc_iterations", "response_max_iterations", "response_max_subspace", "timeout"):
        value = _get(config, section, path, 1)
        if value <= 0:
            report.add(
                "ERROR",
                f"{section}.{path}",
                f"{short} integer option must be positive.",
                value=value,
                action="Use a positive integer.",
            )

    if _get(config, section, "scc_max_step", 0.5) < 0.0:
        report.add(
            "ERROR",
            f"{section}.scc_max_step",
            "SCC mixer max step must be non-negative.",
            value=_get(config, section, "scc_max_step", 0.5),
            action="Use zero to disable the cap or a positive cap.",
        )

    if _get(config, section, "spc", 0.5) < 0.0 and abs(_get(config, section, "spc", 0.5) + 1.0) > 1.0e-12:
        report.add(
            "ERROR",
            f"{section}.spc",
            "MRSF spin-pair coupling must be non-negative, or -1 to inherit the CAM exchange fraction.",
            value=_get(config, section, "spc", 0.5),
            action=f"Use 0.5 for the current {disp} production calibration.",
        )

    if _as_lower(_get(config, section, "scc_mixer", "auto")) not in DFTB_SCC_MIXERS:
        report.add(
            "ERROR",
            f"{section}.scc_mixer",
            f"Unknown {short} SCC mixer.",
            value=_get(config, section, "scc_mixer", "auto"),
            expected=", ".join(sorted(DFTB_SCC_MIXERS)),
            action="Use auto, anderson, broyden, diis, trust/trah, or linear.",
        )

    # QM/MM (ESPF/OpenMM) with the DFTB backend: the openqp-dftb library takes
    # the MM electrostatic potential directly into the SCC Hamiltonian
    # (mol.dftb_external_potential), so the OpenQpQMMM/QMMM_MD drivers support
    # method=dftb with full-ESPF electrostatic or mechanical embedding
    # (including hydrogen link atoms across covalent QM/MM cuts). The legacy
    # 'split' scheme (QM point charges routed through OpenMM) is rejected: it
    # would double-count the coupling already inside the embedded DFTB energy.
    # runtype=namd with method=dftb stays blocked by the generic NAMD gate
    # until a DFTB state-overlap/NACME backend exists.
    if bool(_get(config, "input", "qmmm_flag", False)):
        embedding = str(
            _get(config, "qmmm", "embedding", "electrostatic") or "electrostatic"
        ).strip().lower()
        if embedding == "split":
            report.add(
                "ERROR",
                "qmmm.embedding",
                f"{disp} QM/MM does not support the legacy 'split' embedding.",
                value=embedding,
                expected="electrostatic, espf, or mechanical",
                action="Use [qmmm] embedding=electrostatic (full-ESPF scheme) "
                       f"or mechanical with method={section}.",
            )
        if runtype == "namd":
            soc_val = _get(config, "md", "soc", False)
            soc_on = (soc_val is True) or (str(soc_val).lower() in ("true", "1", "on", "yes"))
            if soc_on:
                report.add(
                    "ERROR",
                    "md.soc",
                    f"{disp} SOC-NAMD with QM/MM embedding is not yet wired: "
                    "the SOC electronic/gradient path still uses the native ESPF "
                    f"hcore/gradient hooks, which the {short} backend does not consume.",
                    value="soc=true with qmmm_flag=true",
                    expected=f"gas-phase SOC-NAMD (qmmm_flag=false), or non-SOC {short} QM/MM NAMD",
                    action=f"Disable [md] soc for {short} QM/MM NAMD, or run SOC-NAMD without QM/MM.",
                )
            if embedding == "mechanical":
                report.add(
                    "ERROR",
                    "qmmm.embedding",
                    f"{disp} QM/MM NAMD requires full-ESPF electrostatic embedding; "
                    "the mechanical path raises at the first electronic step.",
                    value=embedding,
                    expected="electrostatic",
                    action=f"Use [qmmm] embedding=electrostatic for {short} QM/MM NAMD.",
                )

    allowed_runtype = {"energy", "grad", "optimize", "meci", "mep", "data"}
    td_type_for_namd = _as_lower(_get(config, "tdhf", "type", "rpa"))
    if td_type_for_namd == "mrsf" and dftb_type_canon in {"auto", "mrsf"}:
        # MRSF-TDDFTB has a state-overlap (TLF) backend and one-center SOC:
        # overlap-based couplings and surface hopping are available.
        allowed_runtype |= {"nac", "nacme", "namd", "soc"}
    if runtype in {"nac", "nacme", "bp"} and runtype not in allowed_runtype:
        report.add(
            "ERROR",
            "input.runtype",
            f"{disp} NAMD coupling requires the MRSF state-overlap backend of the {short} response.",
            value=runtype,
            expected=f"tdhf.type=mrsf (and {section}.type auto/mrsf) for nac/nacme; bp is not available",
            action="Use tdhf.type=mrsf for DFTB nonadiabatic couplings.",
            wiki=WIKI_HELP[f"{section}.namd"],
        )
    elif runtype not in allowed_runtype:
        report.add(
            "ERROR",
            "input.runtype",
            f"{disp} is currently wired only through energy/gradient-driven workflows.",
            value=runtype,
            expected=", ".join(sorted(allowed_runtype)),
            action=f"Use energy, grad, data, optimize, meci, or mep until Hessian/NAC/SOC {short} hooks are implemented.",
        )

    # [properties] nac also routes into the NAC state-overlap path (e.g. inside
    # runtype=data), not only the standalone nac/nacme runtypes gated above.
    props_nac = _as_list(_get(config, "properties", "nac", []))
    if props_nac and "nac" not in allowed_runtype and runtype not in {"nac", "nacme", "bp"}:
        report.add(
            "ERROR",
            "properties.nac",
            "OpenQP-DFTB NAC requested via [properties] nac requires the "
            "MRSF-TDDFTB state-overlap backend.",
            value=", ".join(str(v) for v in props_nac),
            expected="tdhf.type=mrsf (and dftb.type auto/mrsf)",
            action="Use tdhf.type=mrsf for DFTB nonadiabatic couplings, or drop properties.nac.",
        )

    ground_like = dftb_type in {"ground", "dftb", "dftb0", "ground_noscc", "noscc"}
    if dftb_type == "auto" and runtype in {"optimize", "mep"} and int(istate) == 0:
        ground_like = True

    if not ground_like and td_type not in {"rpa", "tda", "sf", "mrsf"}:
        report.add(
            "ERROR",
            "tdhf.type",
            ("OpenQP-xTB response hooks support TD/TDA, SF, and MRSF responses."
             if is_xtb else
             "OpenQP-DFTB response hooks support TDDFTB/TDA, SF-TDDFTB, and MRSF-TDDFTB."),
            value=td_type,
            expected="rpa/tda, sf, or mrsf",
            action=f"Set [tdhf] type=tda, sf, or mrsf, or set [{section}] type=ground for S0 {short}.",
        )

    if not ground_like:
        try:
            nstate_val = int(nstate)
        except (TypeError, ValueError):
            nstate_val = 0
        if nstate_val < 1:
            report.add(
                "ERROR",
                "tdhf.nstate",
                f"{disp} excited-state response needs at least one root.",
                value=nstate,
                expected="tdhf.nstate >= 1",
                action="Set [tdhf] nstate to the number of response roots (>= 1).",
            )

    if is_xtb:
        model = _as_lower(_get(config, section, "model", "gfn1"))
        if model not in {"gfn1"}:
            report.add(
                "ERROR",
                "xtb.model",
                "Unknown OpenQP-xTB model.",
                value=model,
                expected="gfn1",
                action="Use model=gfn1 (LC-GFN1-xTB); the DFTB2 model lives in the [dftb] backend.",
            )
        spin_scale = _get(config, section, "spin_scale", 1.0)
        if spin_scale <= 0.0:
            report.add(
                "ERROR",
                "xtb.spin_scale",
                "xTB spin-polarization scale must be positive.",
                value=spin_scale,
                action="Use spin_scale > 0 (1.0 = unscaled spGFN1 W constants).",
            )

    requested = []
    if runtype in {"grad", "data"}:
        requested.extend(grad_states)
    if runtype in {"optimize", "mep"}:
        requested.append(istate)
    if runtype == "meci":
        requested.extend([istate, jstate])
    max_requested = _max_state(requested)

    if ground_like and max_requested > 0:
        report.add(
            "ERROR",
            f"{section}.type",
            f"Ground-state {short} can only target state 0.",
            value={"type": dftb_type, "requested_state": max_requested},
            expected="state 0",
            action=f"Use [{section}] type=tddftb/sf/mrsf for excited states.",
        )

    if not ground_like and max_requested > nstate:
        report.add(
            "ERROR",
            "tdhf.nstate",
            f"{short} requested state is above the solved response space.",
            value={"nstate": nstate, "requested_state": max_requested},
            expected=f"nstate >= {max_requested}",
            action="Increase [tdhf] nstate.",
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

    # scf.escalation: optional comma-separated ladder overriding the default
    # DIIS -> SOSCF -> TRAH chain. Each stage must be a concrete converger
    # (not the 'auto'/'ml' manager modes).
    escalation = _as_lower(_get(config, "scf", "escalation", ""))
    escalation_stages = {"diis", "soscf", "trah"}
    for stage in (s.strip() for s in escalation.split(",") if s.strip()):
        if stage not in escalation_stages:
            report.add(
                "ERROR",
                "scf.escalation",
                "Unknown SCF escalation stage.",
                value=stage,
                expected=", ".join(sorted(escalation_stages)),
                action="List concrete convergers, e.g. soscf,trah.",
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


def _check_mp2(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    if method != "mp2":
        return

    functional = _get(config, "input", "functional", "")
    variant = _as_lower(_get(config, "mp2", "variant", "mp2"))
    ss_scale = _get(config, "mp2", "same_spin_scale", 1.0)
    os_scale = _get(config, "mp2", "opposite_spin_scale", 1.0)

    if functional:
        report.add(
            "ERROR",
            "input.functional",
            "Standalone MP2 requires an HF reference, not a DFT functional.",
            value=functional,
            expected="empty functional",
            action="Remove [input] functional for method=mp2.",
            wiki=WIKI_HELP["input.method"],
        )

    if variant not in MP2_VARIANTS:
        report.add(
            "ERROR",
            "mp2.variant",
            "Unknown MP2 spin-scaling variant.",
            value=variant,
            expected=", ".join(sorted(MP2_VARIANTS)),
            action="Choose a named variant or use variant=custom with explicit scales.",
            wiki=WIKI_HELP["mp2.variant"],
        )

    if variant == "custom":
        try:
            finite_scales = math.isfinite(float(ss_scale)) and math.isfinite(float(os_scale))
        except (TypeError, ValueError):
            finite_scales = False
        if not finite_scales:
            report.add(
                "ERROR",
                "mp2.same_spin_scale/mp2.opposite_spin_scale",
                "Custom MP2 scale factors must be finite numbers.",
                value=f"{ss_scale}/{os_scale}",
                action="Set finite same_spin_scale and opposite_spin_scale values.",
                wiki=WIKI_HELP["mp2.variant"],
            )


def _check_fci(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    if method != "fci":
        return

    runtype = _as_lower(_get(config, "input", "runtype", "energy"))
    scf_type = _as_lower(_get(config, "scf", "type", "rhf"))
    multiplicity = _get(config, "scf", "multiplicity", 1)
    functional = _get(config, "input", "functional", "")
    d4_enabled = _get(config, "input", "d4", False)
    nroot = _get(config, "fci", "nroot", 1)
    active_electrons = _get(config, "fci", "active_electrons", 0)
    active_orbitals = _get(config, "fci", "active_orbitals", 0)
    frozen_core = _get(config, "fci", "frozen_core", 0)
    max_det = _get(config, "fci", "max_det", 5000)
    max_memory = _get(config, "fci", "max_memory", 2048)
    eig_tol = _get(config, "fci", "eig_tol", 1.0e-10)
    integral_backend = _check_choice_literal(
        _get(config, "fci", "integral_backend", "native"),
        "fci.integral_backend",
        FCI_INTEGRAL_BACKENDS,
        report,
        message="Unknown FCI integral backend.",
        action="Use [fci] integral_backend=native.",
        fallback="native",
    )
    integral_cutoff = _get(config, "fci", "integral_cutoff", 5.0e-11)
    solver = _check_choice_literal(
        _get(config, "fci", "solver", "auto"),
        "fci.solver",
        FCI_SOLVERS,
        report,
        message="Unknown FCI solver.",
        action="Use [fci] solver=auto, dense, or davidson.",
        fallback="auto",
    )
    davidson_maxiter = _get(config, "fci", "davidson_maxiter", 100)
    davidson_subspace = _get(config, "fci", "davidson_subspace", 0)
    _validate_bool_literal(
        _get(config, "fci", "print_ci_vectors", False),
        "fci.print_ci_vectors",
        report,
        action="Set [fci] print_ci_vectors=true or false.",
    )
    _validate_bool_literal(
        _get(config, "fci", "save_ci_vectors", False),
        "fci.save_ci_vectors",
        report,
        action="Set [fci] save_ci_vectors=true or false.",
    )
    _validate_bool_literal(
        _get(config, "fci", "save_rdm", False),
        "fci.save_rdm",
        report,
        action="Set [fci] save_rdm=true or false.",
    )
    ci_print_threshold = _get(config, "fci", "ci_print_threshold", 5.0e-2)
    target_spin = _as_lower(_get(config, "fci", "target_spin", "any"))

    nroot_valid, nroot_int = _check_int_literal(
        nroot,
        "fci.nroot",
        report,
        minimum=1,
        expected="positive integer",
        action="Set [fci] nroot=1 or higher.",
    )
    if not nroot_valid:
        nroot_int = 1

    active_electrons_valid, active_electrons_total = _check_active_electron_literal(
        active_electrons,
        "fci.active_electrons",
        report,
    )
    active_count_values: dict[str, int] = {
        "active_electrons": active_electrons_total if active_electrons_valid else 0,
    }
    for option, value in (
        ("active_orbitals", active_orbitals),
        ("frozen_core", frozen_core),
    ):
        valid_count, parsed_count = _check_int_literal(
            value,
            f"fci.{option}",
            report,
            expected="non-negative integer",
            action="Use 0 for the default full-space setting or a positive count.",
        )
        active_count_values[option] = parsed_count if valid_count else 0

    _check_int_literal(
        max_det,
        "fci.max_det",
        report,
        minimum=1,
        expected="positive integer",
        action="Increase [fci] max_det.",
    )
    _check_int_literal(
        max_memory,
        "fci.max_memory",
        report,
        minimum=1,
        expected="positive integer MiB",
        action="Set [fci] max_memory to a positive size in MiB.",
    )
    eig_tol_valid, eig_tol_float = _check_float_literal(
        eig_tol,
        "fci.eig_tol",
        report,
        expected="positive finite number",
        action="Set [fci] eig_tol to a positive threshold.",
    )
    integral_cutoff_valid, integral_cutoff_float = _check_float_literal(
        integral_cutoff,
        "fci.integral_cutoff",
        report,
        expected="finite non-negative number",
        action="Set [fci] integral_cutoff to zero or a positive threshold.",
    )
    _check_int_literal(
        davidson_maxiter,
        "fci.davidson_maxiter",
        report,
        minimum=1,
        expected="positive integer",
        action="Set [fci] davidson_maxiter to a positive integer.",
    )
    davidson_subspace_valid, davidson_subspace_int = _check_int_literal(
        davidson_subspace,
        "fci.davidson_subspace",
        report,
        minimum=0,
        expected="0 or positive integer",
        action="Use 0 for the internal default.",
    )
    ci_print_threshold_valid, ci_print_threshold_float = _check_float_literal(
        ci_print_threshold,
        "fci.ci_print_threshold",
        report,
        expected="finite non-negative number",
        action="Set [fci] ci_print_threshold to zero or a positive value.",
    )

    if runtype != "energy":
        report.add(
            "ERROR",
            "input.runtype",
            "FCI is currently implemented only for energy calculations.",
            value=runtype,
            expected="energy",
            action="Set [input] runtype=energy.",
        )

    if scf_type != "rhf" or multiplicity != 1:
        report.add(
            "ERROR",
            "scf.type",
            "FCI MVP requires a closed-shell RHF reference.",
            value=f"{scf_type}/{multiplicity}",
            expected="scf.type=rhf and scf.multiplicity=1",
            action="Use an RHF singlet reference or extend the FCI implementation.",
        )

    if functional:
        report.add(
            "ERROR",
            "input.functional",
            "FCI requires HF one- and two-electron integrals.",
            value=functional,
            action="Unset [input] functional for method=fci.",
        )

    if d4_enabled:
        report.add(
            "ERROR",
            "input.d4",
            "FCI reports the variational electronic energy; D4 dispersion correction is not part of the FCI solver.",
            value=d4_enabled,
            expected=False,
            action="Disable [input] d4 for method=fci.",
        )

    if min(active_count_values.values()) < 0:
        report.add(
            "ERROR",
            "fci.active_space",
            "FCI active-space counts cannot be negative.",
            value={
                "active_electrons": active_electrons,
                "active_orbitals": active_orbitals,
                "frozen_core": frozen_core,
            },
            expected="non-negative integers",
            action="Use 0 for the default full-space setting or a positive count.",
        )

    if eig_tol_valid and eig_tol_float <= 0:
        report.add(
            "ERROR",
            "fci.eig_tol",
            "FCI eigensolver tolerance must be positive.",
            value=eig_tol,
            expected="> 0",
            action="Set [fci] eig_tol to a positive threshold.",
        )

    if integral_cutoff_valid and integral_cutoff_float < 0:
        report.add(
            "ERROR",
            "fci.integral_cutoff",
            "FCI integral cutoff cannot be negative.",
            value=integral_cutoff,
            expected=">= 0",
            action="Set [fci] integral_cutoff to zero or a positive threshold.",
        )

    if integral_backend not in FCI_INTEGRAL_BACKENDS:
        report.add(
            "ERROR",
            "fci.integral_backend",
            "Unknown FCI integral backend.",
            value=integral_backend,
            expected=", ".join(sorted(FCI_INTEGRAL_BACKENDS)),
            action="Use [fci] integral_backend=native.",
        )

    if solver not in FCI_SOLVERS:
        report.add(
            "ERROR",
            "fci.solver",
            "Unknown FCI solver.",
            value=solver,
            expected=", ".join(sorted(FCI_SOLVERS)),
            action="Use [fci] solver=auto, dense, or davidson.",
        )

    if (
        nroot_valid
        and davidson_subspace_valid
        and davidson_subspace_int > 0
        and davidson_subspace_int < nroot_int
    ):
        report.add(
            "ERROR",
            "fci.davidson_subspace",
            "FCI Davidson subspace override must be at least the requested root count.",
            value=davidson_subspace,
            expected=f"0 or >= nroot ({nroot_int})",
            action="Use 0 for the internal default or raise [fci] davidson_subspace.",
        )

    if ci_print_threshold_valid and ci_print_threshold_float < 0:
        report.add(
            "ERROR",
            "fci.ci_print_threshold",
            "FCI CI print threshold cannot be negative.",
            value=ci_print_threshold,
            expected=">= 0",
            action="Set [fci] ci_print_threshold to zero or a positive value.",
        )

    if not _target_spin_is_valid(target_spin):
        report.add(
            "ERROR",
            "fci.target_spin",
            "Unknown target spin label.",
            value=target_spin,
            expected="any, a named spin multiplicity, or a positive integer multiplicity",
            action="Use [fci] target_spin=any, singlet, doublet, triplet, quartet, or a positive integer multiplicity.",
        )


def _target_spin_is_valid(target_spin: Any) -> bool:
    if target_spin is None:
        return True
    if isinstance(target_spin, str):
        target = target_spin.strip().lower()
        if not target or target == "any" or target in TARGET_SPIN_LABELS:
            return True
        valid, parsed = _parse_int_literal(target)
        return valid and parsed >= 1
    valid, parsed = _parse_int_literal(target_spin)
    return valid and parsed >= 1


def _check_casci(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    if method not in {
        "casci", "casscf", "sa-casscf", "sacasscf",
        "caspt2", "ms-caspt2", "mscaspt2", "xms-caspt2", "xmscaspt2",
        "mrmp2", "mcqdpt2", "xmcqdpt2",
    }:
        return

    runtype = _as_lower(_get(config, "input", "runtype", "energy"))
    scf_type = _as_lower(_get(config, "scf", "type", "rhf"))
    multiplicity = _get(config, "scf", "multiplicity", 1)
    functional = _get(config, "input", "functional", "")
    d4_enabled = _get(config, "input", "d4", False)

    active_electrons = _get(config, "cas", "active_electrons", 0)
    active_orbitals = _get(config, "cas", "active_orbitals", 0)
    frozen_core = _get(config, "cas", "frozen_core", 0)
    active_orbital_indices = _as_list(_get(config, "cas", "active_orbital_indices", []))
    core_orbital_indices = _as_list(_get(config, "cas", "core_orbital_indices", []))
    orbital_source = _check_choice_literal(
        _get(config, "cas", "orbital_source", "rhf"),
        "cas.orbital_source",
        CAS_ORBITAL_SOURCES,
        report,
        message="CASCI supports converged RHF orbitals or imported OpenQP JSON MO coefficients.",
        action="Set [cas] orbital_source=rhf or json.",
        fallback="rhf",
    )
    orbital_file_raw = _get(config, "cas", "orbital_file", "")
    localize = _check_choice_literal(
        _get(config, "cas", "localize", "none"),
        "cas.localize",
        CAS_LOCALIZE_OPTIONS,
        report,
        message="Orbital localization is planned but not implemented for CASCI.",
        action="Set [cas] localize=none.",
        fallback="none",
    )
    sort_orbitals = _check_choice_literal(
        _get(config, "cas", "sort_orbitals", "energy"),
        "cas.sort_orbitals",
        CAS_SORT_OPTIONS,
        report,
        message="Only canonical RHF energy ordering is implemented for CASCI.",
        action="Set [cas] sort_orbitals=energy or none.",
        fallback="energy",
    )
    max_det = _get(config, "cas", "max_det", 5000)
    max_memory = _get(config, "cas", "max_memory", 2048)

    nroot = _get(config, "ci", "nroot", 1)
    eig_tol = _get(config, "ci", "eig_tol", 1.0e-10)
    integral_backend = _check_choice_literal(
        _get(config, "ci", "integral_backend", "native"),
        "ci.integral_backend",
        FCI_INTEGRAL_BACKENDS,
        report,
        message="Unknown CASCI integral backend.",
        action="Use [ci] integral_backend=native.",
        fallback="native",
    )
    integral_cutoff = _get(config, "ci", "integral_cutoff", 5.0e-11)
    solver = _check_choice_literal(
        _get(config, "ci", "solver", "auto"),
        "ci.solver",
        FCI_SOLVERS,
        report,
        message="Unknown CASCI solver.",
        action="Use [ci] solver=auto, dense, or davidson.",
        fallback="auto",
    )
    davidson_maxiter = _get(config, "ci", "davidson_maxiter", 100)
    davidson_subspace = _get(config, "ci", "davidson_subspace", 0)
    _validate_bool_literal(
        _get(config, "ci", "print_ci_vectors", False),
        "ci.print_ci_vectors",
        report,
        action="Set [ci] print_ci_vectors=true or false.",
    )
    _validate_bool_literal(
        _get(config, "ci", "save_ci_vectors", False),
        "ci.save_ci_vectors",
        report,
        action="Set [ci] save_ci_vectors=true or false.",
    )
    _validate_bool_literal(
        _get(config, "ci", "save_rdm", False),
        "ci.save_rdm",
        report,
        action="Set [ci] save_rdm=true or false.",
    )
    _, spin_adapted = _validate_bool_literal(
        _get(config, "ci", "spin_adapted", False),
        "ci.spin_adapted",
        report,
        action="Set [ci] spin_adapted=false until spin-adapted CI is implemented.",
    )
    target_spin = _as_lower(_get(config, "ci", "target_spin", "any"))
    root_tracking = _check_choice_literal(
        _get(config, "ci", "root_tracking", "energy"),
        "ci.root_tracking",
        CI_ROOT_TRACKING_OPTIONS,
        report,
        message="CASCI is a single-shot fixed-orbital calculation; overlap root tracking is reserved for future CASSCF iterations.",
        action="Set [ci] root_tracking=energy for CASCI.",
        fallback="energy",
    )
    ci_print_threshold = _get(config, "ci", "ci_print_threshold", 5.0e-2)
    _, state_average_enabled = _validate_bool_literal(
        _get(config, "state_average", "enabled", False),
        "state_average.enabled",
        report,
        action="Set [state_average] enabled=true or false.",
    )
    state_average_weights_raw = _as_list(_get(config, "state_average", "weights", []))
    state_average_nstate = _get(config, "state_average", "nstate", 0)
    state_average_target_roots = _as_list(_get(config, "state_average", "target_roots", []))
    state_average_equal_weights_valid, state_average_equal_weights = _validate_bool_literal(
        _get(config, "state_average", "equal_weights", True),
        "state_average.equal_weights",
        report,
        action="Set [state_average] equal_weights=true or false.",
    )
    _check_choice_literal(
        _get(config, "state_average", "spin_blocks", "diagnostic"),
        "state_average.spin_blocks",
        STATE_AVERAGE_SPIN_BLOCKS,
        report,
        message="State-average spin-block controls are diagnostic-only in this branch.",
        action="Set [state_average] spin_blocks=diagnostic.",
        fallback="diagnostic",
    )
    _check_choice_literal(
        _get(config, "state_average", "root_tracking", "overlap"),
        "state_average.root_tracking",
        STATE_AVERAGE_ROOT_TRACKING,
        report,
        message="State-average root tracking uses conservative RDM-overlap guards across SA-CASSCF orbital iterations.",
        action="Set [state_average] root_tracking=overlap.",
        fallback="overlap",
    )

    if orbital_file_raw is None:
        orbital_file = ""
    elif isinstance(orbital_file_raw, bool) or not isinstance(orbital_file_raw, str):
        report.add(
            "ERROR",
            "cas.orbital_file",
            "CAS orbital file must be a string path.",
            value=orbital_file_raw,
            expected="string path or empty string",
            action="Set [cas] orbital_file to an OpenQP JSON restart file or leave it empty.",
        )
        orbital_file = ""
    else:
        orbital_file = orbital_file_raw.strip()

    active_electrons_valid, active_electrons_total = _check_active_electron_literal(
        active_electrons,
        "cas.active_electrons",
        report,
    )
    active_count_values: dict[str, int] = {
        "active_electrons": active_electrons_total if active_electrons_valid else 0,
    }
    for option, value in (
        ("active_orbitals", active_orbitals),
        ("frozen_core", frozen_core),
    ):
        valid_count, parsed_count = _check_int_literal(
            value,
            f"cas.{option}",
            report,
            expected="non-negative integer",
            action="Use 0 for default inference or positive active-space counts.",
        )
        active_count_values[option] = parsed_count if valid_count else 0
    active_orbitals_int = active_count_values["active_orbitals"]

    parsed_active_orbital_indices: list[int] = []
    active_orbital_indices_valid = True
    for idx in active_orbital_indices:
        valid_index, parsed_index = _check_int_literal(
            idx,
            "cas.active_orbital_indices",
            report,
            minimum=1,
            expected="positive 1-based MO index",
            action="Use 1-based MO labels, e.g. 2,3.",
        )
        if valid_index:
            parsed_active_orbital_indices.append(parsed_index)
        else:
            active_orbital_indices_valid = False

    parsed_core_orbital_indices: list[int] = []
    core_orbital_indices_valid = True
    for idx in core_orbital_indices:
        valid_index, parsed_index = _check_int_literal(
            idx,
            "cas.core_orbital_indices",
            report,
            minimum=1,
            expected="positive 1-based MO index",
            action="Use 1-based MO labels.",
        )
        if valid_index:
            parsed_core_orbital_indices.append(parsed_index)
        else:
            core_orbital_indices_valid = False

    _check_int_literal(
        max_det,
        "cas.max_det",
        report,
        minimum=1,
        expected="positive integer",
        action="Increase [cas] max_det.",
    )
    _check_int_literal(
        max_memory,
        "cas.max_memory",
        report,
        minimum=1,
        expected="positive integer MiB",
        action="Set [cas] max_memory to a positive size in MiB.",
    )
    nroot_valid, nroot_int = _check_int_literal(
        nroot,
        "ci.nroot",
        report,
        minimum=1,
        expected="positive integer",
        action="Set [ci] nroot=1 or higher.",
    )
    if not nroot_valid:
        nroot_int = 1
    eig_tol_valid, eig_tol_float = _check_float_literal(
        eig_tol,
        "ci.eig_tol",
        report,
        expected="positive finite number",
        action="Set [ci] eig_tol to a positive threshold.",
    )
    integral_cutoff_valid, integral_cutoff_float = _check_float_literal(
        integral_cutoff,
        "ci.integral_cutoff",
        report,
        expected="finite non-negative number",
        action="Set [ci] integral_cutoff to zero or a positive threshold.",
    )
    _check_int_literal(
        davidson_maxiter,
        "ci.davidson_maxiter",
        report,
        minimum=1,
        expected="positive integer",
        action="Set [ci] davidson_maxiter to a positive integer.",
    )
    davidson_subspace_valid, davidson_subspace_int = _check_int_literal(
        davidson_subspace,
        "ci.davidson_subspace",
        report,
        minimum=0,
        expected="0 or positive integer",
        action="Use 0 for the internal default.",
    )
    ci_print_threshold_valid, ci_print_threshold_float = _check_float_literal(
        ci_print_threshold,
        "ci.ci_print_threshold",
        report,
        expected="finite non-negative number",
        action="Set [ci] ci_print_threshold to zero or a positive value.",
    )

    # The PT2 family carries central-difference numerical gradients
    # (pt2_numgrad.py) wired through the Gradient seam, so every
    # gradient-consuming driver works for it; casci/casscf remain energy-only.
    pt2_grad_runtypes = {"energy", "grad", "optimize", "meci", "mecp", "ts",
                         "mep", "neb", "irc"}
    if method in PT2_METHOD_ALIASES:
        if runtype not in pt2_grad_runtypes:
            report.add(
                "ERROR",
                "input.runtype",
                "This PT2 runtype is not wired (energies + numerical-gradient "
                "drivers only).",
                value=runtype,
                expected=", ".join(sorted(pt2_grad_runtypes)),
                action="Use energy, grad, or a gradient-driven optimizer runtype "
                       "(optimize/meci/mecp/ts/mep/neb/irc).",
            )
        grad_step = _get(config, "pt2", "grad_step", 1.0e-3)
        valid_gs, gs = _check_float_literal(grad_step, "pt2.grad_step", report)
        if valid_gs and gs <= 0.0:
            report.add(
                "ERROR",
                "pt2.grad_step",
                "The numerical-gradient displacement must be positive.",
                value=grad_step,
                expected="> 0 (Bohr)",
                action="Set [pt2] grad_step to a positive displacement (default 1e-3).",
            )
        grad_guess = str(_get(config, "pt2", "grad_guess", "cold") or "cold").strip().lower()
        if grad_guess not in {"cold", "warm"}:
            report.add(
                "ERROR",
                "pt2.grad_guess",
                "Unknown displaced-point guess policy.",
                value=grad_guess,
                expected="cold or warm",
                action="Use grad_guess=cold (fresh guess, bit-reproducible) or warm "
                       "(reuse in-memory MOs).",
            )
        _check_float_literal(_get(config, "pt2", "grad_gap_warn", 1.0e-5),
                             "pt2.grad_gap_warn", report)
    elif runtype != "energy":
        report.add(
            "ERROR",
            "input.runtype",
            "CASCI is currently implemented only for energy calculations.",
            value=runtype,
            expected="energy",
            action="Set [input] runtype=energy.",
        )

    if scf_type != "rhf" or multiplicity != 1:
        report.add(
            "ERROR",
            "scf.type",
            "CASCI MVP requires a closed-shell RHF reference.",
            value=f"{scf_type}/{multiplicity}",
            expected="scf.type=rhf and scf.multiplicity=1",
            action="Use an RHF singlet reference for method=casci.",
        )

    if functional:
        report.add(
            "ERROR",
            "input.functional",
            "CASCI requires HF one- and two-electron integrals.",
            value=functional,
            action="Unset [input] functional for method=casci.",
        )

    if d4_enabled:
        report.add(
            "ERROR",
            "input.d4",
            "CASCI reports the variational electronic energy; D4 dispersion correction is not part of the CI solver.",
            value=d4_enabled,
            expected=False,
            action="Disable [input] d4 for method=casci.",
        )

    if min(active_count_values.values()) < 0:
        report.add(
            "ERROR",
            "cas.active_space",
            "CAS active-space counts cannot be negative.",
            value={
                "active_electrons": active_electrons,
                "active_orbitals": active_orbitals,
                "frozen_core": frozen_core,
            },
            expected="non-negative integers",
            action="Use 0 for default inference or positive active-space counts.",
        )

    if active_orbital_indices_valid and parsed_active_orbital_indices:
        if len(set(parsed_active_orbital_indices)) != len(parsed_active_orbital_indices):
            report.add("ERROR", "cas.active_orbital_indices", "Active orbital indices must be unique.", value=active_orbital_indices, action="Remove duplicate MO indices.")
        if active_orbitals_int and active_orbitals_int != len(parsed_active_orbital_indices):
            report.add("ERROR", "cas.active_orbitals", "active_orbitals must match the length of active_orbital_indices.", value={"active_orbitals": active_orbitals, "indices": active_orbital_indices}, action="Set active_orbitals to the list length or omit it.")

    if core_orbital_indices_valid and parsed_core_orbital_indices:
        if len(set(parsed_core_orbital_indices)) != len(parsed_core_orbital_indices):
            report.add("ERROR", "cas.core_orbital_indices", "Core orbital indices must be unique.", value=core_orbital_indices, action="Remove duplicate MO indices.")
        if parsed_core_orbital_indices and not parsed_active_orbital_indices:
            report.add("ERROR", "cas.active_orbital_indices", "Explicit core_orbital_indices require active_orbital_indices.", value={"active": active_orbital_indices, "core": core_orbital_indices}, action="Provide explicit active orbital indices or use frozen_core for sequential selection.")
        if parsed_active_orbital_indices and set(parsed_active_orbital_indices) & set(parsed_core_orbital_indices):
            report.add("ERROR", "cas.core_orbital_indices", "Core and active orbital index lists must not overlap.", value={"active": active_orbital_indices, "core": core_orbital_indices}, action="Choose disjoint core and active orbital lists.")

    if orbital_source not in CAS_ORBITAL_SOURCES:
        report.add(
            "ERROR",
            "cas.orbital_source",
            "CASCI supports converged RHF orbitals or imported OpenQP JSON MO coefficients.",
            value=orbital_source,
            expected=", ".join(sorted(CAS_ORBITAL_SOURCES)),
            action="Set [cas] orbital_source=rhf or json.",
        )
    elif orbital_source == "json" and not orbital_file:
        report.add("ERROR", "cas.orbital_file", f"orbital_source={orbital_source} requires an orbital_file path.", value=orbital_file, action="Set [cas] orbital_file to an OpenQP JSON restart file.")

    if localize not in CAS_LOCALIZE_OPTIONS:
        report.add(
            "ERROR",
            "cas.localize",
            "Orbital localization is planned but not implemented for CASCI.",
            value=localize,
            expected=", ".join(sorted(CAS_LOCALIZE_OPTIONS)),
            action="Set [cas] localize=none.",
        )

    if sort_orbitals not in CAS_SORT_OPTIONS:
        report.add(
            "ERROR",
            "cas.sort_orbitals",
            "Only canonical RHF energy ordering is implemented for CASCI.",
            value=sort_orbitals,
            expected=", ".join(sorted(CAS_SORT_OPTIONS)),
            action="Set [cas] sort_orbitals=energy or none.",
        )

    if eig_tol_valid and eig_tol_float <= 0:
        report.add(
            "ERROR",
            "ci.eig_tol",
            "CASCI eigensolver tolerance must be positive.",
            value=eig_tol,
            expected="> 0",
            action="Set [ci] eig_tol to a positive threshold.",
        )

    if integral_cutoff_valid and integral_cutoff_float < 0:
        report.add(
            "ERROR",
            "ci.integral_cutoff",
            "CASCI integral cutoff cannot be negative.",
            value=integral_cutoff,
            expected=">= 0",
            action="Set [ci] integral_cutoff to zero or a positive threshold.",
        )

    if integral_backend not in FCI_INTEGRAL_BACKENDS:
        report.add(
            "ERROR",
            "ci.integral_backend",
            "Unknown CASCI integral backend.",
            value=integral_backend,
            expected=", ".join(sorted(FCI_INTEGRAL_BACKENDS)),
            action="Use [ci] integral_backend=native.",
        )

    if solver not in FCI_SOLVERS:
        report.add(
            "ERROR",
            "ci.solver",
            "Unknown CASCI solver.",
            value=solver,
            expected=", ".join(sorted(FCI_SOLVERS)),
            action="Use [ci] solver=auto, dense, or davidson.",
        )

    if (
        nroot_valid
        and davidson_subspace_valid
        and davidson_subspace_int > 0
        and davidson_subspace_int < nroot_int
    ):
        report.add(
            "ERROR",
            "ci.davidson_subspace",
            "CASCI Davidson subspace override must be at least the requested root count.",
            value=davidson_subspace,
            expected=f"0 or >= nroot ({nroot_int})",
            action="Use 0 for the internal default or raise [ci] davidson_subspace.",
        )

    if spin_adapted:
        report.add(
            "ERROR",
            "ci.spin_adapted",
            "Spin-adapted CI is not implemented; the current solver uses alpha/beta determinants with spin diagnostics.",
            value=spin_adapted,
            expected=False,
            action="Set [ci] spin_adapted=false.",
        )

    if not _target_spin_is_valid(target_spin):
        report.add(
            "ERROR",
            "ci.target_spin",
            "Unknown target spin label.",
            value=target_spin,
            expected="any, a named spin multiplicity, or a positive integer multiplicity",
            action="Use [ci] target_spin=any, singlet, doublet, triplet, quartet, or a positive integer multiplicity.",
        )

    if ci_print_threshold_valid and ci_print_threshold_float < 0:
        report.add(
            "ERROR",
            "ci.ci_print_threshold",
            "CI print threshold cannot be negative.",
            value=ci_print_threshold,
            expected=">= 0",
            action="Set [ci] ci_print_threshold to zero or a positive value.",
        )

    if state_average_enabled:
        state_average_nstate_valid, state_average_nstate_int = _check_int_literal(
            state_average_nstate,
            "state_average.nstate",
            report,
            minimum=0,
            expected="0 or a positive integer",
            action="Set [state_average] nstate to 0 for ci.nroot or to the number of averaged roots.",
        )
        if not state_average_nstate_valid:
            state_average_nstate_int = 0

        parsed_target_roots = []
        target_roots_valid = True
        for root in state_average_target_roots:
            valid_root, parsed_root = _check_int_literal(
                root,
                "state_average.target_roots",
                report,
                minimum=0,
                expected="0-based CI root slots",
                action="Use values such as 0,1 for the first two solved CI roots.",
            )
            if not valid_root:
                target_roots_valid = False
                break
            parsed_target_roots.append(parsed_root)

        if target_roots_valid and len(set(parsed_target_roots)) != len(parsed_target_roots):
            report.add(
                "ERROR",
                "state_average.target_roots",
                "State-average target roots must be unique.",
                value=state_average_target_roots,
                action="Remove duplicate root slots.",
            )
        if target_roots_valid and parsed_target_roots and max(parsed_target_roots) >= nroot_int:
            report.add(
                "ERROR",
                "state_average.target_roots",
                "State-average target roots must refer to solved CI root slots.",
                value=state_average_target_roots,
                expected=f"all roots < ci.nroot ({nroot_int})",
                action="Increase [ci] nroot or lower [state_average] target_roots.",
            )
        if (
            target_roots_valid
            and parsed_target_roots
            and state_average_nstate_int not in (0, len(parsed_target_roots))
        ):
            report.add(
                "ERROR",
                "state_average.nstate",
                "State-average nstate must match target_roots when explicit target roots are provided.",
                value={
                    "nstate": state_average_nstate,
                    "target_roots": state_average_target_roots,
                },
                expected=f"0 or {len(parsed_target_roots)}",
                action="Set nstate to the number of target roots or leave it as 0.",
            )

        averaged_root_count = (
            len(parsed_target_roots)
            if target_roots_valid and parsed_target_roots
            else int(state_average_nstate_int or nroot_int)
        )
        if averaged_root_count < 1:
            report.add(
                "ERROR",
                "state_average.nstate",
                "State averaging requires at least one root.",
                value=averaged_root_count,
                expected=">= 1",
                action="Request at least one CI root and average at least one state.",
            )
        if nroot_valid and averaged_root_count > nroot_int:
            report.add(
                "ERROR",
                "state_average.nstate",
                "State-average root count cannot exceed ci.nroot.",
                value=averaged_root_count,
                expected=f"<= ci.nroot ({nroot_int})",
                action="Increase [ci] nroot or lower [state_average] nstate.",
            )

        if state_average_equal_weights_valid and not state_average_equal_weights:
            parsed_weights = []
            weights_valid = True
            for weight in state_average_weights_raw:
                valid_weight, parsed_weight = _check_float_literal(
                    weight,
                    "state_average.weights",
                    report,
                    expected="finite non-negative number",
                    action="Use non-negative numeric weights.",
                )
                if not valid_weight:
                    weights_valid = False
                    break
                parsed_weights.append(parsed_weight)

            if not weights_valid:
                pass
            elif len(parsed_weights) != averaged_root_count:
                report.add(
                    "ERROR",
                    "state_average.weights",
                    "State-average weights must match the averaged root count.",
                    value=state_average_weights_raw,
                    expected=f"{averaged_root_count} weights",
                    action="Provide one weight for each averaged root or set equal_weights=true.",
                )
            elif not all(math.isfinite(weight) for weight in parsed_weights):
                report.add(
                    "ERROR",
                    "state_average.weights",
                    "State-average weights must be finite.",
                    value=state_average_weights_raw,
                    action="Use finite non-negative weights.",
                )
            elif any(weight < 0.0 for weight in parsed_weights):
                report.add(
                    "ERROR",
                    "state_average.weights",
                    "State-average weights cannot be negative.",
                    value=state_average_weights_raw,
                    action="Use non-negative weights.",
                )
            elif sum(parsed_weights) <= 0.0:
                report.add(
                    "ERROR",
                    "state_average.weights",
                    "State-average weights must have positive total weight.",
                    value=state_average_weights_raw,
                    action="Set at least one positive weight.",
                )


def _check_casscf(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    if method not in {"casscf", "sa-casscf", "sacasscf"}:
        return
    is_sa_alias = method in {"sa-casscf", "sacasscf"}

    casscf_section = config.get("casscf", {})
    raw_max_macro_iterations = _get(config, "casscf", "max_macro_iterations", 20)
    if isinstance(raw_max_macro_iterations, bool):
        report.add(
            "ERROR",
            "casscf.max_macro_iterations",
            "CASSCF macro-iteration count must be an integer, not a boolean.",
            value=raw_max_macro_iterations,
            expected="non-negative integer",
            action="Set [casscf] max_macro_iterations to 0 for the fixed-orbital scaffold or a positive macroiteration budget.",
            wiki=WIKI_HELP["input.method"],
        )
        return

    try:
        max_macro_iterations = int(raw_max_macro_iterations)
    except (TypeError, ValueError):
        report.add(
            "ERROR",
            "casscf.max_macro_iterations",
            "CASSCF macro-iteration count must be an integer.",
            value=raw_max_macro_iterations,
            expected="non-negative integer",
            action="Set [casscf] max_macro_iterations to 0 for the fixed-orbital scaffold or a positive macroiteration budget.",
            wiki=WIKI_HELP["input.method"],
        )
        return

    if max_macro_iterations < 0:
        report.add(
            "ERROR",
            "casscf.max_macro_iterations",
            "CASSCF macro-iteration count must be non-negative.",
            value=max_macro_iterations,
            expected=">= 0",
            action="Use 0 for the fixed-orbital scaffold or a positive value for native CASSCF orbital optimization.",
            wiki=WIKI_HELP["input.method"],
        )
    elif max_macro_iterations == 0:
        report.add(
            "WARNING",
            "input.method",
            "method=casscf with max_macro_iterations=0 runs the fixed-orbital native CASCI path; no orbital optimization is performed.",
            value=method,
            action="Use method=casci for ordinary fixed-orbital active-space CI, or keep this explicit scaffold setting for compatibility tests.",
            wiki=WIKI_HELP["input.method"],
        )

    converger = str(_get(config, "casscf", "converger", "twophase") or "twophase").strip().lower()
    if converger not in {"twophase", "two-phase", "default", "", "ah", "trah",
                         "augmented-hessian", "diis", "auto"}:
        report.add(
            "ERROR",
            "casscf.converger",
            "Unknown CASSCF orbital converger.",
            value=converger,
            expected="twophase, ah, trah, diis, or auto",
            action="Use converger=twophase (default), ah for the trust-region "
                   "augmented Hessian over an assembled orbital Hessian, trah "
                   "for the matrix-free trust-region solver (no Hessian is "
                   "built), diis for orbital-gradient DIIS, or auto for ah "
                   "with a two-phase fallback.",
        )
    hessian_mode = str(_get(config, "casscf", "hessian", "fd") or "fd").strip().lower()
    if hessian_mode not in {"fd", "finite-difference", "finite_difference",
                            "analytic", "analytical", "exact"}:
        report.add(
            "ERROR",
            "casscf.hessian",
            "Unknown CASSCF orbital-Hessian builder.",
            value=hessian_mode,
            expected="fd or analytic",
            action="Use hessian=fd (default finite-difference) or "
                   "hessian=analytic (exact orbital Hessian; ~1 CI solve per "
                   "macroiteration instead of 2*n_par+1).",
        )
    for key, default, minimum in (
        ("ah_start_trust_radius", 0.2, 0.0),
        ("ah_max_trust_radius", 0.0, None),   # <= 0 means auto (max_rotation_norm)
        ("ah_min_trust_radius", 1.0e-6, 0.0),
        ("ah_saddle_curv_tol", 2.5e-2, 0.0),
        ("ah_saddle_egain_tol", 1.0e-3, 0.0),
    ):
        valid_f, numeric = _check_float_literal(
            _get(config, "casscf", key, default), f"casscf.{key}", report)
        if valid_f and minimum is not None and numeric <= minimum:
            report.add(
                "ERROR",
                f"casscf.{key}",
                "CASSCF AH converger tolerance must be positive.",
                value=numeric,
                expected="> 0",
                action=f"Set [casscf] {key} to a positive value.",
            )
    for key, default in (("ah_max_micro", 32), ("ah_max_rejects", 6),
                         ("diis_space", 8), ("diis_start", 2),
                         ("auto_stagnation", 3)):
        _check_int_literal(
            _get(config, "casscf", key, default), f"casscf.{key}", report,
            minimum=1, expected="positive integer",
            action=f"Set [casscf] {key} to a positive integer.",
        )

    root_valid, casscf_root = _check_int_literal(
        _get(config, "casscf", "root", 0),
        "casscf.root",
        report,
        minimum=0,
        expected="0-based non-negative integer root index",
        action="Use 0 for the ground-state CASSCF root.",
    )
    nroot_valid, nroot_int = _parse_int_literal(_get(config, "ci", "nroot", 1))
    if root_valid and nroot_valid and casscf_root >= nroot_int:
        report.add(
            "ERROR",
            "casscf.root",
            "CASSCF root must be smaller than ci.nroot.",
            value=casscf_root,
            expected=f"0 <= root < {nroot_int}",
            action="Increase [ci] nroot or choose a lower CASSCF root.",
            wiki=WIKI_HELP["input.method"],
        )

    optimizer = _as_lower(_get(config, "casscf", "optimizer", "newton")).replace("_", "-")
    if optimizer not in {"newton", "newton-ah", "augmented-hessian", "microiteration",
                         "steepest-descent", "steepest", "diagnostic", "powell", "scipy-powell", "scipy"}:
        report.add(
            "ERROR",
            "casscf.optimizer",
            "CASSCF orbital optimizer is not recognized.",
            value=optimizer,
            expected="newton or powell",
            action="Use optimizer=newton for the native Newton orbital optimizer (default) or optimizer=powell for the gradient-only fallback.",
            wiki=WIKI_HELP["input.method"],
        )

    _check_int_literal(
        _get(config, "casscf", "max_function_evaluations", 0),
        "casscf.max_function_evaluations",
        report,
        minimum=0,
        expected="non-negative integer",
        action="Use 0 for SciPy's default function-evaluation budget, or a positive cap.",
    )

    macro_gradient_valid, macro_gradient_tol = _check_float_literal(
        _get(config, "casscf", "gradient_norm_tol", 1.0e-6),
        "casscf.gradient_norm_tol",
        report,
        expected="positive finite number",
        action="Set a positive finite CASSCF orbital-gradient tolerance.",
    )
    if macro_gradient_valid and macro_gradient_tol <= 0.0:
        report.add(
            "ERROR",
            "casscf.gradient_norm_tol",
            "CASSCF orbital-gradient tolerance must be positive.",
            value=macro_gradient_tol,
            expected="positive finite number",
            action="Set gradient_norm_tol to a positive finite value.",
            wiki=WIKI_HELP["input.method"],
        )

    macro_energy_valid, macro_energy_tol = _check_float_literal(
        _get(config, "casscf", "energy_decrease_tol", 1.0e-10),
        "casscf.energy_decrease_tol",
        report,
        expected="finite non-negative number",
        action="Set a non-negative finite CASSCF energy-decrease tolerance.",
    )
    if macro_energy_valid and macro_energy_tol < 0.0:
        report.add(
            "ERROR",
            "casscf.energy_decrease_tol",
            "CASSCF energy-decrease tolerance must be non-negative.",
            value=macro_energy_tol,
            expected="finite non-negative number",
            action="Set energy_decrease_tol to zero or a positive finite value.",
            wiki=WIKI_HELP["input.method"],
        )

    macro_step_valid, macro_step_tol = _check_float_literal(
        _get(config, "casscf", "step_norm_tol", 1.0e-8),
        "casscf.step_norm_tol",
        report,
        expected="finite non-negative number",
        action="Set a non-negative finite CASSCF step-norm tolerance.",
    )
    if macro_step_valid and macro_step_tol < 0.0:
        report.add(
            "ERROR",
            "casscf.step_norm_tol",
            "CASSCF step-norm tolerance must be non-negative.",
            value=macro_step_tol,
            expected="finite non-negative number",
            action="Set step_norm_tol to zero or a positive finite value.",
            wiki=WIKI_HELP["input.method"],
        )

    macro_rotation_valid, macro_rotation_norm = _check_float_literal(
        _get(config, "casscf", "max_rotation_norm", 5.0e-2),
        "casscf.max_rotation_norm",
        report,
        expected="positive finite number",
        action="Set a positive finite CASSCF maximum orbital-rotation norm.",
    )
    if macro_rotation_valid and macro_rotation_norm <= 0.0:
        report.add(
            "ERROR",
            "casscf.max_rotation_norm",
            "CASSCF maximum orbital-rotation norm must be positive.",
            value=macro_rotation_norm,
            expected="positive finite number",
            action="Set max_rotation_norm to a positive finite value.",
            wiki=WIKI_HELP["input.method"],
        )

    diagnostic_report = _check_bool_literal(
        _get(config, "casscf", "diagnostic_report", False),
        "casscf.diagnostic_report",
        report,
        action="Use true to write a diagnostic-only CASSCF JSON report, or false to disable it.",
    )
    diagnostic_report_file = _get(config, "casscf", "diagnostic_report_file", "")
    if isinstance(diagnostic_report_file, bool) or not isinstance(diagnostic_report_file, str):
        report.add(
            "ERROR",
            "casscf.diagnostic_report_file",
            "CASSCF diagnostic report file must be a filesystem path string.",
            value=diagnostic_report_file,
            expected="string path or blank for the default run-log sidecar",
            action="Set [casscf] diagnostic_report_file to a .json path, or leave it blank.",
            wiki=WIKI_HELP["input.method"],
        )
    elif (
        diagnostic_report
        and diagnostic_report_file.strip()
        and not diagnostic_report_file.strip().lower().endswith(".json")
    ):
        report.add(
            "WARNING",
            "casscf.diagnostic_report_file",
            "CASSCF diagnostic report files should use a .json suffix.",
            value=diagnostic_report_file,
            action="Use a .json file name for portable diagnostic artifacts.",
            wiki=WIKI_HELP["input.method"],
        )

    benchmark_reference_file = _get(config, "casscf", "diagnostic_benchmark_reference_file", "")
    if isinstance(benchmark_reference_file, bool) or not isinstance(benchmark_reference_file, str):
        report.add(
            "ERROR",
            "casscf.diagnostic_benchmark_reference_file",
            "CASSCF diagnostic benchmark reference file must be a filesystem path string.",
            value=benchmark_reference_file,
            expected="string path or blank when no diagnostic benchmark is configured",
            action="Set [casscf] diagnostic_benchmark_reference_file to a .json path, or leave it blank.",
            wiki=WIKI_HELP["input.method"],
        )
        benchmark_reference_path = ""
    else:
        benchmark_reference_path = benchmark_reference_file.strip()
        if benchmark_reference_path and not benchmark_reference_path.lower().endswith(".json"):
            report.add(
                "WARNING",
                "casscf.diagnostic_benchmark_reference_file",
                "CASSCF diagnostic benchmark reference files should use a .json suffix.",
                value=benchmark_reference_file,
                action="Use a .json file name for portable benchmark artifacts.",
                wiki=WIKI_HELP["input.method"],
            )

    benchmark_required = _check_bool_literal(
        _get(config, "casscf", "diagnostic_benchmark_required", False),
        "casscf.diagnostic_benchmark_required",
        report,
        action="Use true only when a diagnostic benchmark reference file is configured.",
    )
    benchmark_tol_valid, benchmark_tol = _check_float_literal(
        _get(config, "casscf", "diagnostic_benchmark_tolerance", 1.0e-5),
        "casscf.diagnostic_benchmark_tolerance",
        report,
        expected="finite non-negative number",
        action="Set a non-negative finite diagnostic benchmark tolerance.",
    )
    if benchmark_tol_valid and benchmark_tol < 0.0:
        report.add(
            "ERROR",
            "casscf.diagnostic_benchmark_tolerance",
            "CASSCF diagnostic benchmark tolerance must be non-negative.",
            value=benchmark_tol,
            expected="finite non-negative number",
            action="Set diagnostic_benchmark_tolerance to zero or a positive finite value.",
            wiki=WIKI_HELP["input.method"],
        )
    if benchmark_required and not benchmark_reference_path:
        report.add(
            "ERROR",
            "casscf.diagnostic_benchmark_reference_file",
            "Required CASSCF diagnostic benchmarks need a reference file.",
            value=benchmark_reference_file,
            expected="path to a CASSCF diagnostic benchmark reference JSON",
            action="Set diagnostic_benchmark_reference_file or disable diagnostic_benchmark_required.",
            wiki=WIKI_HELP["input.method"],
        )
    if (benchmark_required or benchmark_reference_path) and not diagnostic_report:
        report.add(
            "ERROR",
            "casscf.diagnostic_report",
            "CASSCF diagnostic benchmark comparison requires diagnostic_report=true.",
            value=diagnostic_report,
            expected=True,
            action="Enable [casscf] diagnostic_report=true or remove the diagnostic benchmark options.",
            wiki=WIKI_HELP["input.method"],
        )

    root_valid, diagnostic_root = _check_int_literal(
        _get(config, "casscf", "diagnostic_root", 0),
        "casscf.diagnostic_root",
        report,
        minimum=0,
        expected="0-based non-negative integer root index",
        action="Use 0 for the ground-state diagnostic root.",
    )
    nroot_valid, nroot_int = _parse_int_literal(_get(config, "ci", "nroot", 1))
    if root_valid and nroot_valid and diagnostic_root >= nroot_int:
        report.add(
            "ERROR",
            "casscf.diagnostic_root",
            "CASSCF diagnostic root must be smaller than ci.nroot.",
            value=diagnostic_root,
            expected=f"0 <= diagnostic_root < {nroot_int}",
            action="Increase [ci] nroot or choose a lower diagnostic_root.",
            wiki=WIKI_HELP["input.method"],
        )

    _check_int_literal(
        _get(config, "casscf", "diagnostic_max_iterations", 1),
        "casscf.diagnostic_max_iterations",
        report,
        minimum=1,
        expected="positive integer",
        action="Use a small positive diagnostic iteration count.",
    )
    gradient_valid, gradient_tol = _check_float_literal(
        _get(config, "casscf", "diagnostic_gradient_norm_tol", 1.0e-8),
        "casscf.diagnostic_gradient_norm_tol",
        report,
        expected="positive finite number",
        action="Set a positive finite diagnostic gradient-norm tolerance.",
    )
    if gradient_valid and gradient_tol <= 0.0:
        report.add(
            "ERROR",
            "casscf.diagnostic_gradient_norm_tol",
            "CASSCF diagnostic gradient-norm tolerance must be positive.",
            value=gradient_tol,
            expected="positive finite number",
            action="Set diagnostic_gradient_norm_tol to a positive finite value.",
            wiki=WIKI_HELP["input.method"],
        )
    rotation_valid, rotation_norm = _check_float_literal(
        _get(config, "casscf", "diagnostic_max_rotation_norm", 5.0e-2),
        "casscf.diagnostic_max_rotation_norm",
        report,
        expected="positive finite number",
        action="Set a positive finite diagnostic maximum rotation norm.",
    )
    if rotation_valid and rotation_norm <= 0.0:
        report.add(
            "ERROR",
            "casscf.diagnostic_max_rotation_norm",
            "CASSCF diagnostic maximum rotation norm must be positive.",
            value=rotation_norm,
            expected="positive finite number",
            action="Set diagnostic_max_rotation_norm to a positive finite value.",
            wiki=WIKI_HELP["input.method"],
        )
    if diagnostic_report and not casscf_section:
        report.add(
            "ERROR",
            "casscf",
            "CASSCF diagnostic reporting requires an explicit [casscf] section.",
            value=casscf_section,
            expected="[casscf] with diagnostic_report=true",
            action="Add an explicit [casscf] section when requesting a diagnostic report.",
            wiki=WIKI_HELP["input.method"],
        )

    _, state_average_enabled = _parse_bool_literal(
        _get(config, "state_average", "enabled", False)
    )
    if is_sa_alias and not state_average_enabled:
        report.add(
            "ERROR",
            "state_average.enabled",
            "method=sa-casscf requires state averaging to be enabled.",
            value=state_average_enabled,
            expected=True,
            action="Add [state_average] enabled=true with nstate/target_roots/weights for SA-CASSCF.",
            wiki=WIKI_HELP["input.method"],
        )
    elif state_average_enabled and max_macro_iterations > 0:
        report.add(
            "INFO",
            "state_average.enabled",
            "Native SA-CASSCF orbital optimization is enabled; macroiterations use the weighted state-average CI/RDM objective.",
            value=state_average_enabled,
            action="Keep [state_average] target_roots/weights consistent with the solved CI roots.",
            wiki=WIKI_HELP["input.method"],
        )
    elif state_average_enabled:
        report.add(
            "WARNING",
            "state_average.enabled",
            "State-average bookkeeping is enabled with max_macro_iterations=0; no orbital optimization is performed.",
            value=state_average_enabled,
            action="Use a positive max_macro_iterations value for SA-CASSCF orbital optimization, or method=casci for fixed-orbital state averaging.",
            wiki=WIKI_HELP["input.method"],
        )


def _check_pt2(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    pt2 = config.get("pt2", {})
    reference_report = _check_bool_literal(
        _get(config, "pt2", "reference_report", False),
        "pt2.reference_report",
        report,
        action="Use true to write a diagnostic-only CASPT2 reference preflight report, or false to disable it.",
    )
    canonical_method = PT2_METHOD_ALIASES.get(method)
    if canonical_method is None:
        if not reference_report:
            return
        if method not in {"casci", "casscf"}:
            report.add(
                "ERROR",
                "pt2.reference_report",
                "CASPT2 reference preflight reports require method=casci or method=casscf.",
                value=method,
                expected="casci or casscf",
                action="Use method=casci for a CASCI reference preflight, or method=casscf for a CASSCF reference preflight.",
                wiki=WIKI_HELP["input.method"],
            )
            return
        raw_variant = _get(config, "pt2", "variant", "auto")
        if isinstance(raw_variant, str) and raw_variant.strip().lower() in {"", "auto"}:
            canonical_method = "caspt2"
        elif isinstance(raw_variant, str) and raw_variant.strip().lower() in PT2_METHOD_ALIASES:
            canonical_method = PT2_METHOD_ALIASES[raw_variant.strip().lower()]
        else:
            canonical_method = "caspt2"

    ci_nroot = _get(config, "ci", "nroot", 1)
    variant = _check_choice_literal(
        _get(config, "pt2", "variant", "auto"),
        "pt2.variant",
        PT2_VARIANTS,
        report,
        message="Unknown CASPT2-family variant.",
        action="Use auto, caspt2, ms-caspt2, xms-caspt2, mrmp2, mcqdpt2, or xmcqdpt2.",
        fallback="auto",
        default_if_none=True,
        default_if_blank=True,
    )
    reference = _check_choice_literal(
        _get(config, "pt2", "reference", "casscf"),
        "pt2.reference",
        PT2_REFERENCES,
        report,
        message="Unknown PT2 reference type.",
        action="Use reference=casci for fixed-orbital CI or reference=casscf after native CASSCF exists.",
        fallback="casscf",
    )
    contraction = _check_choice_literal(
        _get(config, "pt2", "contraction", "uncontracted"),
        "pt2.contraction",
        PT2_CONTRACTIONS,
        report,
        message="Unknown PT2 contraction mode.",
        action="Use contraction=uncontracted (determinant-space CASPT2/NEVPT2) or "
               "contraction=strong (RDM-based SC-NEVPT2; requires h0=dyall).",
        fallback="uncontracted",
    )
    multistate = _check_choice_literal(
        _get(config, "pt2", "multistate", "auto"),
        "pt2.multistate",
        PT2_MULTISTATE_MODES,
        report,
        message="Unknown PT2 multistate mode.",
        action="Use auto, none, ms, or xms.",
        fallback="auto",
        default_if_none=True,
        default_if_blank=True,
    )
    xms_requested = _check_bool_literal(
        _get(config, "pt2", "xms", False),
        "pt2.xms",
        report,
        action="Set [pt2] xms=true only for XMS-CASPT2, otherwise false.",
    )
    ipea_shift = _get(config, "pt2", "ipea_shift", 0.0)
    imaginary_shift = _get(config, "pt2", "imaginary_shift", 0.0)
    level_shift = _get(config, "pt2", "level_shift", 0.0)
    edshft = _get(config, "pt2", "edshft", 0.0)
    denominator_cutoff = _get(config, "pt2", "denominator_cutoff", 1.0e-10)
    intruder_threshold = _get(config, "pt2", "intruder_threshold", 1.0e-6)
    nroot = _get(config, "pt2", "nroot", 0)
    target_roots = _as_list(_get(config, "pt2", "target_roots", []))
    max_memory = _get(config, "pt2", "max_memory", 2048)
    semi_canonical = _check_bool_literal(
        _get(config, "pt2", "semi_canonical", True),
        "pt2.semi_canonical",
        report,
        action="Keep [pt2] semi_canonical=true until noncanonical PT2 is implemented.",
    )
    save_amplitudes = _check_bool_literal(
        _get(config, "pt2", "save_amplitudes", False),
        "pt2.save_amplitudes",
        report,
        action="Set [pt2] save_amplitudes=true or false.",
    )
    print_amplitudes = _check_bool_literal(
        _get(config, "pt2", "print_amplitudes", False),
        "pt2.print_amplitudes",
        report,
        action="Set [pt2] print_amplitudes=true or false.",
    )
    amplitude_threshold = _get(config, "pt2", "amplitude_threshold", 1.0e-2)
    benchmark_reference_file_raw = _get(config, "pt2", "benchmark_reference_file", "")
    benchmark_tolerance = _get(config, "pt2", "benchmark_tolerance", 1.0e-4)
    benchmark_required = _check_bool_literal(
        _get(config, "pt2", "benchmark_required", False),
        "pt2.benchmark_required",
        report,
        action="Set [pt2] benchmark_required=true only when benchmark_reference_file is configured.",
    )
    reference_report_file_raw = _get(config, "pt2", "reference_report_file", "")
    if isinstance(reference_report_file_raw, bool) or not isinstance(reference_report_file_raw, str):
        report.add(
            "ERROR",
            "pt2.reference_report_file",
            "CASPT2 reference preflight report file must be a filesystem path string.",
            value=reference_report_file_raw,
            expected="string path or blank for the default run-log sidecar",
            action="Set [pt2] reference_report_file to a .json path, or leave it blank.",
            wiki=WIKI_HELP["input.method"],
        )
    elif (
        reference_report
        and reference_report_file_raw.strip()
        and not reference_report_file_raw.strip().lower().endswith(".json")
    ):
        report.add(
            "WARNING",
            "pt2.reference_report_file",
            "CASPT2 reference preflight report files should use a .json suffix.",
            value=reference_report_file_raw,
            action="Use a .json file name for portable diagnostic artifacts.",
            wiki=WIKI_HELP["input.method"],
        )

    if variant not in PT2_VARIANTS:
        report.add(
            "ERROR",
            "pt2.variant",
            "Unknown CASPT2-family variant.",
            value=variant,
            expected=", ".join(sorted(PT2_VARIANTS)),
            action="Use auto, caspt2, ms-caspt2, xms-caspt2, mrmp2, mcqdpt2, or xmcqdpt2.",
        )
        normalized_variant = canonical_method
    elif variant == "auto":
        normalized_variant = canonical_method
    else:
        normalized_variant = variant

    if variant not in {"auto", canonical_method}:
        report.add(
            "ERROR",
            "pt2.variant",
            "PT2 variant must agree with the requested method label.",
            value={"method": canonical_method, "variant": variant},
            expected=f"auto or {canonical_method}",
            action="Make [pt2] variant match [input] method.",
        )

    if canonical_method in PT2_QDPT_METHODS:
        h0_mode = str(_get(config, "pt2", "h0", "fock") or "fock").strip().lower()
        if h0_mode not in {"fock", "caspt2", ""}:
            report.add(
                "ERROR",
                "pt2.h0",
                "MRMP2/MCQDPT2/XMCQDPT2 use the diagonal-Fock (MP-like) zeroth order only.",
                value=h0_mode,
                expected="fock (or leave unset)",
                action="Drop [pt2] h0 for QDPT-family methods, or use the caspt2 methods for h0=dyall.",
            )

    engine_mode = str(_get(config, "pt2", "engine", "auto") or "auto").strip().lower()
    if engine_mode not in {"auto", "direct", "python", "fortran", "dense"}:
        report.add(
            "ERROR",
            "pt2.engine",
            "Unknown PT2 engine.",
            value=engine_mode,
            expected="auto, fortran, direct (python), or dense",
            action="Use engine=auto (matrix-free NumPy streaming), fortran "
                   "(liboqp hash kernel, explicit opt-in), direct, or dense.",
        )
    elif engine_mode in {"direct", "python", "fortran"} \
            and canonical_method not in PT2_QDPT_METHODS:
        report.add(
            "ERROR",
            "pt2.engine",
            "The matrix-free direct engine requires a diagonal H0 (QDPT family).",
            value={"method": canonical_method, "engine": engine_mode},
            expected="method=mrmp2/mcqdpt2/xmcqdpt2 for engine=direct/fortran",
            action="Use the QDPT-family methods, or engine=dense for CASPT2/NEVPT2.",
        )
    _check_int_literal(
        _get(config, "pt2", "max_terms", 30000000),
        "pt2.max_terms",
        report,
        minimum=1,
        expected="positive streamed-term guard",
        action="Set [pt2] max_terms to a positive integer.",
    )
    _check_int_literal(
        _get(config, "pt2", "nproc", 1),
        "pt2.nproc",
        report,
        minimum=1,
        expected="positive worker count",
        action="Set [pt2] nproc to 1 (serial) or the number of worker processes.",
    )

    if reference not in PT2_REFERENCES:
        report.add(
            "ERROR",
            "pt2.reference",
            "Unknown PT2 reference type.",
            value=reference,
            expected=", ".join(sorted(PT2_REFERENCES)),
            action="Use reference=casci for fixed-orbital CI or reference=casscf after native CASSCF exists.",
        )
    elif reference_report and method in {"casci", "casscf"} and reference != method:
        report.add(
            "ERROR",
            "pt2.reference",
            "CASPT2 reference preflight source must match the native reference method.",
            value={"method": method, "reference": reference},
            expected=method,
            action=f"Set [pt2] reference={method} for this reference preflight.",
            wiki=WIKI_HELP["input.method"],
        )

    if contraction not in PT2_CONTRACTIONS:
        report.add(
            "ERROR",
            "pt2.contraction",
            "Unknown PT2 contraction mode.",
            value=contraction,
            expected=", ".join(sorted(PT2_CONTRACTIONS)),
            action="Use contraction=uncontracted; internally contracted CASPT2 is not implemented in this native helper.",
        )

    if contraction in PT2_STRONG_CONTRACTIONS:
        h0_mode = str(_get(config, "pt2", "h0", "fock") or "fock").strip().lower()
        if h0_mode not in {"dyall", "nevpt2", "nevpt"}:
            report.add(
                "ERROR",
                "pt2.contraction",
                "Strong contraction (SC-NEVPT2) is defined for Dyall's H0 only.",
                value=f"contraction={contraction}, h0={h0_mode}",
                expected="h0=dyall",
                action="Set [pt2] h0=dyall together with contraction=strong.",
            )

    if multistate not in PT2_MULTISTATE_MODES:
        report.add(
            "ERROR",
            "pt2.multistate",
            "Unknown PT2 multistate mode.",
            value=multistate,
            expected=", ".join(sorted(PT2_MULTISTATE_MODES)),
            action="Use auto, none, ms, or xms.",
        )
        normalized_multistate = "auto"
    elif multistate == "auto":
        if normalized_variant in PT2_SINGLE_STATE_METHODS:
            normalized_multistate = "none"
        elif normalized_variant in PT2_MS_METHODS:
            normalized_multistate = "ms"
        else:
            normalized_multistate = "xms"
    else:
        normalized_multistate = multistate
    xms = normalized_variant in PT2_XMS_METHODS or xms_requested

    if normalized_variant in PT2_SINGLE_STATE_METHODS and normalized_multistate != "none":
        report.add(
            "ERROR",
            "pt2.multistate",
            "Single-state PT2 (caspt2/mrmp2) cannot request an MS/XMS multistate mode.",
            value=multistate,
            expected="auto or none",
            action="Use method=ms-caspt2/mcqdpt2 or xms-caspt2/xmcqdpt2 for multistate PT2.",
        )
    if normalized_variant in PT2_MS_METHODS and normalized_multistate != "ms":
        report.add(
            "ERROR",
            "pt2.multistate",
            "MS-CASPT2/MCQDPT2 require multistate=ms.",
            value=multistate,
            expected="auto or ms",
            action="Use method=caspt2/mrmp2 for single-state PT2 or method=xms-caspt2/xmcqdpt2 for the extended variants.",
        )
    if normalized_variant in PT2_XMS_METHODS:
        if normalized_multistate != "xms":
            report.add(
                "ERROR",
                "pt2.multistate",
                "XMS-CASPT2/XMCQDPT2 require multistate=xms.",
                value=multistate,
                expected="auto or xms",
                action="Use method=ms-caspt2/mcqdpt2 for the non-extended multistate variants.",
            )
    elif xms_requested:
        report.add(
            "ERROR",
            "pt2.xms",
            "xms=true is only consistent with method=xms-caspt2 or xmcqdpt2.",
            value=xms,
            expected=False,
            action="Use method=xms-caspt2 or xmcqdpt2 for the extended multistate variants.",
        )

    nroot_valid, nroot_int = _check_int_literal(
        nroot,
        "pt2.nroot",
        report,
        minimum=0,
        expected="0 or a positive integer",
        action="Use 0 to inherit ci.nroot or set a positive PT2 root count.",
    )
    if not nroot_valid:
        nroot_int = 0

    try:
        ci_nroot_int = int(ci_nroot)
    except (TypeError, ValueError):
        ci_nroot_int = 1
    retained_roots = nroot_int or ci_nroot_int
    if retained_roots < 1:
        report.add(
            "ERROR",
            "pt2.nroot",
            "PT2 requires at least one retained reference root.",
            value=retained_roots,
            expected=">= 1",
            action="Increase [ci] nroot or [pt2] nroot.",
        )
    if normalized_variant in (PT2_MS_METHODS | PT2_XMS_METHODS) and retained_roots < 2:
        report.add(
            "ERROR",
            "pt2.nroot",
            "Multistate PT2 requires at least two retained reference roots.",
            value=retained_roots,
            expected=">= 2",
            action="Set [ci] nroot=2 or higher, or set [pt2] nroot=2 or higher.",
        )

    parsed_target_roots = []
    for root in target_roots:
        valid_root, parsed_root = _check_int_literal(
            root,
            "pt2.target_roots",
            report,
            minimum=0,
            expected="0-based PT2 root slots",
            action="Use values such as 0,1 for the first two retained roots.",
        )
        if not valid_root:
            parsed_target_roots = []
            break
        parsed_target_roots.append(parsed_root)
    if len(set(parsed_target_roots)) != len(parsed_target_roots):
        report.add(
            "ERROR",
            "pt2.target_roots",
            "PT2 target roots must be unique.",
            value=target_roots,
            action="Remove duplicate PT2 root slots.",
        )
    if parsed_target_roots and max(parsed_target_roots) >= retained_roots:
        report.add(
            "ERROR",
            "pt2.target_roots",
            "PT2 target roots must refer to retained PT2 root slots.",
            value=target_roots,
            expected=f"all roots < pt2.nroot ({retained_roots})",
            action="Increase [pt2] nroot or lower [pt2] target_roots.",
        )

    scalar_nonnegative = {
        "pt2.ipea_shift": ipea_shift,
        "pt2.imaginary_shift": imaginary_shift,
        "pt2.level_shift": level_shift,
        "pt2.edshft": edshft,
        "pt2.intruder_threshold": intruder_threshold,
        "pt2.amplitude_threshold": amplitude_threshold,
    }
    parsed_scalars = {}
    for path, value in scalar_nonnegative.items():
        valid_numeric, numeric = _check_float_literal(value, path, report)
        if not valid_numeric:
            continue
        parsed_scalars[path] = numeric
        if numeric < 0.0:
            report.add(
                "ERROR",
                path,
                "PT2 scalar option cannot be negative.",
                value=value,
                expected=">= 0",
                action="Use zero or a positive value.",
            )
    if parsed_scalars.get("pt2.edshft", 0.0) and (
        parsed_scalars.get("pt2.level_shift", 0.0)
        or parsed_scalars.get("pt2.imaginary_shift", 0.0)
    ):
        report.add(
            "ERROR",
            "pt2.edshft",
            "The ISA denominator shift (edshft, GAMESS EDSHFT) is mutually exclusive "
            "with the real/imaginary level shifts.",
            value={"edshft": edshft, "level_shift": level_shift,
                   "imaginary_shift": imaginary_shift},
            expected="only one regularization scheme",
            action="Use either [pt2] edshft OR level_shift/imaginary_shift, not both.",
        )

    valid_cutoff, denominator_cutoff_float = _check_float_literal(
        denominator_cutoff,
        "pt2.denominator_cutoff",
        report,
        expected="positive finite number",
        action="Use a small positive threshold.",
    )
    if not valid_cutoff:
        denominator_cutoff_float = 1.0
    if valid_cutoff and denominator_cutoff_float <= 0.0:
        report.add(
            "ERROR",
            "pt2.denominator_cutoff",
            "PT2 denominator cutoff must be positive and finite.",
            value=denominator_cutoff,
            expected="> 0",
            action="Use a small positive threshold.",
        )

    _check_int_literal(
        max_memory,
        "pt2.max_memory",
        report,
        minimum=1,
        expected="positive integer MiB",
        action="Set [pt2] max_memory to a positive MiB value.",
    )

    if not semi_canonical:
        report.add(
            "ERROR",
            "pt2.semi_canonical",
            "Native CASPT2 requires semi-canonical inactive/active/virtual orbitals.",
            value=semi_canonical,
            expected=True,
            action="Keep [pt2] semi_canonical=true until noncanonical PT2 is implemented.",
        )

    if isinstance(benchmark_reference_file_raw, bool) or not isinstance(benchmark_reference_file_raw, str):
        report.add(
            "ERROR",
            "pt2.benchmark_reference_file",
            "PT2 benchmark reference file must be a string path to a same-setup JSON report.",
            value=benchmark_reference_file_raw,
            expected="string path or empty string",
            action="Use a JSON reference report path generated from the same CASPT2 setup.",
        )
        benchmark_reference_file = ""
    else:
        benchmark_reference_file = benchmark_reference_file_raw.strip()

    valid_benchmark_tolerance, benchmark_tolerance_float = _check_float_literal(
        benchmark_tolerance,
        "pt2.benchmark_tolerance",
        report,
        expected="finite non-negative number",
        action="Use an absolute energy tolerance such as 1.0e-4.",
    )
    if valid_benchmark_tolerance and benchmark_tolerance_float < 0.0:
        report.add(
            "ERROR",
            "pt2.benchmark_tolerance",
            "PT2 benchmark tolerance cannot be negative.",
            value=benchmark_tolerance,
            expected=">= 0",
            action="Use zero for exact matching or a positive absolute energy tolerance.",
        )

    if benchmark_required and not benchmark_reference_file:
        report.add(
            "ERROR",
            "pt2.benchmark_reference_file",
            "PT2 benchmark_required needs an explicit same-setup reference report.",
            value=benchmark_reference_file_raw,
            expected="non-empty JSON reference report path",
            action="Set [pt2] benchmark_reference_file or disable benchmark_required.",
        )

    if pt2:
        report.add(
            "INFO",
            "pt2",
            (
                "XMS-CASPT2 options were parsed and validated for native determinant-space extended multistate PT2 dispatch."
                if canonical_method == "xms-caspt2"
                else (
                    "MS-CASPT2 options were parsed and validated for native determinant-space multistate PT2 dispatch."
                    if canonical_method == "ms-caspt2"
                    else "CASPT2 options were parsed and validated for native determinant-space PT2 dispatch."
                )
            ),
            value=pt2,
            action=(
                "Use method=xms-caspt2 for native XMS state rotation plus determinant-space multistate PT2."
                if canonical_method == "xms-caspt2"
                else (
                    "Use method=ms-caspt2 for native multistate PT2."
                    if canonical_method == "ms-caspt2"
                    else "Use method=caspt2 for diagonal state-specific native PT2."
                )
            ),
        )


def _check_planned_wavefunction_methods(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    if method not in PLANNED_WAVEFUNCTION_METHODS:
        return

    canonical = {
        "sacasscf": "sa-casscf",
        "mscaspt2": "ms-caspt2",
        "xmscaspt2": "xms-caspt2",
    }.get(method, method)
    report.add(
        "ERROR",
        "input.method",
        f"method={canonical} is recognized but not implemented in this milestone.",
        value=method,
        expected="casci, fci, casscf, sa-casscf, caspt2, ms-caspt2, or xms-caspt2",
        action="Use method=caspt2, method=ms-caspt2, or method=xms-caspt2 for native determinant-space PT2.",
        wiki=WIKI_HELP["input.method"],
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

    if method in {
        "fci", "casci", "casscf", "sa-casscf", "sacasscf",
        "caspt2", "ms-caspt2", "mscaspt2", "xms-caspt2", "xmscaspt2",
    } and scf_prop:
        method_label = method.upper()
        report.add(
            "WARNING",
            "properties.scf_prop",
            f"SCF properties requested with method={method} are RHF reference-density properties, not {method_label}-density properties.",
            value=scf_prop,
            action=f"Set [properties] scf_prop empty for a pure {method_label} energy job, or interpret these as RHF-reference diagnostics only.",
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
    if runtype in {"optimize", "mep", "ts", "irc", "neb"}:
        requested.append(_get(config, "optimize", "istate", 0))
    if runtype in {"meci", "mecp"}:
        search = _as_lower(_get(config, "optimize", "meci_search", "auto"))
        multistates = _as_list(_get(config, "optimize", "states", []))
        if runtype == "meci" and search in {"auto", "baeka"} and multistates:
            requested.extend(multistates)
        else:
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

    if method == "mp2" and runtype != "energy":
        report.add(
            "ERROR",
            "input.runtype",
            "MP2 currently supports energy-only calculations.",
            value=runtype,
            expected="energy",
            action="Use runtype=energy until MP2 gradients and derivative workflows are implemented.",
            wiki=WIKI_HELP["input.method"],
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

    if runtype == "namd":
        td_type = _as_lower(_get(config, "tdhf", "type", "rpa"))
        if method not in {"tdhf", "dftb", "xtb"} or td_type != "mrsf":
            report.add(
                "ERROR",
                "input.runtype",
                "NAMD (surface hopping) currently supports only MRSF-TDDFT.",
                value=f"{method}/{td_type}",
                expected="input.method=tdhf and tdhf.type=mrsf",
                action="Set [input] method=tdhf and [tdhf] type=mrsf for NAMD.",
                wiki=WIKI_HELP["tdhf.type"],
            )
        nstate = int(_get(config, "tdhf", "nstate", 1))
        active = int(_get(config, "md", "active", 1))
        soc_val = _get(config, "md", "soc", False)
        soc = (soc_val is True) or (str(soc_val).lower() in ("true", "1", "on", "yes"))
        if soc:
            # SOC-NAMD hops on the spin-adiabatic manifold: ns singlets +
            # 3*nt triplet Ms sublevels (ns = nt = tdhf.nstate).
            amax = nstate + 3 * nstate
            if active < 1 or active > amax:
                report.add(
                    "ERROR",
                    "md.active",
                    "Initial active state must be a valid spin-adiabatic state.",
                    value=active,
                    expected=f"1 <= md.active <= ns+3*nt ({amax}) for SOC-NAMD",
                    action="Set [md] active within the spin-adiabatic manifold (4*nstate states), or raise [tdhf] nstate.",
                )
        elif active < 1 or active > nstate:
            report.add(
                "ERROR",
                "md.active",
                "Initial active state must be a valid excited state.",
                value=active,
                expected=f"1 <= md.active <= tdhf.nstate ({nstate})",
                action="Set [md] active within the number of excited states, and raise [tdhf] nstate if needed.",
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

    if method in ("dftb", "xtb") and runtype in {"nac", "bp"}:
        # numerical NAC vectors / branching-plane are not wired for DFTB; the
        # DFTB runtype gate in _check_dftb already reports the unsupported case.
        # Still validate the NAC worker count so an invalid [nac] nproc is caught
        # in preflight rather than crashing multiprocessing.Pool(processes=0).
        if runtype == "nac":
            nac_nproc = _get(config, "nac", "nproc", 1)
            try:
                nac_nproc_int = int(nac_nproc)
            except (TypeError, ValueError):
                nac_nproc_int = 0
            if nac_nproc_int < 1:
                report.add(
                    "ERROR",
                    "nac.nproc",
                    "Number of NAC workers must be positive.",
                    value=nac_nproc,
                    expected="nac.nproc >= 1",
                    action="Set [nac] nproc to a positive integer.",
                )
        return
    if method in ("dftb", "xtb") and runtype == "nacme":
        # DFTB nacme IS wired (TLF state-overlap backend), but it still needs the
        # previous-step geometry; run only that portion of the nacme validation.
        _check_dftb_nacme_previous_geometry(config, report)
        return

    if runtype == "neb":
        _check_neb(config, report, input_dir)

    if runtype in {"nac", "bp"}:
        _check_nac(config, report)

    if runtype == "nacme":
        _check_nacme(config, report)

    if runtype == "soc":
        _check_soc(config, report)

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
    meci_search = _as_lower(_get(config, "optimize", "meci_search", "auto"))
    meci_states = _as_list(_get(config, "optimize", "states", []))

    if bool(_get(config, "input", "qmmm_flag", False)):
        report.add(
            "ERROR",
            "input.qmmm_flag",
            "Geometry and reaction-path drivers are not connected to the active QM/MM force backend.",
            value=f"qmmm_flag=true/runtype={runtype}",
            expected="a supported QM/MM energy, md, or namd workflow",
            action="Disable qmmm_flag for this geometry job; do not run a gas-phase optimizer on embedded coordinates.",
        )

    if lib not in OPT_LIBS:
        report.add(
            "ERROR",
            "optimize.lib",
            "Unknown optimization library.",
            value=lib,
            expected=", ".join(sorted(OPT_LIBS)),
            action="Use oqp (recommended), scipy, or geometric.",
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

    if lib == "oqp":
        auto_recovery = _get(config, "oqp", "auto_recovery", True)
        if not isinstance(auto_recovery, bool):
            report.add(
                "ERROR", "oqp.auto_recovery",
                "Native optimizer recovery must be a boolean.",
                value=auto_recovery, expected="true or false",
                action="Set [oqp] auto_recovery=true to enable the safe restart.",
            )
        recovery_maxit = _get(config, "oqp", "recovery_maxit", 30)
        try:
            valid_recovery_maxit = int(recovery_maxit) > 0
        except (TypeError, ValueError):
            valid_recovery_maxit = False
        if not valid_recovery_maxit:
            report.add(
                "ERROR", "oqp.recovery_maxit",
                "Native optimizer recovery needs a positive step budget.",
                value=recovery_maxit, expected="a positive integer",
                action="Set [oqp] recovery_maxit=30 or another positive integer.",
            )
        recovery_trust = _get(config, "oqp", "recovery_trust", 0.02)
        try:
            valid_recovery_trust = (
                math.isfinite(float(recovery_trust))
                and float(recovery_trust) > 0.0
            )
        except (TypeError, ValueError):
            valid_recovery_trust = False
        if not valid_recovery_trust:
            report.add(
                "ERROR", "oqp.recovery_trust",
                "Native optimizer recovery trust radius must be positive.",
                value=recovery_trust, expected="a finite positive number",
                action="Set [oqp] recovery_trust=0.02.",
            )

        freeze = str(_get(config, "oqp", "freeze", "") or "").strip()
        if freeze:
            if runtype != "optimize":
                report.add(
                    "ERROR", "oqp.freeze",
                    "Frozen-distance constraints are currently supported only for a native minimum search.",
                    value=runtype, expected="runtype=optimize",
                    action="Use a separate opt job for the constrained minimum.",
                )
            pattern = re.compile(
                r"(?:distance|r)\(\s*(\d+)\s*,\s*(\d+)\s*\)", re.I
            )
            items = [item.strip() for item in freeze.split(";") if item.strip()]
            valid = bool(items)
            for item in items:
                match = pattern.fullmatch(item)
                if not match or int(match.group(1)) < 1 or int(match.group(2)) < 1 \
                        or int(match.group(1)) == int(match.group(2)):
                    valid = False
                    break
            if not valid:
                report.add(
                    "ERROR", "oqp.freeze",
                    "Invalid native frozen-distance expression.",
                    value=freeze, expected="distance(1,2);distance(2,3)",
                    action="Use one-based distinct atom pairs separated by semicolons.",
                )

        init_hessian = _as_lower(_get(config, "oqp", "init_hessian", "model"))
        if init_hessian not in {"model", "numerical", "analytical"}:
            report.add(
                "ERROR", "oqp.init_hessian",
                "Unknown native initial-Hessian policy.",
                value=init_hessian,
                expected="model, numerical, or analytical",
                action="Use model unless a fresh numerical or analytical TS Hessian is required.",
            )
        elif init_hessian != "model" and runtype != "ts":
            report.add(
                "ERROR", "oqp.init_hessian",
                "A real native initial Hessian is currently consumed only by runtype=ts.",
                value=f"{runtype}/{init_hessian}",
                expected="runtype=ts",
                action="Use oqp.init_hessian=model or switch to a native TS search.",
            )
        elif init_hessian == "analytical":
            capability, reason = analytic_hessian_capability(config)
            if capability != "supported":
                report.add(
                    "ERROR", "oqp.init_hessian",
                    "Analytical native TS initialization is unavailable for this method/state.",
                    value=init_hessian,
                    expected="a supported HF/DFT ground-state analytic Hessian",
                    action=reason,
                )
            max_l = _basis_max_angular_momentum(config)
            if max_l is not None and max_l >= 4:
                report.add(
                    "ERROR", "input.basis",
                    "Analytical native TS initialization supports basis angular momentum only up to L=3.",
                    value=f"max L={max_l}",
                    expected="max L <= 3",
                    action="Use hessian=numerical for this TS search, or choose a basis without g/higher functions.",
                )

        hess_restart = _get(config, "hess", "restart", False)
        if (runtype == "irc" and hess_restart) or (
                runtype == "ts" and init_hessian != "model" and hess_restart):
            report.add(
                "ERROR", "hess.restart",
                "Native TS/IRC Hessian restart artifacts are not yet signed to a geometry and electronic model.",
                value=hess_restart,
                expected="False",
                action="Set [hess] restart=False; standalone numerical Hessian restart remains available.",
            )

        irc_direction = _as_lower(_get(config, "oqp", "irc_direction", "forward"))
        if runtype == "irc" and irc_direction not in {"forward", "backward", "reverse"}:
            report.add(
                "ERROR", "oqp.irc_direction",
                "Unknown native IRC direction.",
                value=irc_direction,
                expected="forward or backward",
                action="Set irc_direction=forward or irc_direction=backward.",
            )

        if runtype in {"irc", "mep"}:
            path_gtol = _get(config, "oqp", "path_gtol", 1.0e-4)
            try:
                path_gtol_value = float(path_gtol)
                valid_path_gtol = math.isfinite(path_gtol_value) and path_gtol_value > 0
            except (TypeError, ValueError):
                valid_path_gtol = False
            if not valid_path_gtol:
                report.add(
                    "ERROR", "oqp.path_gtol",
                    "Native path gradient tolerance must be positive.",
                    value=path_gtol,
                    expected="> 0",
                    action="Set path_gtol to a positive value such as 1e-4.",
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

    if method == "tdhf" and runtype in {"optimize", "mep", "ts", "irc", "neb"} and istate == 0:
        report.add(
            "ERROR",
            "optimize.istate",
            "TDHF optimization cannot target state 0 in the current implementation.",
            value=istate,
            expected=">= 1",
            action="Set istate to a TDHF excited state index.",
        )

    if runtype == "meci":
        if method not in {"tdhf", "dftb", "xtb"}:
            report.add(
                "ERROR",
                "input.method",
                "MECI optimization requires an excited-state method.",
                value=method,
                expected="tdhf, dftb, or xtb",
                action="Set [input] method=tdhf, dftb, or xtb and configure the response state space.",
            )
        if (meci_search not in {"auto", "baeka"} or not meci_states) \
                and jstate <= istate:
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
                action="Use auto (recommended), penalty, ubp, hybrid, or baeka.",
            )
        elif meci_search in {"auto", "baeka"}:
            selected = meci_states or [istate, jstate]
            try:
                selected = [int(state) for state in selected]
            except (TypeError, ValueError):
                selected = []
            if len(selected) < 2:
                report.add(
                    "ERROR", "optimize.states",
                    "BaekA MECI requires at least two response roots.",
                    value=meci_states,
                    expected="two or more increasing roots",
                    action="Set states=1,2 or a longer consecutive list.",
                )
            elif any(state < 1 for state in selected):
                report.add(
                    "ERROR", "optimize.states",
                    "BaekA states must be positive response roots; root 0 is the internal reference.",
                    value=selected,
                    expected="1,2 or 1,2,3,...",
                    action="Use positive response-root indices only.",
                )
            elif any(right != left + 1
                     for left, right in zip(selected, selected[1:])):
                report.add(
                    "ERROR", "optimize.states",
                    "BaekA MECI requires distinct, increasing, consecutive response roots.",
                    value=selected,
                    expected="1,2 or 1,2,3,...",
                    action="Sort the roots and include every intervening state.",
                )
            if lib != "oqp" and (
                    meci_search == "baeka" or len(selected) > 2):
                report.add(
                    "ERROR", "optimize.lib",
                    "BaekA multistate MECI is implemented by the native oqp optimizer.",
                    value=lib, expected="oqp",
                    action="Set [optimize] lib=oqp.",
                )

            scalar_controls = {
                "pen_sigma": (1.0, "initial sigma"),
                "pen_alpha": (0.0, "alpha in Hartree"),
                "pen_delta": (0.025, "additive delta_beta"),
                "gap_weight": (1.0, "gap weight"),
                "energy_gap": (1.0e-4, "outer-span gap tolerance"),
            }
            for key, (default, label) in scalar_controls.items():
                raw = _get(config, "optimize", key, default)
                try:
                    value = float(raw)
                    if key == "gap_weight":
                        valid = math.isfinite(value) and value == 1.0
                    else:
                        valid = math.isfinite(value) and (
                            value > 0 or (key == "pen_alpha" and value == 0)
                        )
                except (TypeError, ValueError):
                    valid = False
                if not valid:
                    if key == "gap_weight":
                        report.add(
                            "ERROR", "optimize.gap_weight",
                            "BaekA fixes gap_weight at 1.0; sigma is the published penalty multiplier.",
                            value=raw, expected="1.0",
                            action="Set gap_weight=1.0 and adjust pen_sigma if needed.",
                        )
                        continue
                    report.add(
                        "ERROR", f"optimize.{key}",
                        f"BaekA {label} must be positive.",
                        value=raw, expected="> 0",
                        action=("Use pen_alpha=0.02; alpha has energy units. "
                                "The legacy zero sentinel also selects 0.02 for BaekA."
                                if key == "pen_alpha" else "Use a positive value."),
                    )

            jumps = _get(
                config, "optimize", "pen_jump",
                [10, 10, 25, 25, 100, 100, 1000, 1000, 3000],
            )
            if isinstance(jumps, str):
                jumps = [item.strip() for item in jumps.split(",") if item.strip()]
            elif not isinstance(jumps, (list, tuple)):
                jumps = [jumps]
            try:
                jumps = [float(item) for item in jumps]
                valid_jumps = bool(jumps) and all(
                    math.isfinite(item) and item > 0 for item in jumps
                )
            except (TypeError, ValueError):
                valid_jumps = False
            if not valid_jumps:
                report.add(
                    "ERROR", "optimize.pen_jump",
                    "BaekA beta jump schedule must contain positive values.",
                    value=jumps, expected="10,10,25,...",
                    action="Provide a non-empty comma-separated positive schedule.",
                )
        elif meci_states:
            report.add(
                "ERROR", "optimize.states",
                "The multistate roots list is used only by meci_search=baeka.",
                value=meci_states,
                expected="meci_search=baeka",
                action="Remove states or select the BaekA algorithm.",
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
            expected="oqp or geometric",
            action="Use [optimize] lib=oqp (recommended) or lib=geometric for runtype=ts/irc.",
        )

    if runtype == "neb" and lib not in {"geometric", "oqp"}:
        report.add(
            "ERROR",
            "optimize.lib",
            "NEB is wired through the native oqp optimizer and the optional geomeTRIC backend.",
            value=lib,
            expected="oqp or geometric",
            action="Set [optimize] lib=oqp (recommended) or lib=geometric for runtype=neb.",
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
    lib = _as_lower(_get(config, "optimize", "lib", "oqp"))
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

    if lib == "oqp":
        # Keep sectioned legacy inputs as strict as the concise .oqp parser.
        # These controls are consumed only by the native NEB implementation.
        numeric_controls = {
            "spring": (0.05, 0.0, True),
            "fmax": (2.0e-3, 0.0, False),
            "frms": (2.0e-3, 0.0, False),
            "climb_fmax": (0.05, 0.0, False),
            "neb_dt": (0.5, 0.0, False),
            "maxmove": (0.2, 0.0, False),
            "end_fmax": (1.0e-3, 0.0, False),
        }
        for key, (default, lower, inclusive) in numeric_controls.items():
            value = _get(config, "oqp", key, default)
            try:
                number = float(value)
                valid = math.isfinite(number) and (
                    number >= lower if inclusive else number > lower
                )
            except (TypeError, ValueError):
                valid = False
            if not valid:
                relation = ">=" if inclusive else ">"
                report.add(
                    "ERROR",
                    f"oqp.{key}",
                    "Invalid native NEB control.",
                    value=value,
                    expected=f"a finite number {relation} {lower:g}",
                    action=f"Set [oqp] {key} to a valid native NEB value.",
                )

        for key in ("climb", "align", "opt_ends"):
            value = _get(config, "oqp", key, True)
            if not isinstance(value, bool):
                report.add(
                    "ERROR",
                    f"oqp.{key}",
                    "Invalid native NEB switch.",
                    value=value,
                    expected="true or false",
                    action=f"Set [oqp] {key}=true or false.",
                )

        climb = _get(config, "oqp", "climb", True)
        try:
            fmax = float(_get(config, "oqp", "fmax", 2.0e-3))
            climb_fmax = float(_get(config, "oqp", "climb_fmax", 0.05))
        except (TypeError, ValueError):
            fmax = climb_fmax = float("nan")
        if (isinstance(climb, bool) and climb
                and math.isfinite(fmax) and math.isfinite(climb_fmax)
                and climb_fmax < fmax):
            report.add(
                "ERROR",
                "oqp.climb_fmax",
                "Climbing-image activation must occur before final NEB convergence.",
                value=climb_fmax,
                expected="climb_fmax >= fmax when climb=true",
                action="Increase climb_fmax or disable climbing-image NEB.",
            )


def _check_soc(config: dict[str, Any], report: CheckReport) -> None:
    method = _as_lower(_get(config, "input", "method", "hf"))
    td_type = _as_lower(_get(config, "tdhf", "type", "rpa"))
    scf_type = _as_lower(_get(config, "scf", "type", "rhf"))
    scf_mult = _get(config, "scf", "multiplicity", 1)

    if method in ("dftb", "xtb"):
        # MRSF one-center SOC: the ROKS triplet reference is implied
        # by the TB backend (the adapter enforces it), so only the response
        # type needs checking here.
        if td_type != "mrsf":
            report.add(
                "ERROR",
                "tdhf.type",
                "SOC requires tdhf.type=mrsf.",
                value=td_type,
                expected="mrsf",
                action="Set [tdhf] type=mrsf.",
            )
        return

    if method != "tdhf":
        report.add(
            "ERROR",
            "input.method",
            "SOC requires method=tdhf.",
            value=method,
            expected="tdhf",
            action="Set [input] method=tdhf.",
        )

    if td_type != "mrsf":
        report.add(
            "ERROR",
            "tdhf.type",
            "SOC requires tdhf.type=mrsf.",
            value=td_type,
            expected="mrsf",
            action="Set [tdhf] type=mrsf.",
        )

    if scf_type != "rohf":
        report.add(
            "ERROR",
            "scf.type",
            "SOC requires an ROHF reference.",
            value=scf_type,
            expected="rohf",
            action="Set [scf] type=rohf.",
        )

    if scf_mult != 3:
        report.add(
            "ERROR",
            "scf.multiplicity",
            "SOC requires a triplet ROHF reference (multiplicity=3).",
            value=scf_mult,
            expected="3",
            action="Set [scf] multiplicity=3.",
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

    if hess_type not in {"numerical", "analytical"}:
        report.add(
            "ERROR",
            "hess.type",
            "Unknown Hessian type.",
            value=hess_type,
            expected="numerical or analytical",
            action="Set [hess] type=numerical or type=analytical.",
        )
        return

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
    nac_type = _as_lower(_get(config, "nac", "type", "numerical"))

    if nac_type != "numerical":
        report.add(
            "ERROR",
            "nac.type",
            "Analytical NAC vectors are not available.",
            value=nac_type,
            expected="numerical",
            action="Set [nac] type=numerical.",
        )

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


def _check_dftb_nacme_previous_geometry(config: dict[str, Any], report: CheckReport) -> None:
    """DFTB nacme reuses the previous-geometry validation of _check_nacme
    (presence, system2 path existence, non-empty inline geometry, nac.align)
    without the all-electron method=tdhf gate of _check_nac."""
    guess_file2 = _get(config, "guess", "file2", "")
    system2 = _get(config, "input", "system2", "")

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

    align = _as_lower(_get(config, "nac", "align", "reorder"))
    if align not in {"no", "reorder"}:
        report.add(
            "WARNING",
            "nac.align",
            "Only align=no and align=reorder are clearly handled by the Python layer.",
            value=align,
            action="Use align=reorder unless you have implemented another mode downstream.",
        )


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
    # Expand top-level DFTB response shorthands (method=mrsf-tddftb / sf-tddftb /
    # tddftb) to the canonical method=dftb form before any validation runs.
    expand_dftb_method_alias(config)
    method = _as_lower(_get(config, "input", "method", "hf"))

    if method not in METHODS and method not in PLANNED_WAVEFUNCTION_METHODS:
        report.add(
            "ERROR",
            "input.method",
            "Unknown electronic structure method.",
            value=method,
            expected=", ".join(sorted(METHODS)),
            action="Choose hf, tdhf, mp2, dftb, or xtb (or the DFTB shortcuts "
                   "mrsf-tddftb, sf-tddftb, tddftb), or a wavefunction method: "
                   "fci, casci, casscf, sa-casscf, caspt2, ms-caspt2, xms-caspt2, "
                   "mrmp2, mcqdpt2, xmcqdpt2.",
            wiki=WIKI_HELP["input.method"],
        )

    _check_system(config, report)
    _check_basis(config, report)
    _check_guess(config, report)
    _check_pcm(config, report)
    _check_dftb(config, report)
    _check_xtb(config, report)
    _check_scf(config, report)
    _check_symmetry(config, report)
    _check_tdhf(config, report)
    _check_mp2(config, report)
    _check_fci(config, report)
    _check_casci(config, report)
    _check_casscf(config, report)
    _check_pt2(config, report)
    _check_planned_wavefunction_methods(config, report)
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
