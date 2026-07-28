"""Semantic ``.oqp`` input parser and deterministic request compiler.

The existing INI-style ``.inp`` reader remains the execution format used by
the native/Python layers.  This module adds a small, user-facing language that
describes *physical* methods and states and lowers them to that legacy format.
In particular, users request ``mrsf ... opt(S0)``; they do not need to know that
OpenQP represents the MRSF singlet ground state as response root 1 on top of a
triplet ROHF reference.

This module deliberately uses only the Python standard library.  It can be
loaded before :mod:`oqp` initializes the native library, which is useful for
CLI pre-processing and lightweight editor/web validation.
"""

from __future__ import annotations

import ast
import difflib
import json
import math
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


class OQPInputError(ValueError):
    """Raised when semantic input is ambiguous, contradictory, or invalid."""


@dataclass(frozen=True)
class StateRef:
    """A physical state label, or an explicit implementation root.

    ``root=N`` is intentionally available for SF-TDDFT, whose spin-adapted
    physical state labels cannot be known reliably before diagonalization.
    """

    label: Optional[str] = None
    root: Optional[int] = None

    @property
    def multiplicity(self) -> Optional[int]:
        if self.label is None:
            return None
        return {"S": 1, "T": 3, "Q": 5}[self.label[0]]

    @property
    def physical_index(self) -> Optional[int]:
        return int(self.label[1:]) if self.label is not None else None

    def __str__(self) -> str:
        return self.label if self.label is not None else "root=%d" % self.root


@dataclass(frozen=True)
class CallSpec:
    """A driver or section-style call from canonical input."""

    name: str
    args: Tuple[Any, ...] = ()
    kwargs: Mapping[str, Any] = field(default_factory=dict)
    explicit: bool = True


@dataclass(frozen=True)
class CalculationSpec:
    """Parsed, implementation-independent calculation request."""

    model: str
    functional: str
    basis: str
    model_options: Mapping[str, Any]
    options: Mapping[str, Any]
    driver: CallSpec
    modifiers: Tuple[CallSpec, ...]
    source_text: str = ""

    @property
    def reference_method(self) -> str:
        if self.model in {
            "mrsf", "mrsf-hf", "mrsf-dftb", "sf", "sf-hf", "sf-dftb"
        }:
            return "ROHF"
        if self.model in {"umrsf", "umrsf-hf"}:
            return "UHF"
        if self.model in {"rks", "uks", "roks"}:
            return {"rks": "RHF", "uks": "UHF", "roks": "ROHF"}[self.model]
        if self.model in {"rhf", "uhf", "rohf"}:
            return self.model.upper()
        if self.model == "mp2" and "reference" in self.model_options:
            return str(self.model_options["reference"]).upper()
        mult = int(self.options.get("mult", 1))
        return "RHF" if mult == 1 else "UHF"

    @property
    def reference_multiplicity(self) -> int:
        if self.model in {
            "mrsf", "mrsf-hf", "mrsf-dftb", "umrsf", "umrsf-hf", "sf", "sf-hf",
            "sf-dftb",
        }:
            return 3
        return int(self.options.get("mult", 1))

    @property
    def physical_method(self) -> str:
        labels = {
            "mrsf": "MRSF-TDDFT" if self.functional else "MRSF-TDHF",
            "umrsf": "UMRSF-TDDFT" if self.functional else "UMRSF-TDHF",
            "sf": "SF-TDDFT" if self.functional else "SF-TDHF",
            "tddft": "TDDFT",
            "tda": "TDA-TDDFT" if self.functional else "TDA-TDHF",
            "tdhf": "TDHF",
            "dft": "DFT",
            "rks": "RKS-DFT",
            "uks": "UKS-DFT",
            "roks": "ROKS-DFT",
            "hf": "HF",
            "rhf": "RHF",
            "uhf": "UHF",
            "rohf": "ROHF",
            "mp2": "MP2",
            "dftb": "DFTB",
            "dftb0": "DFTB0",
            "tddftb": "TD-DFTB (TDA)",
            "tda-dftb": "TD-DFTB (TDA)",
            "sf-dftb": "SF-TDDFTB",
            "mrsf-dftb": "MRSF-TDDFTB",
            "mrsf-hf": "MRSF-TDHF",
            "umrsf-hf": "UMRSF-TDHF",
            "sf-hf": "SF-TDHF",
            "tda-hf": "TDA-TDHF (CIS)",
        }
        return labels[self.model]


@dataclass(frozen=True)
class OQPResolution:
    """Result of classifying, compiling, reparsing, and lowering an input."""

    spec: CalculationSpec
    canonical_text: str
    legacy_config: Mapping[str, Mapping[str, str]]
    was_natural: bool
    source_path: Optional[Path] = None
    resolved_path: Optional[Path] = None


def _normalize_default_section_call(call: CallSpec) -> Optional[CallSpec]:
    """Remove concise section options whose runtime defaults say the same thing."""

    if call.name == "afqmc":
        aliases = {
            "output": "output_dir",
            "walkers": "nwalkers",
            "steps": "nsteps",
            "dt": "timestep",
            "threads": "omp_threads",
            "chol": "chol_tol",
        }
        kwargs = {}
        for key, value in call.kwargs.items():
            canonical = aliases.get(key, key)
            if canonical in kwargs:
                raise OQPInputError("Duplicate AFQMC option: %s" % canonical)
            kwargs[canonical] = value
        return CallSpec(call.name, call.args, kwargs, call.explicit)
    if call.name != "dftb":
        return call
    kwargs = dict(call.kwargs)
    if str(kwargs.get("backend", "")).strip().lower() == "native":
        kwargs.pop("backend")
    if "parameter_path" in kwargs and not str(kwargs["parameter_path"]).strip():
        kwargs.pop("parameter_path")
    if not call.args and not kwargs:
        return None
    return CallSpec(call.name, call.args, kwargs, call.explicit)


def _normalize_geometry_owned_section_call(
    call: CallSpec,
    geometry: Any,
) -> CallSpec:
    """Remove section values already supplied by the compact geometry."""

    if call.name != "qmmm" or "pdb_file" not in call.kwargs:
        return call
    geometry_pdb = _qmmm_pdb_from_geometry(geometry, source_dir=None)
    explicit_pdb = call.kwargs["pdb_file"]
    if not isinstance(geometry_pdb, str) or not isinstance(explicit_pdb, str):
        return call
    normalized_geometry = os.path.normcase(
        os.path.normpath(os.path.expanduser(geometry_pdb.strip()))
    )
    normalized_explicit = os.path.normcase(
        os.path.normpath(os.path.expanduser(explicit_pdb.strip()))
    )
    if normalized_geometry != normalized_explicit:
        return call
    kwargs = dict(call.kwargs)
    kwargs.pop("pdb_file")
    return CallSpec(call.name, call.args, kwargs, call.explicit)


DEFAULT_SINGLET_MODELS = {
    "dft", "rks", "uks", "roks", "hf", "rhf", "uhf", "rohf", "mp2",
    "tddft", "tda", "tdhf", "tda-hf",
}


# Calls in this set are mutually exclusive.  Everything else is a modifier or
# a direct legacy-section call.  Aliases are normalized before counting, so
# ``grad(...) opt(...)`` is always rejected rather than silently picking one.
PRIMARY_ALIASES = {
    "energy": "energy",
    "sp": "energy",
    "grad": "grad",
    "gradient": "grad",
    "opt": "optimize",
    "optimize": "optimize",
    "optimization": "optimize",
    "meci": "meci",
    "mecp": "mecp",
    "tci": "tci",
    "mep": "mep",
    "ts": "ts",
    "irc": "irc",
    "neb": "neb",
    "hess": "hess",
    "hessian": "hess",
    "nac": "nac",
    "bp": "bp",
    "nacme": "nacme",
    "soc": "soc",
    "md": "md",
    "namd": "namd",
    "ekt": "ekt",
    "thermo": "thermo",
    "prop": "prop",
    "data": "data",
}

# These modifiers have a complete, useful meaning without arguments.  Accept
# ``nmr`` and ``pcm`` just as we accept ``opt`` and ``soc``; rendering still
# has one deterministic spelling for every accepted form.
BARE_MODIFIER_CALLS = {"pcm", "nmr", "ir", "raman", "d4"}

SECTION_NAMES = {
    "input", "afqmc", "mp2", "guess", "pcm", "dftb", "symmetry", "scf",
    "dftgrid", "tdhf", "ekt", "properties", "optimize", "geometric",
    "oqp", "neb", "hess", "nac", "md", "qmmm", "json", "tests",
}


def _keys(words: str) -> frozenset[str]:
    """Build a compact, immutable keyword vocabulary from plain text."""

    return frozenset(words.split())


# Static mirror of OQP_CONFIG_SCHEMA, partitioned by ownership.  This module is
# intentionally stdlib-only and can be imported before the native ``oqp``
# package, so it cannot import oqpdata.py directly.  A pure-Python CI test
# compares this manifest with the authoritative schema via AST and fails when a
# schema keyword is added without an explicit owner here.
GENERIC_SCHEMA_KEYS = {
    "afqmc": _keys("""
        output_dir chol_tol trial trial_file nwalkers nsteps timestep seed
        omp_threads stabilize_every population_control_every estimate_every
        accumulate_after force_bias_bound
    """),
    "qmmm": _keys("""
        forcefield nonbondedmethod constraints rigidwater nsteps n_steps
        timestep pdb_file forcefield_files qm_atoms cutoff embedding
        temperature ensemble friction pressure barostat_interval
        trajectory_format trajectory_file log_file report_interval energy_file
        qm_atoms_xyz qm_list frontier_scheme
    """),
    "input": _keys("library perf ispher d4 qmmm_flag soc_2e omp_threads"),
    "mp2": _keys("variant same_spin_scale opposite_spin_scale"),
    "guess": _keys("type file file2 save_mol continue_geom swapmo"),
    "pcm": _keys("enabled backend mode model solvent epsilon radii"),
    "dftb": _keys("""
        backend parameter_path library_path executable scc_tolerance scc_mixer
        scc_mixing scc_history scc_max_step max_scc_iterations
        response_tolerance response_max_iterations response_max_subspace
        response_solver print_level state_to_state_spectrum
        spc spc_coco spc_ovov spc_coov omega cam_alpha
        cam_beta lc_gamma lc_ground_state zvector spin_complete mrsf_shift_oo
        mrsf_shift_co mrsf_shift_ov mrsf_shift_cv timeout
        model c_mrsf c_mrsf_oo response_global_hybrid onsite_exchange_scale
        w_scale response_w_scale response_omega response_cam_alpha
        response_cam_beta onsite_ss onsite_sp onsite_pp
    """),
    "symmetry": _keys("""
        enabled point_group subgroup label_mo label_states label_modes
        use_integral_symmetry use_response_symmetry tolerance strict
    """),
    "scf": _keys("""
        maxit forced_attempt maxdiis diis_reset_mod diis_reset_conv diis_type
        cdiis_switch vdiis_vshift_switch vshift mom mom_switch pfon
        pfon_start_temp pfon_cooling_rate pfon_nsmear conv incremental pscreen
        pscreen_k pscreen_cap pscreen_tight pscreen_xc_dcut pscreen_xc_aocut
        pscreen_grid_rad pscreen_grid_ang xc_c2f xc_phi_cache xc_incdft
        grad_cutoff init_scf init_basis init_library init_it init_conv
        init_converger save_molden rstctmo converger_type scal_rel stability
        soscf_lvl_shift alternative_scf escalation verbose trh_stab trh_ls
        trh_sub_solver trh_nrtv trh_r0 trh_jd_start trh_nmic trh_gred
        trh_lred trh_impl
    """),
    "dftgrid": _keys("""
        hfscale cam_flag cam_alpha cam_beta cam_mu rad_type rad_npts ang_npts
        partfun pruned grid_ao_pruned grid_ao_threshold grid_ao_sparsity_ratio
    """),
    "tdhf": _keys("""
        maxit maxit_zv conv nstate_s nstate_t zvconv nvdav tlf hfscale
        cam_alpha cam_beta cam_mu spc_coco spc_ovov spc_coov conf_threshold
        ixcore z_solver gmres_dim resp_cutoff fp32 zv_warmstart
    """),
    "ekt": frozenset(),
    "properties": _keys("scf_prop nmr_gauge td_prop nac export title back_door"),
    "optimize": frozenset(),
    "geometric": frozenset(),
    "oqp": frozenset(),
    "neb": frozenset(),
    "hess": frozenset(),
    "nac": frozenset(),
    "md": frozenset(),
    "json": _keys("scf_type basis library do_init"),
    "tests": _keys("exception"),
}

ROUTE_DRIVER_SCHEMA_KEYS = {
    "input": _keys("charge basis functional method runtype system system2"),
    "dftb": _keys("type reference_multiplicity target_multiplicity"),
    "scf": _keys("type multiplicity"),
    "tdhf": _keys("type multiplicity nstate target"),
    "ekt": _keys("ip ea"),
    "properties": _keys("grad"),
    "optimize": _keys("""
        lib maxit rmsd_grad rmsd_step max_grad max_step istate jstate kstate states
        imult jmult energy_shift energy_gap meci_search pen_sigma pen_alpha
        pen_incre pen_delta pen_jump gap_weight init_scf
    """),
    "neb": _keys("product nimage"),
    "oqp": _keys("""
        coordsys trust trust_max auto_recovery recovery_maxit recovery_trust
        freeze follow init_hessian spring climb fmax frms climb_fmax neb_dt
        maxmove align opt_ends end_fmax neb_output irc_step irc_direction
        mep_step path_gtol
    """),
    "hess": _keys("type state dx nproc read restart temperature clean"),
    "nac": _keys("type dt dx bp nproc restart clean states align"),
    "md": _keys("""
        nstep dt active substep decoherence edc_c thrshe tdc trivial
        trivial_thresh init_temp velocity seed restart soc soc_basis
        soc_du_dt_corr soc_tdc_grad_corr grad_wthr init_state econs
        dt_adaptive dt_min dx_max
    """),
}

# Traditional sectioned ``.inp`` files and the Python workflow API retain
# backend selection and geomeTRIC/SciPy controls.  The concise ``.oqp`` format
# intentionally hides those implementation details and always selects the
# native OpenQP geometry engine.  Keeping these keys in the inventory preserves
# the legacy schema without exposing them as canonical top-level options.
LEGACY_ONLY_SCHEMA_KEYS = {
    "optimize": _keys("optimizer step_size step_tol mep_maxit"),
    "geometric": _keys("""
        coordsys trust tmax convergence_set prefix hessian irc_direction
        constraints_file enforce conmethod
    """),
    "neb": _keys("k maxg avgg climb align optep"),
}

INTENTIONALLY_FORBIDDEN_SCHEMA_KEYS = {
    # Retained in the legacy schema for the disconnected libopenmm driver.
    # Active qmmm_md/NAMD paths do not read this numeric selector; physical
    # state-aware workflows use their normal driver state instead.
    "qmmm": _keys("istate"),
}

SCHEMA_KEY_OWNERS: Dict[str, Dict[str, str]] = {name: {} for name in SECTION_NAMES}
for _owner, _manifest in (
    ("generic", GENERIC_SCHEMA_KEYS),
    ("route_driver", ROUTE_DRIVER_SCHEMA_KEYS),
    ("legacy_only", LEGACY_ONLY_SCHEMA_KEYS),
    ("intentional_forbidden", INTENTIONALLY_FORBIDDEN_SCHEMA_KEYS),
):
    for _section, _section_keys in _manifest.items():
        if _section not in SCHEMA_KEY_OWNERS:
            raise RuntimeError("Unknown schema section in .oqp manifest: %s" % _section)
        for _key in _section_keys:
            if _key in SCHEMA_KEY_OWNERS[_section]:
                raise RuntimeError("Duplicate .oqp owner for %s.%s" % (_section, _key))
            SCHEMA_KEY_OWNERS[_section][_key] = _owner

OQP_SCHEMA_KEYS = {
    section: frozenset(owners) for section, owners in SCHEMA_KEY_OWNERS.items()
}

MODEL_ALIASES = {
    "mrsf": "mrsf",
    "mrsf-tddft": "mrsf",
    "mrsftddft": "mrsf",
    "mrsf-tdhf": "mrsf-hf",
    "umrsf": "umrsf",
    "umrsf-tddft": "umrsf",
    "umrsftddft": "umrsf",
    "umrsf-tdhf": "umrsf-hf",
    "sf": "sf",
    "sf-tddft": "sf",
    "sftddft": "sf",
    "sf-tdhf": "sf-hf",
    "tddft": "tddft",
    "td-dft": "tddft",
    "tda": "tda",
    "tda-tdhf": "tda-hf",
    "cis": "tda-hf",
    "tdhf": "tdhf",
    "dft": "dft",
    "ks-dft": "dft",
    "rks": "rks",
    "uks": "uks",
    "roks": "roks",
    "hf": "hf",
    "rhf": "rhf",
    "uhf": "uhf",
    "rohf": "rohf",
    "mp2": "mp2",
    "dftb": "dftb",
    "dftb0": "dftb0",
    "dftb-noscc": "dftb0",
    "dftb-nonscc": "dftb0",
    "tddftb": "tddftb",
    "td-dftb": "tddftb",
    "tda-tddftb": "tda-dftb",
    "tda-dftb": "tda-dftb",
    "tddftb-tda": "tda-dftb",
    "sf-tddftb": "sf-dftb",
    "sf-dftb": "sf-dftb",
    "sftddftb": "sf-dftb",
    "mrsf-tddftb": "mrsf-dftb",
    "mrsf-dftb": "mrsf-dftb",
}

TOP_OPTION_ALIASES = {
    "geom": "geom",
    "geometry": "geom",
    "system": "geom",
    "geom2": "geom2",
    "geometry2": "geom2",
    "system2": "geom2",
    "charge": "charge",
    "mult": "mult",
    "multiplicity": "mult",
    "basis": "basis",
    "functional": "functional",
    "library": "library",
    "ispher": "ispher",
    "perf": "perf",
    "d4": "d4",
    "qmmm": "qmmm_flag",
    "qmmm_flag": "qmmm_flag",
    "omp_threads": "omp_threads",
    "threads": "omp_threads",
}

# Public compact-input driver signatures. State selectors (positional S0/S1/T0 or
# ``root=N`` for SF) are handled separately; the sets below are the concise
# workflow options that may live directly in the primary call.
_OPT_OPTIONS = {
    "maxit", "rmsd_grad", "rmsd_step", "max_grad", "max_step",
    "energy_shift", "energy_gap", "meci_search", "pen_sigma",
    "pen_alpha", "pen_incre", "pen_delta", "pen_jump", "gap_weight", "init_scf",
}
_NATIVE_ENGINE_OPTIONS = {
    "coordsys", "trust", "trust_max",
    "auto_recovery", "recovery_maxit", "recovery_trust",
}
_NATIVE_CONSTRAINT_OPTIONS = {"freeze"}
_GEOMETRIC_ENGINE_OPTIONS = {
    "coordsys", "trust", "tmax", "convergence_set", "hessian",
}
_OPTIMIZER_BACKEND_OPTIONS = (
    {"lib"} | _NATIVE_ENGINE_OPTIONS | _NATIVE_CONSTRAINT_OPTIONS
    | _GEOMETRIC_ENGINE_OPTIONS
)
_MECI_PUBLIC_OPTIONS = {
    "algorithm", "sigma", "alpha", "delta_beta", "beta_schedule", "gap",
}
_MECI_OPTION_ALIASES = {
    "algorithm": "meci_search",
    "sigma": "pen_sigma",
    "alpha": "pen_alpha",
    "delta_beta": "pen_delta",
    "beta_schedule": "pen_jump",
    "gap": "energy_gap",
}
_TCI_OPTIONS = set(_OPT_OPTIONS) - {"meci_search", "pen_delta", "pen_jump"}
DRIVER_OPTIONS = {
    "energy": set(),
    "grad": {"td_prop", "export", "title"},
    "optimize": set(_OPT_OPTIONS) | set(_OPTIMIZER_BACKEND_OPTIONS),
    "meci": set(_OPT_OPTIONS) | set(_MECI_PUBLIC_OPTIONS) | set(_NATIVE_ENGINE_OPTIONS),
    "mecp": set(_OPT_OPTIONS) | set(_NATIVE_ENGINE_OPTIONS),
    "tci": set(_TCI_OPTIONS) | set(_NATIVE_ENGINE_OPTIONS),
    "mep": {"maxit", "points", "step", "mep_step", "gtol"},
    "ts": set(_OPT_OPTIONS) | set(_NATIVE_ENGINE_OPTIONS) | {"follow", "hessian"},
    "irc": {"maxit", "direction", "step", "irc_step", "hessian", "gtol"},
    "neb": {
        "maxit",
        "product", "images", "nimage", "spring", "climb", "fmax",
        "frms", "climb_fmax", "dt", "neb_dt", "maxmove", "align",
        "opt_ends", "end_fmax", "output",
    },
    "hess": {"type", "dx", "nproc", "read", "restart", "temperature", "clean"},
    "nac": {"type", "dx", "nproc", "restart", "clean", "align"},
    "bp": {"type", "dx", "nproc", "restart", "clean", "align"},
    "nacme": {"dt", "align"},
    "soc": {"soc_2e", "ns", "nt"},
    # QM/MM MD controls live in qmmm(...).  The similarly named [md]
    # controls below belong to the nonadiabatic-dynamics driver only.
    "md": set(),
    "namd": {
        "nstep", "dt", "active", "substep", "decoherence", "edc_c",
        "thrshe", "tdc", "trivial", "trivial_thresh", "init_temp",
        "velocity", "seed", "restart", "soc", "soc_basis",
        "soc_du_dt_corr", "soc_tdc_grad_corr", "grad_wthr", "init_state",
        "econs", "dt_adaptive", "dt_min", "dx_max",
    },
    "ekt": {"ip", "ea"},
    "thermo": {"type", "dx", "nproc", "read", "restart", "temperature", "clean"},
    "prop": {"scf_prop", "nmr_gauge", "td_prop", "export", "title"},
    "data": {"scf_prop", "nmr_gauge", "td_prop", "export", "title"},
}

# Exact ``oqp(...)`` section calls are accepted only when every key is consumed
# by the selected native workflow.  This prevents valid-looking but ignored
# controls such as ``energy oqp(init_hessian=analytical)``.
OQP_DRIVER_OPTIONS = {
    "optimize": set(_NATIVE_ENGINE_OPTIONS) | set(_NATIVE_CONSTRAINT_OPTIONS),
    "meci": set(_NATIVE_ENGINE_OPTIONS),
    "mecp": set(_NATIVE_ENGINE_OPTIONS),
    "tci": set(_NATIVE_ENGINE_OPTIONS),
    "ts": set(_NATIVE_ENGINE_OPTIONS) | {"follow", "init_hessian"},
    "mep": {"mep_step", "path_gtol"},
    "irc": {"irc_step", "irc_direction", "path_gtol"},
    "neb": {
        "spring", "climb", "fmax", "frms", "climb_fmax", "neb_dt",
        "maxmove", "align", "opt_ends", "end_fmax", "neb_output",
    },
}

# Route parentheses describe the electronic model, not an alternate spelling
# for every legacy section.  Keeping this surface deliberately small prevents
# route options from silently overriding section-owned controls.  All advanced
# settings remain available through exact calls such as ``tdhf(nvdav=30)``.
RESPONSE_MODEL_OPTIONS = {"nstate"}
MP2_MODEL_OPTIONS = {
    "reference", "variant", "same_spin_scale", "opposite_spin_scale",
}
RESPONSE_MODELS = {
    "mrsf", "mrsf-hf", "mrsf-dftb", "umrsf", "umrsf-hf", "sf", "sf-hf",
    "sf-dftb", "tddft", "tda", "tda-hf", "tdhf", "tddftb", "tda-dftb"
}
SF_STATE_AWARE_DRIVERS = {
    "grad", "optimize", "mep", "ts", "irc", "neb", "hess", "thermo",
    "md", "namd", "data",
}

ACTIVE_QMMM_DRIVERS = {"energy", "md", "namd"}

_STATE_RE = re.compile(r"^([STQ])(\d+)$", re.IGNORECASE)
_IDENT_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_-]*$")


def _is_integer(value: Any) -> bool:
    """Return whether *value* is an actual integer, excluding booleans."""

    return isinstance(value, int) and not isinstance(value, bool)


def _is_positive_integer(value: Any) -> bool:
    """Return whether a parsed value is an actual positive integer.

    ``int(value)`` is intentionally not used here: it would silently accept
    values such as ``2.5``, ``"2"``, or ``true`` in count fields.
    """

    return _is_integer(value) and value > 0


def _is_non_negative_integer(value: Any) -> bool:
    """Return whether *value* is an integer greater than or equal to zero."""

    return _is_integer(value) and value >= 0


def _strip_comments(text: str) -> str:
    """Remove ``#``/``!`` comments outside quotes.

    A marker is never required.  Supporting comments remains convenient for
    generated and hand-written multi-line files.
    """

    out: List[str] = []
    for line in text.splitlines():
        quote: Optional[str] = None
        escaped = False
        kept: List[str] = []
        for index, char in enumerate(line):
            if escaped:
                kept.append(char)
                escaped = False
                continue
            if char == "\\" and quote:
                kept.append(char)
                escaped = True
                continue
            if char in {'"', "'"}:
                if quote == char:
                    quote = None
                elif quote is None:
                    quote = char
                kept.append(char)
                continue
            if quote is None and char in {"#", "!"} and (
                index == 0 or line[index - 1].isspace()
            ):
                break
            kept.append(char)
        out.append("".join(kept))
    return "\n".join(out)


def _compact_syntactic_whitespace(text: str) -> str:
    """Accept harmless spaces around ``=`` and before a call's ``(``.

    Quoted values are copied byte-for-byte, so a geometry path or title may
    contain the same characters without being rewritten.
    """

    out: List[str] = []
    plain: List[str] = []
    quote: Optional[str] = None
    escaped = False

    def flush_plain() -> None:
        segment = "".join(plain)
        plain.clear()
        segment = re.sub(r"\s*=\s*", "=", segment)
        segment = re.sub(
            r"\b([A-Za-z_][A-Za-z0-9_-]*)\s+\(", r"\1(", segment
        )
        out.append(segment)

    for char in text:
        if quote is None:
            if char in {'"', "'"}:
                flush_plain()
                quote = char
                out.append(char)
            else:
                plain.append(char)
            continue
        out.append(char)
        if escaped:
            escaped = False
        elif char == "\\":
            escaped = True
        elif char == quote:
            quote = None
    flush_plain()
    return "".join(out)


def _split_top_level(text: str, delimiter: Optional[str] = None) -> List[str]:
    """Split on a delimiter (or whitespace) outside calls and quotes."""

    result: List[str] = []
    buf: List[str] = []
    depth = 0
    quote: Optional[str] = None
    escaped = False
    for char in text:
        if escaped:
            buf.append(char)
            escaped = False
            continue
        if char == "\\" and quote:
            buf.append(char)
            escaped = True
            continue
        if char in {'"', "'"}:
            if quote == char:
                quote = None
            elif quote is None:
                quote = char
            buf.append(char)
            continue
        if quote is None:
            if char == "(":
                depth += 1
            elif char == ")":
                depth -= 1
                if depth < 0:
                    raise OQPInputError("Unmatched ')' in canonical .oqp input")
            should_split = (
                depth == 0
                and ((delimiter is None and char.isspace()) or char == delimiter)
            )
            if should_split:
                item = "".join(buf).strip()
                if item:
                    result.append(item)
                buf = []
                continue
        buf.append(char)
    if quote:
        raise OQPInputError("Unterminated quoted value in canonical .oqp input")
    if depth:
        raise OQPInputError("Unclosed '(' in canonical .oqp input")
    item = "".join(buf).strip()
    if item:
        result.append(item)
    return result


def _split_assignment(text: str) -> Optional[Tuple[str, str]]:
    pieces = _split_top_level(text, "=")
    if len(pieces) == 1:
        return None
    if len(pieces) != 2:
        raise OQPInputError("Invalid assignment: %s" % text)
    return pieces[0].strip(), pieces[1].strip()


def _parse_value(text: str) -> Any:
    text = text.strip()
    if not text:
        raise OQPInputError("Missing value after '='")
    state = _STATE_RE.fullmatch(text)
    if state:
        return StateRef(label=state.group(1).upper() + state.group(2))
    if text[0:1] in {'"', "'"}:
        try:
            value = ast.literal_eval(text)
        except (SyntaxError, ValueError) as exc:
            raise OQPInputError("Invalid quoted value: %s" % text) from exc
        if not isinstance(value, str):
            raise OQPInputError("Quoted .oqp values must be strings")
        return value
    lower = text.lower()
    if lower == "true":
        return True
    if lower == "false":
        return False
    if lower == "null":
        return None
    if re.fullmatch(r"[+-]?\d+", text):
        return int(text)
    if re.fullmatch(
        r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+-]?\d+)", text
    ) or re.fullmatch(r"[+-]?(?:\d+\.\d*|\.\d+)", text):
        return float(re.sub("[dD]", "e", text))
    return text


def _parse_call(text: str) -> CallSpec:
    match = re.fullmatch(r"([A-Za-z_][A-Za-z0-9_-]*)\s*\((.*)\)", text, re.DOTALL)
    if not match:
        if _IDENT_RE.fullmatch(text):
            name = text.lower().replace("-", "_")
            if name in PRIMARY_ALIASES:
                return CallSpec(PRIMARY_ALIASES[name])
            if name in BARE_MODIFIER_CALLS:
                return CallSpec(name)
        raise OQPInputError("Expected a call such as opt(S0): %s" % text)
    name = match.group(1).lower().replace("-", "_")
    raw_args = _split_top_level(match.group(2), ",")
    args: List[Any] = []
    kwargs: Dict[str, Any] = {}
    keyword_seen = False
    for raw in raw_args:
        assignment = _split_assignment(raw)
        if assignment is None:
            if keyword_seen:
                raise OQPInputError(
                    "Positional arguments must precede keyword arguments in %s" % name
                )
            args.append(_parse_value(raw))
            continue
        keyword_seen = True
        key, value = assignment
        key = key.strip().lower().replace("-", "_")
        if not _IDENT_RE.fullmatch(key):
            raise OQPInputError("Invalid keyword name in %s: %s" % (name, key))
        if key in kwargs:
            raise OQPInputError("Duplicate keyword in %s: %s" % (name, key))
        kwargs[key] = _parse_value(value)
    return CallSpec(name=name, args=tuple(args), kwargs=kwargs)


def _parse_route(route: str) -> Tuple[str, Dict[str, Any], str, str]:
    parts = _split_top_level(route, "/")
    if not parts or len(parts) > 3:
        raise OQPInputError(
            "Route must be model[/functional][/basis], got: %s" % route
        )
    model_call = _parse_call(parts[0]) if "(" in parts[0] else CallSpec(parts[0])
    alias = model_call.name.lower().replace("_", "-")
    if alias not in MODEL_ALIASES:
        close = difflib.get_close_matches(alias, MODEL_ALIASES, n=1, cutoff=0.6)
        hint = " Did you mean '%s'?" % close[0] if close else ""
        raise OQPInputError(
            "Unknown electronic-structure model: %s.%s" % (model_call.name, hint)
        )
    if model_call.args:
        raise OQPInputError("Model options must use names, for example mrsf(nstate=5)")
    model = MODEL_ALIASES[alias]
    model_options = dict(model_call.kwargs)
    allowed_model_options = (
        RESPONSE_MODEL_OPTIONS
        if model in RESPONSE_MODELS
        else MP2_MODEL_OPTIONS if model == "mp2" else set()
    )
    unknown_model_options = set(model_options) - allowed_model_options
    if unknown_model_options:
        key = sorted(unknown_model_options)[0]
        if model in {
            "mrsf", "mrsf-hf", "mrsf-dftb", "umrsf", "umrsf-hf"
        } and key in {
            "mult", "multiplicity", "spin"
        }:
            raise OQPInputError(
                "%s(...) does not take %s. Select the target manifold in the driver, "
                "for example opt(S0) or opt(T0)." % (alias, key)
            )
        section = "mp2" if model == "mp2" else "tdhf" if allowed_model_options else model
        raise OQPInputError(
            "%s is not a route option for %s; use the exact section call %s(%s=...)"
            % (key, alias, section, key)
        )
    functional = ""
    basis = ""
    if model in {
        "dftb", "dftb0", "tddftb", "tda-dftb", "sf-dftb", "mrsf-dftb"
    } and len(parts) != 1:
        raise OQPInputError("%s does not take a functional or basis route component" % alias)
    if model in {
        "hf", "rhf", "uhf", "rohf", "mp2", "tdhf", "mrsf-hf",
        "umrsf-hf", "sf-hf", "tda-hf",
    } and len(parts) == 3:
        raise OQPInputError("%s route is model/basis; it does not take a functional" % alias)
    if len(parts) == 2:
        if model in {
            "hf", "rhf", "uhf", "rohf", "mp2", "tdhf", "mrsf-hf",
            "umrsf-hf", "sf-hf", "tda-hf", "dftb", "dftb0", "tddftb",
            "tda-dftb", "sf-dftb", "mrsf-dftb",
        }:
            basis = parts[1]
        else:
            functional = parts[1]
    elif len(parts) == 3:
        functional, basis = parts[1], parts[2]
    if model in {"dft", "rks", "uks", "roks", "tddft", "tda", "mrsf", "umrsf", "sf"} and not functional:
        raise OQPInputError("%s requires a functional in the route" % model)
    return model, model_options, functional.lower(), basis.lower()


def looks_canonical(text: str) -> bool:
    """Return whether malformed text should be treated as canonical, not prose."""

    raw = text.lstrip()
    if raw.startswith("#") and "/" in raw.splitlines()[0]:
        # A Gaussian-style route marker is a canonical-input mistake, not
        # prose to be reinterpreted by the optional correction helper.
        return True
    content = _strip_comments(text).lstrip()
    # A user may reuse a valid route prefix while describing the rest of the
    # job in prose (for example ``... h2o.xyz에서 triplet energy
    # calculation``).  File-name particles and complete action phrases are
    # strong prose signals; quoted Korean file names alone are not.
    if re.search(r"\.(?:xyz|pdb)(?:에서|를|을|으로|로)", content, re.I):
        return False
    if re.search(
        r"\b(?:energy calculation|calculate\s+(?:the\s+)?energy|"
        r"geometry optimization|optimi[sz]e\s+(?:the\s+)?geometry)\b",
        content,
        re.I,
    ):
        return False
    stripped = content.strip()
    if not stripped:
        return False
    if stripped.startswith("["):
        return True
    first = stripped.split(None, 1)[0]
    first_name = re.split(r"[/()]", first, maxsplit=1)[0].lower().replace("_", "-")
    if first_name in MODEL_ALIASES:
        return True
    if "/" in first and ("(" in first or re.match(r"^[A-Za-z][\w-]*/", first)):
        return True
    if re.search(r"\b(?:geom|geometry|charge|basis|functional)\s*=", stripped, re.I):
        return True
    if re.search(r"\b(?:energy|grad|opt|meci|mecp|neb|hess)\s*\(", stripped, re.I):
        return True
    return False


def parse_canonical_oqp(text: str) -> CalculationSpec:
    """Parse and semantically validate markerless canonical ``.oqp`` text."""

    cleaned = _compact_syntactic_whitespace(_strip_comments(text)).strip()
    if not cleaned:
        single = text.strip()
        if single.startswith("#") and "/" in single:
            raise OQPInputError(
                "A leading '#' route marker is not used in .oqp. Remove '#'; "
                "the route may be followed by options and calls on the same or later lines."
            )
        raise OQPInputError("Empty .oqp input")
    if cleaned.startswith("["):
        raise OQPInputError(
            "Sectioned [input]/[scf] syntax belongs in a legacy .inp file. "
            "Use compact .oqp syntax such as dft/pbe0/def2-svp, energy, and "
            "geom=\"h2o.xyz\" on one or more lines."
        )
    tokens = _split_top_level(cleaned)
    model, model_options, functional, basis = _parse_route(tokens[0])
    options: Dict[str, Any] = {}
    calls: List[CallSpec] = []
    for token_index, token in enumerate(tokens[1:], start=1):
        assignment = _split_assignment(token)
        if token_index == 1 and assignment is None and "(" not in token:
            positional_geom = _parse_value(token)
            if (
                isinstance(positional_geom, str)
                and Path(positional_geom).suffix.lower() in {".xyz", ".pdb"}
            ):
                options["geom"] = positional_geom
                continue
        if assignment is not None and "(" not in token:
            key, value = assignment
            key = key.lower().replace("-", "_")
            if key not in TOP_OPTION_ALIASES:
                if key in {"istate", "jstate", "kstate"}:
                    raise OQPInputError(
                        "Raw %s is an internal legacy keyword. Put physical labels in the driver, "
                        "for example opt(S0) or meci(S1,S2)." % key
                    )
                close = difflib.get_close_matches(key, TOP_OPTION_ALIASES, n=1, cutoff=0.65)
                suggestion = TOP_OPTION_ALIASES[close[0]] if close else ""
                hint = " Did you mean '%s'?" % suggestion if suggestion else ""
                raise OQPInputError(
                    "Unknown top-level option '%s'; use a section call for advanced options.%s"
                    % (key, hint)
                )
            canonical_key = TOP_OPTION_ALIASES[key]
            if canonical_key in options:
                raise OQPInputError("Duplicate top-level option: %s" % canonical_key)
            parsed_value = _parse_value(value)
            if canonical_key in {"geom", "geom2"} and value.lstrip().startswith(
                ('"""', "'''")
            ):
                # Triple-quoted geometry blocks are formatted with the opening
                # and closing delimiter on their own lines. Those presentation
                # newlines are not part of the molecular coordinates.
                if parsed_value.startswith("\n"):
                    parsed_value = parsed_value[1:]
                if parsed_value.endswith("\n"):
                    parsed_value = parsed_value[:-1]
            options[canonical_key] = parsed_value
        else:
            calls.append(_parse_call(token))

    # Explicit route components win only if they agree with a flat spelling.
    if "functional" in options:
        flat = str(options.pop("functional")).lower()
        if functional and flat != functional:
            raise OQPInputError("Conflicting functional values in route and options")
        functional = flat
    if "basis" in options:
        flat = str(options.pop("basis")).lower()
        if basis and flat != basis:
            raise OQPInputError("Conflicting basis values in route and options")
        basis = flat
    if _is_integer(options.get("charge")) and options["charge"] == 0:
        options.pop("charge")
    if (
        model in DEFAULT_SINGLET_MODELS
        and _is_integer(options.get("mult"))
        and options["mult"] == 1
    ):
        options.pop("mult")

    primary: List[CallSpec] = []
    modifiers: List[CallSpec] = []
    for call in calls:
        normalized = call.name.lower().replace("-", "_")
        if normalized in PRIMARY_ALIASES:
            primary.append(CallSpec(
                PRIMARY_ALIASES[normalized], call.args, dict(call.kwargs), call.explicit
            ))
        elif normalized in SECTION_NAMES or normalized in {"nmr", "ir", "raman", "d4"}:
            normalized_call = _normalize_default_section_call(
                CallSpec(normalized, call.args, dict(call.kwargs))
            )
            if normalized_call is not None:
                modifiers.append(normalized_call)
        else:
            choices = list(PRIMARY_ALIASES) + list(SECTION_NAMES) + ["nmr", "ir", "raman", "d4"]
            close = difflib.get_close_matches(normalized, choices, n=1, cutoff=0.6)
            hint = " Did you mean '%s(...)'?" % close[0] if close else ""
            raise OQPInputError("Unknown .oqp call: %s.%s" % (call.name, hint))
    if len(primary) > 1:
        names = ", ".join(call.name for call in primary)
        raise OQPInputError(
            "Exactly one primary driver is allowed: choose one calculation such as "
            "energy, grad, or opt per .oqp file. Found %s; put each calculation in "
            "its own .oqp file."
            % names
        )
    driver = primary[0] if primary else CallSpec("energy", explicit=False)
    driver = _normalize_driver_defaults(model, driver)
    spec = CalculationSpec(
        model=model,
        functional=functional,
        basis=basis,
        model_options=model_options,
        options=options,
        driver=driver,
        modifiers=tuple(modifiers),
        source_text=text,
    )
    _validate_semantics(spec)
    return spec


def _state_from_value(value: Any, *, keyword: str = "state") -> StateRef:
    if isinstance(value, StateRef):
        return value
    if keyword == "root":
        if not _is_positive_integer(value):
            raise OQPInputError("root must be a positive integer")
        return StateRef(root=value)
    if isinstance(value, str) and _STATE_RE.fullmatch(value):
        return StateRef(label=value.upper())
    raise OQPInputError("Expected a physical state label such as S0, S1, or T0")


def _driver_states(driver: CallSpec) -> Tuple[StateRef, ...]:
    states: List[StateRef] = []
    for arg in driver.args:
        states.append(_state_from_value(arg))
    keys = ["state"]
    if driver.name in {"meci", "mecp", "nac", "bp", "nacme"}:
        keys.extend(("state1", "state2", "i", "j"))
    elif driver.name == "tci":
        keys.extend(("state1", "state2", "state3", "i", "j", "k"))
    for key in keys:
        if key in driver.kwargs:
            states.append(_state_from_value(driver.kwargs[key]))
    if "root" in driver.kwargs:
        states.append(_state_from_value(driver.kwargs["root"], keyword="root"))
    return tuple(states)


def _normalize_ground_state_driver(model: str, driver: CallSpec) -> CallSpec:
    """Make an explicit HF/DFT ``S0`` selector equal to the omitted default.

    Ground-state HF and DFT have only one physical surface, so ``grad(S0)`` and
    ``opt(S0)`` carry no more information than ``grad`` and ``opt``. Removing
    the redundant selector here gives both spellings the same canonical form,
    while leaving response-method state labels untouched.
    """

    if model not in {"dft", "rks", "uks", "roks", "hf", "rhf", "uhf", "rohf"}:
        return driver
    if driver.name not in {"grad", "optimize"}:
        return driver
    states = _driver_states(driver)
    if len(states) != 1 or states[0].label != "S0":
        return driver

    kwargs = dict(driver.kwargs)
    if driver.args:
        args: Tuple[Any, ...] = ()
    elif "state" in kwargs:
        kwargs.pop("state")
        args = driver.args
    else:
        return driver
    return CallSpec(driver.name, args, kwargs, driver.explicit)


def _normalize_driver_defaults(model: str, driver: CallSpec) -> CallSpec:
    """Remove driver options that exactly reproduce public runtime defaults."""

    driver = _normalize_ground_state_driver(model, driver)
    maxit = driver.kwargs.get("maxit")
    if (
        driver.name != "optimize"
        or not _is_integer(maxit)
        or maxit != 30
    ):
        return driver
    kwargs = dict(driver.kwargs)
    kwargs.pop("maxit")
    return CallSpec(driver.name, driver.args, kwargs, driver.explicit)


def _validate_semantics(spec: CalculationSpec) -> None:
    model = spec.model
    driver = spec.driver
    states = _driver_states(driver)
    if "geom" not in spec.options:
        raise OQPInputError(
            "Add a geometry after the route, for example h2o.xyz or geom=\"h2o.xyz\""
        )
    if not spec.basis and model not in {
        "dftb", "dftb0", "tddftb", "tda-dftb", "sf-dftb", "mrsf-dftb"
    }:
        raise OQPInputError("The route requires a basis set")

    if "charge" in spec.options and not _is_integer(spec.options["charge"]):
        raise OQPInputError("charge must be an integer")
    if "mult" in spec.options and not _is_positive_integer(spec.options["mult"]):
        raise OQPInputError("mult must be a positive integer")
    if "omp_threads" in spec.options and not _is_positive_integer(
            spec.options["omp_threads"]):
        raise OQPInputError("omp_threads must be a positive integer")

    if model in {"sf", "sf-hf", "sf-dftb"} and not states \
            and driver.name in SF_STATE_AWARE_DRIVERS:
        raise OQPInputError(
            "SF states cannot be labeled or inferred before diagonalization. "
            "Specify an implementation root explicitly, for example %s(root=1)."
            % ({"optimize": "opt"}.get(driver.name, driver.name))
        )

    if model in {
        "mrsf", "mrsf-hf", "mrsf-dftb", "umrsf", "umrsf-hf", "sf", "sf-hf",
        "sf-dftb",
    } and "mult" in spec.options:
        raise OQPInputError(
            "Top-level mult is not used for %s. The internal high-spin reference is "
            "automatic; select the target manifold with S0, T0, or Q0 in the driver."
            % model
        )
    if model == "rhf" and spec.options.get("mult", 1) != 1:
        raise OQPInputError("rhf requires mult=1; use rohf or uhf for an open-shell reference")
    if model == "rohf" and spec.options.get("mult", 1) < 2:
        raise OQPInputError("rohf requires an open-shell mult value such as mult=2 or mult=3")
    if model == "rks" and spec.options.get("mult", 1) != 1:
        raise OQPInputError("rks requires mult=1; use roks or uks for open-shell DFT")
    if model == "roks" and spec.options.get("mult", 1) < 2:
        raise OQPInputError("roks requires an open-shell mult value such as mult=2 or mult=3")
    if model == "mp2" and driver.name != "energy":
        raise OQPInputError("MP2 currently supports energy() only")
    if model == "mp2" and "reference" in spec.model_options:
        reference = str(spec.model_options["reference"]).strip().lower()
        if reference not in {"rhf", "rohf", "uhf"}:
            raise OQPInputError("MP2 reference must be rhf, rohf, or uhf")
        mult = spec.options.get("mult", 1)
        if reference == "rhf" and mult != 1:
            raise OQPInputError("MP2 reference=rhf requires mult=1")
        if reference == "rohf" and mult < 2:
            raise OQPInputError(
                "MP2 reference=rohf requires an open-shell mult value such as mult=2 or mult=3"
            )
    if model in {"umrsf", "umrsf-hf"} and driver.name != "energy":
        raise OQPInputError("UMRSF currently supports energy() only")
    if driver.name == "ekt" and model != "mrsf":
        raise OQPInputError("ekt currently requires an MRSF-TDDFT route")
    if driver.name in {"soc", "namd"} and model not in {
        "mrsf", "mrsf-hf", "mrsf-dftb"
    }:
        raise OQPInputError("%s currently requires an MRSF route" % driver.name)

    expected_counts = {"mecp": 2, "tci": 3, "nac": 2, "bp": 2, "nacme": 2}
    if driver.name == "meci" and len(states) < 2:
        raise OQPInputError("meci requires at least two physical states")
    if driver.name in expected_counts and len(states) != expected_counts[driver.name]:
        raise OQPInputError(
            "%s requires exactly %d physical states" % (driver.name, expected_counts[driver.name])
        )
    if driver.name not in expected_counts and driver.name != "meci" and len(states) > 1:
        raise OQPInputError("%s accepts at most one target state" % driver.name)
    if driver.name == "meci":
        roots_with_spin = [
            (state.multiplicity, _internal_root(model, state))
            for state in states
        ]
        if len(set(roots_with_spin)) != len(roots_with_spin):
            raise OQPInputError("meci requires distinct states")
    if driver.name == "energy" and states:
        if model not in RESPONSE_MODELS or any(
            state.label is None or state.physical_index != 0 for state in states
        ):
            raise OQPInputError(
                "energy() accepts only S0/T0/Q0 to select a response spin manifold"
            )
    if driver.name == "soc" and states:
        raise OQPInputError("%s does not take a target state" % driver.name)
    if driver.name == "prop" and model not in {"mrsf", "mrsf-hf"}:
        raise OQPInputError("prop([STATE]) currently requires MRSF-TDDFT or MRSF-TDHF")

    options = (_normalized_meci_options(driver)
               if driver.name == "meci" else _driver_options(driver))
    if driver.name == "meci":
        algorithm = str(options.get(
            "meci_search", "baeka" if len(states) > 2 else "auto"
        )).strip().lower()
        if algorithm not in {"auto", "penalty", "ubp", "hybrid", "baeka"}:
            raise OQPInputError(
                "meci algorithm must be auto, penalty, ubp, hybrid, or baeka"
            )
        if len(states) > 2 and algorithm not in {"auto", "baeka"}:
            raise OQPInputError(
                "MECI calculations with more than two states require "
                "algorithm=auto or algorithm=baeka"
            )
        if algorithm in {"auto", "baeka"}:
            roots = sorted(_internal_root(model, state) for state in states)
            if any(right != left + 1 for left, right in zip(roots, roots[1:])):
                raise OQPInputError(
                    "BaekA requires consecutive response roots in one spin manifold"
                )
            positive = {
                "pen_sigma": "sigma",
                "pen_alpha": "alpha",
                "pen_delta": "delta_beta",
                "energy_gap": "gap",
            }
            baeka_defaults = {
                "pen_sigma": 1.0,
                "pen_alpha": 0.02,
                "pen_delta": 0.025,
                "energy_gap": 1.0e-4,
            }
            for key, label in positive.items():
                try:
                    value = float(options.get(key, baeka_defaults[key]))
                    if not math.isfinite(value) or value <= 0:
                        raise ValueError
                except (TypeError, ValueError) as exc:
                    raise OQPInputError(
                        "BaekA %s must be a positive number" % label
                    ) from exc
            try:
                gap_weight = float(options.get("gap_weight", 1.0))
            except (TypeError, ValueError) as exc:
                raise OQPInputError(
                    "BaekA gap_weight is fixed at 1.0"
                ) from exc
            if not math.isfinite(gap_weight) or gap_weight != 1.0:
                raise OQPInputError(
                    "BaekA gap_weight is fixed at 1.0; use sigma for the "
                    "published penalty multiplier"
                )
            # ``auto`` starts with the historical penalty objective, where
            # pen_incre and the ordinary geometry convergence thresholds are
            # still meaningful.  They become invalid only for a direct BaekA
            # request, which has no conventional penalty phase.
            if algorithm == "baeka" and "pen_incre" in options:
                raise OQPInputError(
                    "BaekA uses additive delta_beta, not legacy pen_incre"
                )
            unused_convergence = (
                {"max_grad", "rmsd_step", "max_step"}.intersection(options)
                if algorithm == "baeka" else set()
            )
            if unused_convergence:
                raise OQPInputError(
                    "BaekA convergence does not use %s; remove it and use "
                    "energy_shift, gap, and rmsd_grad"
                    % ", ".join(sorted(unused_convergence))
                )
            schedule = options.get(
                "pen_jump", "10,10,25,25,100,100,1000,1000,3000"
            )
            if isinstance(schedule, str):
                schedule = [item.strip() for item in schedule.split(",") if item.strip()]
            elif not isinstance(schedule, (list, tuple)):
                schedule = [schedule]
            try:
                jumps = [float(item) for item in schedule]
            except (TypeError, ValueError) as exc:
                raise OQPInputError(
                    "BaekA beta_schedule must be a comma-separated list of positive numbers"
                ) from exc
            if not jumps or any(not math.isfinite(item) or item <= 0 for item in jumps):
                raise OQPInputError(
                    "BaekA beta_schedule must be a comma-separated list of positive numbers"
                )
        elif {"pen_delta", "pen_jump"}.intersection(options):
            raise OQPInputError(
                "delta_beta and beta_schedule are used only with "
                "algorithm=auto or algorithm=baeka"
            )
    legacy_backend_options = {
        "optimizer", "step_size", "step_tol", "mep_maxit",
        "k", "maxg", "avgg", "optep",
    }.intersection(options)
    if legacy_backend_options:
        key = sorted(legacy_backend_options)[0]
        if key in {"optimizer", "step_size", "step_tol", "mep_maxit"}:
            raise OQPInputError(
                "%s is a legacy SciPy-backend option. Concise .oqp workflows "
                "support lib=oqp or lib=geometric." % key
            )
        replacement = {
            "k": "spring", "maxg": "fmax", "avgg": "frms", "optep": "opt_ends",
        }[key]
        raise OQPInputError(
            "%s is a legacy geomeTRIC NEB option with different units or semantics; "
            "use the native .oqp option %s instead." % (key, replacement)
        )

    unknown_driver_options = set(options) - DRIVER_OPTIONS[driver.name]
    if unknown_driver_options:
        key = sorted(unknown_driver_options)[0]
        raise OQPInputError(
            "%s does not define option '%s'; use the exact legacy section call for advanced options"
            % (driver.name, key)
        )
    backend = str(options.get("lib", "oqp")).strip().lower()
    if driver.name == "optimize":
        if backend not in {"oqp", "geometric"}:
            raise OQPInputError(
                "opt lib must be oqp or geometric in concise .oqp input"
            )
        if backend == "geometric":
            native_only = {
                "trust_max", "auto_recovery", "recovery_maxit",
                "recovery_trust", "freeze",
            }.intersection(options)
            if native_only:
                raise OQPInputError(
                    "%s is a native lib=oqp option; use geomeTRIC tmax or a "
                    "traditional constraints file as appropriate"
                    % sorted(native_only)[0]
                )
        else:
            geometric_only = {
                "tmax", "convergence_set", "hessian",
            }.intersection(options)
            if geometric_only:
                raise OQPInputError(
                    "%s is available only with opt(...,lib=geometric)"
                    % sorted(geometric_only)[0]
                )
    if driver.name in {"optimize", "meci", "mecp", "tci", "mep", "ts", "irc", "neb"}:
        if "maxit" in options:
            if not _is_positive_integer(options["maxit"]):
                raise OQPInputError("maxit must be a positive integer")
        coordsys = str(options.get("coordsys", "tric")).strip().lower()
        if "coordsys" in options and coordsys not in {
            "auto", "tric", "dlc", "ric", "internal", "cart", "cartesian",
        }:
            raise OQPInputError(
                "coordsys must be auto, tric, dlc, ric, internal, cart, or cartesian"
            )
        if driver.name == "optimize" and backend == "geometric":
            try:
                trust = float(options.get("trust", 0.1))
                tmax = float(options.get("tmax", 0.3))
            except (TypeError, ValueError) as exc:
                raise OQPInputError(
                    "geomeTRIC trust and tmax must be positive numbers"
                ) from exc
            if (not math.isfinite(trust) or not math.isfinite(tmax)
                    or trust <= 0 or tmax <= 0 or trust > tmax):
                raise OQPInputError(
                    "geomeTRIC trust and tmax require 0 < trust <= tmax"
                )
        elif {"trust", "trust_max"}.intersection(options):
            try:
                trust = float(options.get("trust", 0.2))
                trust_max = float(options.get("trust_max", 0.5))
            except (TypeError, ValueError) as exc:
                raise OQPInputError("trust and trust_max must be positive numbers") from exc
            if (not math.isfinite(trust) or not math.isfinite(trust_max)
                    or trust <= 0 or trust_max <= 0 or trust > trust_max):
                raise OQPInputError("trust and trust_max require 0 < trust <= trust_max")
        if "auto_recovery" in options and not isinstance(
                options["auto_recovery"], bool):
            raise OQPInputError("auto_recovery must be true or false")
        if "recovery_maxit" in options and not _is_positive_integer(
                options["recovery_maxit"]):
            raise OQPInputError("recovery_maxit must be a positive integer")
        if "recovery_trust" in options:
            try:
                recovery_trust = float(options["recovery_trust"])
                if not math.isfinite(recovery_trust) or recovery_trust <= 0:
                    raise ValueError
            except (TypeError, ValueError) as exc:
                raise OQPInputError(
                    "recovery_trust must be a positive number"
                ) from exc
        if "freeze" in options:
            if driver.name != "optimize":
                raise OQPInputError("freeze is currently supported only by opt(...)")
            _validate_freeze_spec(options["freeze"])
        if driver.name == "mep":
            if "points" in options and "maxit" in options:
                raise OQPInputError("mep points and maxit specify the same native path limit; use one")
            if "step" in options and "mep_step" in options:
                raise OQPInputError("mep step and mep_step specify the same native step; use one")
            if "points" in options:
                if not _is_positive_integer(options["points"]):
                    raise OQPInputError("mep points must be a positive integer")
            for key in ("step", "mep_step"):
                if key in options:
                    try:
                        value = float(options[key])
                        if not math.isfinite(value) or value <= 0:
                            raise ValueError
                    except (TypeError, ValueError) as exc:
                        raise OQPInputError("mep %s must be a positive number" % key) from exc
            if "gtol" in options:
                try:
                    value = float(options["gtol"])
                    if not math.isfinite(value) or value <= 0:
                        raise ValueError
                except (TypeError, ValueError) as exc:
                    raise OQPInputError("mep gtol must be a positive number") from exc
        if driver.name == "ts":
            if "follow" in options:
                if not _is_non_negative_integer(options["follow"]):
                    raise OQPInputError("ts follow must be a non-negative mode index")
            if "hessian" in options:
                hessian = str(options["hessian"]).strip().lower()
                if hessian not in {"model", "numerical", "analytical"}:
                    raise OQPInputError(
                        "ts hessian must be model, numerical, or analytical"
                    )
        if driver.name == "irc":
            if "step" in options and "irc_step" in options:
                raise OQPInputError("irc step and irc_step specify the same native step; use one")
            for key in ("step", "irc_step"):
                if key in options:
                    try:
                        value = float(options[key])
                        if not math.isfinite(value) or value <= 0:
                            raise ValueError
                    except (TypeError, ValueError) as exc:
                        raise OQPInputError("irc %s must be a positive number" % key) from exc
            if "direction" in options and str(options["direction"]).strip().lower() not in {
                "forward", "backward", "reverse",
            }:
                raise OQPInputError("irc direction must be forward or backward")
            if "hessian" in options:
                hessian = str(options["hessian"]).strip().lower()
                if hessian not in {"numerical", "analytical"}:
                    raise OQPInputError(
                        "irc hessian must be numerical or analytical"
                    )
            if "gtol" in options:
                try:
                    value = float(options["gtol"])
                    if not math.isfinite(value) or value <= 0:
                        raise ValueError
                except (TypeError, ValueError) as exc:
                    raise OQPInputError("irc gtol must be a positive number") from exc
    if driver.name == "neb" and "product" not in driver.kwargs:
        raise OQPInputError('neb requires product="product.xyz"')
    if driver.name == "neb":
        if "images" in options and "nimage" in options:
            raise OQPInputError("neb images and nimage specify the same image count; use one")
        for key in ("images", "nimage"):
            if key in options:
                if not _is_integer(options[key]) or options[key] < 3:
                    raise OQPInputError(
                        "neb %s must be an integer of at least 3" % key
                    )
        if "dt" in options and "neb_dt" in options:
            raise OQPInputError("neb dt and neb_dt specify the same FIRE time step; use one")
        for key in ("climb", "align", "opt_ends"):
            if key in options and not isinstance(options[key], bool):
                raise OQPInputError("neb %s must be true or false" % key)
        for key in ("spring", "fmax", "frms", "climb_fmax", "dt", "neb_dt", "maxmove", "end_fmax"):
            if key in options:
                try:
                    value = float(options[key])
                    if not math.isfinite(value) or value <= 0:
                        raise ValueError
                except (TypeError, ValueError) as exc:
                    raise OQPInputError("neb %s must be a positive number" % key) from exc
        if (options.get("climb", True)
                and float(options.get("climb_fmax", 0.05))
                < float(options.get("fmax", 2.0e-3))):
            raise OQPInputError(
                "neb climb_fmax must be greater than or equal to fmax when climb=true"
            )
    if driver.name == "soc":
        soc_options = _driver_options(driver)
        has_ns = "ns" in soc_options
        has_nt = "nt" in soc_options
        if has_ns != has_nt:
            raise OQPInputError("soc requires ns and nt together")
        if has_ns and "nstate" in spec.model_options:
            raise OQPInputError(
                "Do not combine route nstate with soc(ns=...,nt=...); choose equal or explicit counts"
            )
        for key in ("ns", "nt"):
            if key in soc_options:
                if not _is_positive_integer(soc_options[key]):
                    raise OQPInputError("soc %s must be a positive integer" % key)
    if driver.name in {"hess", "thermo"} and "type" in driver.kwargs:
        hess_type = str(driver.kwargs["type"]).strip().lower()
        if hess_type not in {"numerical", "analytical"}:
            raise OQPInputError("%s type must be numerical or analytical" % driver.name)
    if driver.name in {"nac", "bp", "nacme"} and "type" in driver.kwargs:
        nac_type = str(driver.kwargs["type"]).strip().lower()
        if nac_type != "numerical":
            raise OQPInputError(
                "%s currently supports type=numerical only; analytical NAC is unavailable"
                % driver.name
            )
    qmmm_section = next(
        (call for call in spec.modifiers if call.name == "qmmm"), None
    )
    input_section = next(
        (call for call in spec.modifiers if call.name == "input"), None
    )
    qmmm_flag_sources: List[Tuple[str, Any]] = []
    if "qmmm_flag" in spec.options:
        qmmm_flag_sources.append(("top-level qmmm", spec.options["qmmm_flag"]))
    if input_section is not None and "qmmm_flag" in input_section.kwargs:
        qmmm_flag_sources.append(("input(qmmm_flag=...)", input_section.kwargs["qmmm_flag"]))
    if qmmm_section is not None:
        qmmm_flag_sources.append(("qmmm(...)", True))
    if len(qmmm_flag_sources) > 1:
        sources = ", ".join(source for source, _ in qmmm_flag_sources)
        raise OQPInputError(
            "QM/MM activation is specified more than once (%s). "
            "Use qmmm(...) alone when supplying QM/MM settings." % sources
        )

    def requested(value: Any) -> bool:
        if isinstance(value, str):
            return value.strip().lower() not in {"", "0", "false", "no", "off"}
        return bool(value)

    has_qmmm = qmmm_section is not None or any(
        requested(value) for _, value in qmmm_flag_sources
    )
    if driver.name == "md":
        if not has_qmmm:
            raise OQPInputError("md(...) is the QM/MM molecular-dynamics driver and requires qmmm(...)")
    if has_qmmm and driver.name not in ACTIVE_QMMM_DRIVERS:
        raise OQPInputError(
            "The active QM/MM backend supports energy, md, and namd. "
            "%s is not connected and would otherwise run without the requested QM/MM forces."
            % driver.name
        )
    if states and set(_driver_options(driver)).intersection({"active", "init_state"}):
        key = sorted(set(_driver_options(driver)).intersection({"active", "init_state"}))[0]
        raise OQPInputError(
            "%s is an internal state selector and cannot be combined with a physical driver state"
            % key
        )
    if driver.name in {"nac", "bp", "nacme"}:
        if model not in {"mrsf", "mrsf-hf", "mrsf-dftb"}:
            raise OQPInputError("%s currently requires an MRSF route" % driver.name)
        if len({state.multiplicity for state in states}) != 1:
            raise OQPInputError(
                "%s requires two states in the same spin manifold; use mecp or soc for mixed spin"
                % driver.name
            )
        if driver.name == "bp" and model == "mrsf-dftb":
            raise OQPInputError("bp is not available for MRSF-TDDFTB")
    if driver.name == "nacme":
        has_previous = "geom2" in spec.options or any(
            call.name == "guess" and bool(call.kwargs.get("file2"))
            for call in spec.modifiers
        )
        if not has_previous:
            raise OQPInputError('nacme requires geom2="previous.xyz" or guess(file2="previous.json")')

    for state in states:
        if model in {"sf", "sf-hf", "sf-dftb"} and state.label is not None:
            raise OQPInputError(
                "SF-TDDFT states are not safely labelable before the run; use root=N"
            )
        if model in {"tddftb", "tda-dftb"} and state.label and \
                state.label.startswith("T"):
            raise OQPInputError(
                "Conventional TD-DFTB currently supports singlet targets only; "
                "use sf-tddftb or mrsf-tddftb for triplet-state calculations"
            )
        if model in {
            "dft", "rks", "uks", "roks", "hf", "rhf", "uhf", "rohf", "mp2",
            "dftb", "dftb0",
        } and state.label != "S0":
            raise OQPInputError("%s only supports the ground-state label S0" % model)
        if model in {
            "tddft", "tda", "tda-hf", "tdhf", "tddftb", "tda-dftb"
        } and state.label:
            if state.label.startswith("Q"):
                raise OQPInputError("Conventional TD calculations accept singlet or triplet labels")
        if model == "mrsf-dftb" and state.label and state.label.startswith("Q"):
            raise OQPInputError("MRSF-TDDFTB supports singlet or triplet labels, not quintet")

    if driver.name == "meci" and len({state.multiplicity for state in states}) != 1:
        raise OQPInputError("MECI states must have the same multiplicity; use mecp for different spins")
    if driver.name == "mecp" and len({state.multiplicity for state in states}) != 2:
        raise OQPInputError("MECP states must have different multiplicities")
    if driver.name == "tci" and len({state.multiplicity for state in states}) != 1:
        raise OQPInputError("TCI states must have the same multiplicity")
    if model in {
        "tddft", "tda", "tda-hf", "tdhf", "tddftb", "tda-dftb"
    } and driver.name in {
        "meci", "mecp", "tci"
    } and any(_internal_root(model, state) == 0 for state in states):
        raise OQPInputError(
            "%s cannot use conventional-TD ground root S0; use positive excited states or MRSF"
            % driver.name
        )
    if driver.name in {"meci", "tci", "nac", "bp", "nacme"}:
        roots_with_spin = [
            (state.multiplicity, _internal_root(model, state)) for state in states
        ]
        if len(set(roots_with_spin)) != len(roots_with_spin):
            raise OQPInputError("%s requires distinct states" % driver.name)

    if "nstate" in spec.model_options:
        nstate = spec.model_options["nstate"]
        if not _is_positive_integer(nstate):
            raise OQPInputError("nstate must be a positive integer")
        roots = [_internal_root(spec.model, state) for state in states]
        if roots and max(roots) > nstate:
            raise OQPInputError(
                "Requested state needs response root %d, but nstate=%d" % (max(roots), nstate)
            )

    # Generic section calls preserve current keywords, but semantic/internal
    # keys cannot be used to contradict the physical route/state spelling.
    reserved = {
        "input": {"method", "runtype", "system", "system2", "basis", "functional", "charge"},
        "tdhf": {"type", "multiplicity", "nstate", "target"},
        "properties": {"grad"},
        "optimize": {"istate", "jstate", "kstate", "imult", "jmult"},
        "hess": {"state"},
        "nac": {"states"},
        "md": {"active", "init_state"},
        "qmmm": {"istate"},
        "scf": {"type", "multiplicity"},
    }
    if model in {
        "dftb", "dftb0", "tddftb", "tda-dftb", "sf-dftb", "mrsf-dftb"
    }:
        reserved["dftb"] = {"type", "target_multiplicity", "reference_multiplicity"}
    if model == "mp2":
        reserved["mp2"] = set(spec.model_options)
    modifier_names: set[str] = set()
    for call in spec.modifiers:
        if call.name in modifier_names:
            raise OQPInputError("Duplicate modifier/section call: %s" % call.name)
        modifier_names.add(call.name)
    modifier_by_name = {call.name: call for call in spec.modifiers}

    tdhf_call = modifier_by_name.get("tdhf")
    if tdhf_call is not None:
        for key in ("nstate_s", "nstate_t"):
            if key in tdhf_call.kwargs and not _is_positive_integer(
                tdhf_call.kwargs[key]
            ):
                raise OQPInputError(
                    "tdhf %s must be a positive integer" % key
                )
    if driver.name == "soc" and tdhf_call is not None:
        split_keys = {"nstate_s", "nstate_t"}.intersection(tdhf_call.kwargs)
        if len(split_keys) == 1:
            raise OQPInputError(
                "SOC tdhf(...) requires nstate_s and nstate_t together"
            )
        public_split = {"ns", "nt"}.intersection(_driver_options(driver))
        if split_keys and public_split:
            raise OQPInputError(
                "SOC state counts are specified twice; use soc(ns=...,nt=...) "
                "or tdhf(nstate_s=...,nstate_t=...), not both"
            )
        if split_keys and "nstate" in spec.model_options:
            raise OQPInputError(
                "Do not combine route nstate with tdhf.nstate_s/nstate_t for SOC; "
                "choose equal or explicit per-manifold counts"
            )

    input_call = modifier_by_name.get("input")
    if driver.name == "soc" and "soc_2e" in _driver_options(driver) \
            and input_call is not None and "soc_2e" in input_call.kwargs:
        raise OQPInputError(
            "SOC two-electron treatment is specified twice; use soc(soc_2e=...) "
            "or input(soc_2e=...), not both"
        )

    nmr_call = modifier_by_name.get("nmr")
    properties_call = modifier_by_name.get("properties")
    if nmr_call is not None:
        property_conflicts = {"scf_prop", "nmr_gauge"}
        if properties_call is not None:
            conflict = property_conflicts.intersection(properties_call.kwargs)
            if conflict:
                raise OQPInputError(
                    "nmr and properties(%s=...) specify the same NMR request; use one form"
                    % sorted(conflict)[0]
                )
        if driver.name in {"prop", "data"}:
            conflict = property_conflicts.intersection(_driver_options(driver))
            if conflict:
                raise OQPInputError(
                    "nmr and %s(%s=...) specify the same property; use one form"
                    % (driver.name, sorted(conflict)[0])
                )

    if "d4" in modifier_by_name and input_call is not None \
            and "d4" in input_call.kwargs:
        raise OQPInputError("DFT-D4 is specified twice; use d4 once")

    for call in spec.modifiers:
        legacy_keys = LEGACY_ONLY_SCHEMA_KEYS.get(call.name, frozenset()).intersection(
            call.kwargs
        )
        if call.name == "geometric" or legacy_keys:
            key_text = ".%s" % sorted(legacy_keys)[0] if legacy_keys else ""
            raise OQPInputError(
                "%s%s is available only in traditional sectioned .inp files. "
                "Concise .oqp geometry workflows use the native OpenQP optimizer "
                "automatically." % (call.name, key_text)
            )
        if call.name in SECTION_NAMES:
            unknown = set(call.kwargs) - OQP_SCHEMA_KEYS[call.name]
            if unknown:
                key = sorted(unknown)[0]
                close = difflib.get_close_matches(
                    key, OQP_SCHEMA_KEYS[call.name], n=1, cutoff=0.6
                )
                hint = " Did you mean '%s'?" % close[0] if close else ""
                raise OQPInputError(
                    "Unknown option '%s.%s'.%s" % (call.name, key, hint)
                )
        if call.name in {"ir", "raman"}:
            if call.args or call.kwargs:
                raise OQPInputError("%s does not accept options" % call.name)
            if spec.driver.name not in {"hess", "thermo"}:
                raise OQPInputError("%s requires hess(...) or thermo() as the primary driver" % call.name)
        if call.name == "d4" and (call.args or call.kwargs):
            raise OQPInputError("d4() does not accept options")
        if call.name == "d4" and "d4" in spec.options:
            raise OQPInputError("DFT-D4 is specified twice; use d4() once")
        if call.name in SECTION_NAMES and call.args:
            if call.name != "pcm":
                raise OQPInputError("Section call %s accepts keyword arguments only" % call.name)
        if (
            call.name == "qmmm"
            and driver.name == "md"
            and {"nsteps", "n_steps"}.issubset(call.kwargs)
        ):
            raise OQPInputError(
                "qmmm nsteps and n_steps specify the same QM/MM MD step count; use one"
            )
        if call.name == "pcm" and call.args and "solvent" in call.kwargs:
            raise OQPInputError("PCM solvent is specified twice; use pcm(water) or pcm(solvent=water)")
        conflict = reserved.get(call.name, set()).intersection(call.kwargs)
        if conflict:
            if call.name == "qmmm" and "istate" in conflict:
                raise OQPInputError(
                    "qmmm.istate is an obsolete numeric selector used only by the "
                    "disconnected legacy libopenmm driver. Active QM/MM energy uses "
                    "the method route; state-aware QM/MM gradients/optimizations are "
                    "not currently connected. Remove qmmm(istate=...)."
                )
            raise OQPInputError(
                "%s(%s=...) is route-owned; use physical method/state syntax instead"
                % (call.name, sorted(conflict)[0])
            )
        if call.name == "input":
            optional_input_owners = {
                "library": "library", "ispher": "ispher", "perf": "perf",
                "d4": "d4", "qmmm_flag": "qmmm_flag", "omp_threads": "omp_threads",
            }
            duplicate = [
                key for key, option in optional_input_owners.items()
                if key in call.kwargs and option in spec.options
            ]
            if duplicate:
                raise OQPInputError(
                    "input(%s=...) duplicates the top-level %s option"
                    % (duplicate[0], optional_input_owners[duplicate[0]])
                )
        overlap = set(call.kwargs).intersection(_driver_options(driver))
        if call.name == {
            "grad": "properties", "optimize": "optimize", "meci": "optimize",
            "mecp": "optimize", "tci": "optimize", "mep": "optimize",
            "ts": "optimize", "irc": "optimize", "neb": "optimize",
            "hess": "hess", "nac": "nac", "bp": "nac", "nacme": "nac",
            "md": "md", "namd": "md", "ekt": "ekt", "prop": "properties",
            "data": "properties", "thermo": "hess",
        }.get(driver.name) and overlap:
            raise OQPInputError(
                "Option '%s' is specified in both %s(...) and %s(...)"
                % (sorted(overlap)[0], driver.name, call.name)
            )

    oqp_call = modifier_by_name.get("oqp")
    if oqp_call is not None:
        allowed_oqp = OQP_DRIVER_OPTIONS.get(driver.name, set())
        ignored = set(oqp_call.kwargs) - allowed_oqp
        if ignored:
            key = sorted(ignored)[0]
            raise OQPInputError(
                "oqp.%s is not used by the native %s workflow"
                % (key, driver.name)
            )
        if "init_hessian" in oqp_call.kwargs:
            policy = str(oqp_call.kwargs["init_hessian"]).strip().lower()
            if policy not in {"model", "numerical", "analytical"}:
                raise OQPInputError(
                    "oqp.init_hessian must be model, numerical, or analytical"
                )
        if "irc_direction" in oqp_call.kwargs:
            direction = str(oqp_call.kwargs["irc_direction"]).strip().lower()
            if direction not in {"forward", "backward", "reverse"}:
                raise OQPInputError(
                    "oqp.irc_direction must be forward or backward"
                )
        if "coordsys" in oqp_call.kwargs:
            coordsys = str(oqp_call.kwargs["coordsys"]).strip().lower()
            if coordsys not in {
                "auto", "tric", "dlc", "ric", "internal", "cart", "cartesian",
            }:
                raise OQPInputError(
                    "oqp.coordsys must be auto, tric, dlc, ric, internal, cart, or cartesian"
                )
        if {"trust", "trust_max"}.intersection(oqp_call.kwargs):
            try:
                trust = float(oqp_call.kwargs.get("trust", 0.2))
                trust_max = float(oqp_call.kwargs.get("trust_max", 0.5))
            except (TypeError, ValueError) as exc:
                raise OQPInputError("oqp.trust and oqp.trust_max must be positive numbers") from exc
            if (not math.isfinite(trust) or not math.isfinite(trust_max)
                    or trust <= 0 or trust_max <= 0 or trust > trust_max):
                raise OQPInputError("oqp trust controls require 0 < trust <= trust_max")
        if "freeze" in oqp_call.kwargs:
            _validate_freeze_spec(oqp_call.kwargs["freeze"])
        if "follow" in oqp_call.kwargs:
            if not _is_non_negative_integer(oqp_call.kwargs["follow"]):
                raise OQPInputError("oqp.follow must be a non-negative mode index")
        for key in {
            "spring", "fmax", "frms", "climb_fmax", "neb_dt", "maxmove",
            "end_fmax", "irc_step", "mep_step", "path_gtol",
        }.intersection(oqp_call.kwargs):
            try:
                value = float(oqp_call.kwargs[key])
                if not math.isfinite(value) or value <= 0:
                    raise ValueError
            except (TypeError, ValueError) as exc:
                raise OQPInputError("oqp.%s must be a positive number" % key) from exc
        for key in {"climb", "align", "opt_ends"}.intersection(oqp_call.kwargs):
            if not isinstance(oqp_call.kwargs[key], bool):
                raise OQPInputError("oqp.%s must be true or false" % key)
        native_aliases = {
            "coordsys": "coordsys", "trust": "trust", "trust_max": "trust_max",
        }
        if driver.name == "optimize":
            native_aliases.update({"freeze": "freeze"})
        elif driver.name == "ts":
            native_aliases.update({"follow": "follow", "hessian": "init_hessian"})
        elif driver.name == "mep":
            native_aliases.update({"step": "mep_step", "mep_step": "mep_step"})
        elif driver.name == "irc":
            native_aliases.update({
                "direction": "irc_direction", "step": "irc_step", "irc_step": "irc_step",
            })
        elif driver.name == "neb":
            native_aliases.update({
                "spring": "spring", "climb": "climb", "fmax": "fmax",
                "frms": "frms", "climb_fmax": "climb_fmax", "dt": "neb_dt",
                "neb_dt": "neb_dt", "maxmove": "maxmove", "align": "align",
                "opt_ends": "opt_ends", "end_fmax": "end_fmax",
                "output": "neb_output",
            })
        for public_key, section_key in native_aliases.items():
            if public_key in options and section_key in oqp_call.kwargs:
                raise OQPInputError(
                    "%s and oqp.%s specify the same native control; use one form"
                    % (public_key, section_key)
                )
        if driver.name == "neb":
            effective_climb = oqp_call.kwargs.get(
                "climb", options.get("climb", True)
            )
            effective_fmax = float(oqp_call.kwargs.get(
                "fmax", options.get("fmax", 2.0e-3)
            ))
            effective_climb_fmax = float(oqp_call.kwargs.get(
                "climb_fmax", options.get("climb_fmax", 0.05)
            ))
            if effective_climb and effective_climb_fmax < effective_fmax:
                raise OQPInputError(
                    "neb climb_fmax must be greater than or equal to fmax when climb=true"
                )


def _internal_root(model: str, state: StateRef) -> int:
    if state.root is not None:
        return state.root
    if state.label is None:
        raise OQPInputError("Missing target state")
    index = state.physical_index
    if model in {"mrsf", "mrsf-hf", "mrsf-dftb", "umrsf", "umrsf-hf"}:
        return index + 1
    if model in {
        "tddft", "tda", "tda-hf", "tdhf", "tddftb", "tda-dftb"
    } and state.label[0] in {"T", "Q"}:
        return index + 1
    return index


def _default_state(spec: CalculationSpec) -> Optional[StateRef]:
    if spec.driver.name in {
        "grad", "optimize", "mep", "ts", "irc", "neb", "hess", "thermo", "ekt", "prop"
    }:
        return StateRef(label="S0")
    if spec.driver.name == "md":
        return StateRef(label="S0")
    ground_models = {
        "dft", "rks", "uks", "roks", "hf", "rhf", "uhf", "rohf", "mp2",
        "dftb", "dftb0",
    }
    if spec.driver.name == "namd" and {
        "active", "init_state"
    }.intersection(_driver_options(spec.driver)):
        # Numeric ``active`` addresses the NAMD propagation manifold (and, for
        # SOC, the spin-adiabatic surfaces), while ``init_state`` selects MCH
        # character.  Either is already a complete selector and must not be
        # replaced by the default physical state below.
        return None
    if spec.driver.name in {"namd", "data"} and spec.model not in ground_models:
        if spec.model in {"sf", "sf-hf", "sf-dftb"}:
            return StateRef(root=1)
        return StateRef(label="S1")
    return None


def _driver_options(driver: CallSpec) -> Dict[str, Any]:
    state_keys = {"state", "root"}
    if driver.name in {"meci", "mecp", "nac", "bp", "nacme"}:
        state_keys.update({"state1", "state2", "i", "j"})
    elif driver.name == "tci":
        state_keys.update({"state1", "state2", "state3", "i", "j", "k"})
    return {key: value for key, value in driver.kwargs.items() if key not in state_keys}


def _normalized_meci_options(driver: CallSpec) -> Dict[str, Any]:
    """Return MECI options using their legacy/configuration key names.

    The concise surface uses readable names such as ``algorithm`` and
    ``delta_beta``.  Existing spellings such as ``meci_search`` and
    ``pen_delta`` remain accepted, but specifying both forms is an error rather
    than a last-one-wins ambiguity.
    """

    options = _driver_options(driver)
    normalized: Dict[str, Any] = {}
    sources: Dict[str, str] = {}
    for key, value in options.items():
        target = _MECI_OPTION_ALIASES.get(key, key)
        if target in normalized:
            raise OQPInputError(
                "MECI option '%s' is specified twice through %s and %s"
                % (target, sources[target], key)
            )
        normalized[target] = value
        sources[target] = key
    return normalized


def _validate_freeze_spec(value: Any) -> str:
    """Validate the concise native frozen-distance expression."""

    text = str(value).strip()
    if not text:
        raise OQPInputError("freeze requires distance(i,j) with one-based atom indices")
    items = [item.strip() for item in text.split(";") if item.strip()]
    pattern = re.compile(r"(?:distance|r)\(\s*(\d+)\s*,\s*(\d+)\s*\)", re.I)
    if not items:
        raise OQPInputError("freeze requires distance(i,j) with one-based atom indices")
    for item in items:
        match = pattern.fullmatch(item)
        if not match or int(match.group(1)) < 1 or int(match.group(2)) < 1 \
                or match.group(1) == match.group(2):
            raise OQPInputError(
                "freeze accepts semicolon-separated distance(i,j) pairs with distinct positive indices"
            )
    return ";".join(items)


def _as_config_string(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (list, tuple)):
        if value and all(isinstance(item, (list, tuple)) for item in value):
            return ", ".join(" ".join(str(x) for x in item) for item in value)
        return ",".join(str(item) for item in value)
    return str(value)


def _resolve_path(value: Any, source_dir: Optional[Path]) -> Any:
    if source_dir is None or not isinstance(value, str) or "\n" in value:
        return value
    # Inline geometry is a whitespace-separated atom table, not a filename.
    # File-backed geometries/products are normally .xyz/.pdb; retain arbitrary
    # one-token names so the existing checker can give its usual file error.
    stripped = value.strip()
    # An empty value means "unset" (for example dftb.parameter_path resolving
    # the bundled default); Path("") would otherwise lower it to the .oqp
    # source directory itself.
    if not stripped:
        return value
    # QM/MM compact geometry references append atom selectors after the PDB
    # filename (``system.pdb 1 2 5-8``).  Resolve only the path prefix against
    # the .oqp source directory and preserve the selector suffix verbatim.
    pdb_end = stripped.lower().find(".pdb")
    if pdb_end >= 0:
        pdb_end += len(".pdb")
        selector_suffix = stripped[pdb_end:]
        if selector_suffix.strip():
            candidate = Path(stripped[:pdb_end]).expanduser()
            resolved = candidate if candidate.is_absolute() else source_dir / candidate
            return str(resolved.resolve()) + selector_suffix
    if any(char.isspace() for char in stripped) and not stripped.lower().endswith((".xyz", ".pdb")):
        return value
    candidate = Path(value).expanduser()
    if candidate.is_absolute():
        return str(candidate)
    return str((source_dir / candidate).resolve())


def _normalize_geometry(value: Any, source_dir: Optional[Path]) -> Any:
    """Resolve geometry files and normalize compact inline atom tables."""

    value = _resolve_path(value, source_dir)
    if not isinstance(value, str) or "\n" in value:
        return value
    tokens = value.split()
    if len(tokens) < 8 or len(tokens) % 4:
        return value
    try:
        for offset in range(0, len(tokens), 4):
            float(tokens[offset + 1])
            float(tokens[offset + 2])
            float(tokens[offset + 3])
    except ValueError:
        return value
    return "\n".join(
        " ".join(tokens[offset:offset + 4])
        for offset in range(0, len(tokens), 4)
    )


def _qmmm_pdb_from_geometry(value: Any, source_dir: Optional[Path]) -> Optional[str]:
    """Extract the PDB path from ``geom="system.pdb ATOM-SELECTION"``."""

    if not isinstance(value, str) or "\n" in value:
        return None
    stripped = value.strip()
    pdb_end = stripped.lower().find(".pdb")
    if pdb_end < 0:
        return None
    pdb_path = stripped[:pdb_end + len(".pdb")]
    return str(_resolve_path(pdb_path, source_dir))


def _resolve_search_path_list(value: Any, source_dir: Optional[Path]) -> Any:
    """Resolve explicit/local files while preserving package search names."""

    if source_dir is None or not isinstance(value, str):
        return value
    resolved: List[str] = []
    for raw in value.split(","):
        item = raw.strip()
        if not item:
            continue
        candidate = Path(item).expanduser()
        local = candidate if candidate.is_absolute() else (source_dir / candidate).resolve()
        explicit = item.startswith(("./", "../", "~")) or candidate.is_absolute()
        # OpenMM package identifiers such as amber14/tip3p.xml must remain
        # searchable unless a same-named local file actually exists.
        resolved.append(str(local) if explicit or local.exists() else item)
    return ",".join(resolved)


def lower_to_legacy(
    spec: CalculationSpec, *, source_dir: Optional[Path] = None
) -> Dict[str, Dict[str, str]]:
    """Lower physical semantics to a configparser-compatible legacy dict.

    Relative geometry/product paths are resolved against the source ``.oqp``
    directory rather than the process working directory.
    """

    source_dir = Path(source_dir).resolve() if source_dir is not None else None
    config: Dict[str, Dict[str, str]] = {}

    def put(section: str, key: str, value: Any) -> None:
        config.setdefault(section, {})[key] = _as_config_string(value)

    put("input", "system", _normalize_geometry(spec.options["geom"], source_dir))
    if "geom2" in spec.options:
        put("input", "system2", _normalize_geometry(spec.options["geom2"], source_dir))
    put("input", "charge", spec.options.get("charge", 0))
    if spec.basis:
        put("input", "basis", spec.basis)
    if spec.functional:
        put("input", "functional", spec.functional)
    for key in ("library", "ispher", "perf", "d4", "qmmm_flag", "omp_threads"):
        if key in spec.options:
            put("input", key, spec.options[key])

    if spec.model in {"dft", "rks", "uks", "roks", "hf", "rhf", "uhf", "rohf"}:
        method = "hf"
    elif spec.model == "mp2":
        method = "mp2"
    elif spec.model in {
        "dftb", "dftb0", "tddftb", "tda-dftb", "sf-dftb", "mrsf-dftb"
    }:
        method = "dftb"
    else:
        method = "tdhf"

    states = list(_driver_states(spec.driver))
    if not states:
        default = _default_state(spec)
        if default is not None:
            states.append(default)

    # A conventional TD S0 derivative is a ground-state DFT/HF request.  MRSF
    # S0 remains response root 1 and therefore stays on the TDHF driver.
    if (
        spec.model in {"tddft", "tda", "tdhf"}
        and states
        and all(state.label == "S0" for state in states)
        and spec.driver.name in {"grad", "optimize", "mep", "ts", "irc", "neb", "hess"}
    ):
        method = "hf"
    put("input", "method", method)
    # Branching-plane analysis is implemented as a NAC run with nac.bp=True;
    # there is no separate Python run-function dispatch for legacy runtype=bp.
    runtime_name = {
        "bp": "nac",
        # Hessian evaluation already performs and prints thermochemistry at
        # each requested temperature; legacy runtype=thermo is not accepted by
        # the checker, so the semantic alias uses the supported Hessian path.
        "thermo": "hess",
    }.get(spec.driver.name, spec.driver.name)
    put("input", "runtype", runtime_name)

    if spec.model in {
        "mrsf", "mrsf-hf", "mrsf-dftb", "sf", "sf-hf", "sf-dftb"
    }:
        put("scf", "type", "rohf")
        put("scf", "multiplicity", 3)
    elif spec.model in {"umrsf", "umrsf-hf"}:
        put("scf", "type", "uhf")
        put("scf", "multiplicity", 3)
    else:
        ref_mult = spec.options.get("mult", 1)
        if spec.model == "mp2" and "reference" in spec.model_options:
            ref_type = str(spec.model_options["reference"]).strip().lower()
        else:
            ref_type = {
                "rhf": "rhf", "uhf": "uhf", "rohf": "rohf",
                "rks": "rhf", "uks": "uhf", "roks": "rohf",
            }.get(spec.model, "rhf" if ref_mult == 1 else "uhf")
        put("scf", "type", ref_type)
        put("scf", "multiplicity", ref_mult)

    if spec.model == "mp2":
        for key, value in spec.model_options.items():
            if key != "reference":
                put("mp2", key, value)

    if spec.model in RESPONSE_MODELS:
        td_type = {
            "mrsf": "mrsf", "mrsf-hf": "mrsf", "mrsf-dftb": "mrsf",
            "umrsf": "umrsf", "umrsf-hf": "umrsf",
            "sf": "sf", "sf-hf": "sf", "sf-dftb": "sf",
            "tddft": "rpa", "tdhf": "rpa",
            "tda": "tda", "tda-hf": "tda", "tda-dftb": "tda",
            "tddftb": "tda",
        }[spec.model]
        put("tdhf", "type", td_type)
        multiplicities = [state.multiplicity for state in states if state.multiplicity]
        target_mult = multiplicities[0] if multiplicities else 1
        put("tdhf", "multiplicity", target_mult)
        roots = [_internal_root(spec.model, state) for state in states]
        explicit_nstate = spec.model_options.get("nstate")
        nstate = explicit_nstate if explicit_nstate is not None else max(roots or [1])
        put("tdhf", "nstate", nstate)
        for key, value in spec.model_options.items():
            if key != "nstate":
                put("tdhf", key, value)

    if spec.model in {
        "dftb", "dftb0", "tddftb", "tda-dftb", "sf-dftb", "mrsf-dftb"
    }:
        put("dftb", "type", {
            "dftb": "ground", "dftb0": "ground_noscc", "tddftb": "tddftb",
            "tda-dftb": "tddftb", "sf-dftb": "sf", "mrsf-dftb": "mrsf",
        }[spec.model])
        if spec.model in {"dftb", "dftb0", "tddftb", "tda-dftb"}:
            reference_mult = spec.options.get("mult", 1)
            if "mult" in spec.options or reference_mult != 1:
                put("dftb", "reference_multiplicity", reference_mult)
        if spec.model in {"tddftb", "tda-dftb", "sf-dftb", "mrsf-dftb"}:
            put("dftb", "target_multiplicity", target_mult)

    # Apply advanced section calls.  Semantic/internal state keys were rejected
    # by _validate_semantics, so these cannot undo the safe route mapping.
    for call in spec.modifiers:
        if call.name == "d4":
            put("input", "d4", True)
            continue
        if call.name in {"ir", "raman"}:
            # Hessian workflows already compute and log both intensities.  The
            # modifier is retained in canonical text to document user intent.
            continue
        if call.name == "nmr":
            props = config.setdefault("properties", {})
            current = [item.strip() for item in props.get("scf_prop", "").split(",") if item.strip()]
            if "nmr" not in current:
                current.append("nmr")
            props["scf_prop"] = ",".join(current)
            gauge = call.kwargs.get("gauge", call.kwargs.get("nmr_gauge"))
            # The new concise surface defaults to gauge-origin-independent
            # shielding. Legacy .inp retains its established CGO default.
            props["nmr_gauge"] = _as_config_string(gauge if gauge is not None else "giao")
            unknown = set(call.kwargs) - {"gauge", "nmr_gauge"}
            if unknown or call.args:
                raise OQPInputError("nmr accepts only gauge=...")
            continue
        if call.name == "pcm":
            put("pcm", "enabled", True)
            if call.args:
                if len(call.args) != 1:
                    raise OQPInputError("pcm accepts at most one positional solvent")
                put("pcm", "solvent", call.args[0])
            for key, value in call.kwargs.items():
                put("pcm", key, value)
            continue
        for key, value in call.kwargs.items():
            target_key = key
            if (call.name, key) == ("input", "system2"):
                value = _normalize_geometry(value, source_dir)
            elif (call.name, key) in {
                ("neb", "product"), ("guess", "file"), ("guess", "file2"),
                ("dftb", "parameter_path"), ("dftb", "library_path"),
                ("geometric", "constraints_file"), ("oqp", "neb_output"),
            }:
                value = _resolve_path(value, source_dir)
            elif call.name == "qmmm" and key in {
                "pdb_file", "qm_atoms_xyz", "trajectory_file", "log_file", "energy_file"
            }:
                value = _resolve_path(value, source_dir)
            elif call.name == "qmmm" and key in {"forcefield", "forcefield_files"}:
                value = _resolve_search_path_list(value, source_dir)
            if call.name == "qmmm" and spec.driver.name == "md" and key == "nsteps":
                target_key = "n_steps"
            put(call.name, target_key, value)
        if call.name == "qmmm":
            put("input", "qmmm_flag", True)
            if "pdb_file" not in call.kwargs:
                inferred_pdb = _qmmm_pdb_from_geometry(
                    spec.options.get("geom"), source_dir
                )
                if inferred_pdb is not None:
                    put("qmmm", "pdb_file", inferred_pdb)

    roots = [_internal_root(spec.model, state) for state in states]
    driver_options = _driver_options(spec.driver)
    name = spec.driver.name
    if name == "grad":
        put("properties", "grad", roots[0] if roots else 0)
        for key, value in driver_options.items():
            put("properties", key, value)
    elif name in {"optimize", "mep", "ts", "irc", "neb"}:
        put("optimize", "istate", roots[0] if roots else 0)
        optimizer_backend = str(
            driver_options.get("lib", "oqp")
        ).strip().lower()
        for key, value in driver_options.items():
            if name == "optimize" and key == "lib":
                put("optimize", "lib", optimizer_backend)
            elif (name == "optimize" and optimizer_backend == "geometric"
                  and key in _GEOMETRIC_ENGINE_OPTIONS):
                put("geometric", key, value)
            elif name == "neb" and key in {"product", "images", "nimage"}:
                target_key = "nimage" if key in {"images", "nimage"} else "product"
                if target_key == "product":
                    value = _resolve_path(value, source_dir)
                put("neb", target_key, value)
            elif name == "neb" and key in {
                "spring", "climb", "fmax", "frms", "climb_fmax", "neb_dt",
                "maxmove", "align", "opt_ends", "end_fmax",
            }:
                put("oqp", key, value)
            elif name == "neb" and key == "dt":
                put("oqp", "neb_dt", value)
            elif name == "neb" and key == "output":
                put("oqp", "neb_output", _resolve_path(value, source_dir))
            elif name == "irc" and key in {"direction", "step", "irc_step", "hessian", "gtol"}:
                if key == "hessian":
                    put("hess", "type", value)
                elif key == "gtol":
                    put("oqp", "path_gtol", value)
                elif key in {"step", "irc_step"}:
                    put("oqp", "irc_step", value)
                else:
                    direction = str(value).strip().lower()
                    put("oqp", "irc_direction", "backward" if direction == "reverse" else direction)
            elif name == "mep" and key in {"points", "step", "mep_step", "gtol"}:
                if key == "points":
                    put("optimize", "maxit", value)
                elif key in {"step", "mep_step"}:
                    put("oqp", "mep_step", value)
                elif key == "gtol":
                    put("oqp", "path_gtol", value)
            elif name == "ts" and key in {"coordsys", "trust", "trust_max", "follow", "hessian"}:
                put("oqp", "init_hessian" if key == "hessian" else key, value)
            elif name == "optimize" and key in (
                    _NATIVE_ENGINE_OPTIONS | _NATIVE_CONSTRAINT_OPTIONS):
                put("oqp", key, value)
            else:
                put("optimize", key, value)
    elif name in {"meci", "mecp", "tci"}:
        if name in {"meci", "tci"}:
            # State order is not physical for same-spin intersections, while
            # all current optimizers require ascending internal roots.  Thus
            # user-friendly meci(S2,S1) is equivalent to meci(S1,S2).
            roots = sorted(roots)
        for key, root in zip(("istate", "jstate", "kstate"), roots):
            put("optimize", key, root)
        if name == "mecp":
            put("optimize", "imult", states[0].multiplicity)
            put("optimize", "jmult", states[1].multiplicity)
        if name == "meci":
            driver_options = _normalized_meci_options(spec.driver)
            algorithm = str(driver_options.get(
                "meci_search", "baeka" if len(roots) > 2 else "auto"
            )).strip().lower()
            driver_options["meci_search"] = algorithm
            if algorithm in {"auto", "baeka"}:
                put("optimize", "states", roots)
                driver_options.setdefault("pen_sigma", 1.0)
                driver_options.setdefault("pen_alpha", 0.02)
                driver_options.setdefault("pen_delta", 0.025)
                driver_options.setdefault(
                    "pen_jump", "10,10,25,25,100,100,1000,1000,3000"
                )
                driver_options.setdefault("energy_gap", 1.0e-4)
        for key, value in driver_options.items():
            put("oqp" if key in _NATIVE_ENGINE_OPTIONS else "optimize", key, value)
    elif name in {"hess", "thermo"}:
        put("hess", "state", roots[0] if roots else 0)
        for key, value in driver_options.items():
            put("hess", key, value)
    elif name in {"nac", "bp", "nacme"}:
        put("nac", "states", [[roots[0], roots[1]]])
        if name == "bp":
            put("nac", "bp", True)
        for key, value in driver_options.items():
            put("nac", key, value)
    elif name == "md":
        put("properties", "grad", roots[0] if roots else 0)
    elif name == "namd":
        if roots:
            if driver_options.get("soc") and states[0].label:
                put("md", "init_state", states[0].label)
            else:
                put("md", "active", roots[0])
        for key, value in driver_options.items():
            put("md", key, value)
    elif name == "ekt":
        if roots:
            put("tdhf", "target", roots[0])
        for key, value in driver_options.items():
            put("ekt", key, value)
    elif name == "soc":
        for key, value in driver_options.items():
            if key == "ns":
                put("tdhf", "nstate_s", value)
            elif key == "nt":
                put("tdhf", "nstate_t", value)
            else:
                put("input", key, value)
    elif name == "prop":
        section = "properties"
        put(section, "grad", roots[0] if roots else 1)
        for key, value in driver_options.items():
            put(section, key, value)
    elif name == "data":
        if roots:
            put("properties", "grad", roots[0])
        for key, value in driver_options.items():
            put("properties", key, value)
    elif driver_options:
        raise OQPInputError("%s does not accept driver options" % name)

    # Defensive final pass: ConfigParser.load_dict requires strings.
    return {
        section: {key: _as_config_string(value) for key, value in values.items()}
        for section, values in config.items()
    }


def _render_value(value: Any) -> str:
    if isinstance(value, StateRef):
        return str(value) if value.label is not None else str(value.root)
    if isinstance(value, str):
        if _IDENT_RE.fullmatch(value) and not _STATE_RE.fullmatch(value):
            return value
        return json.dumps(value, ensure_ascii=False)
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "null"
    return str(value)


def _render_geometry_value(value: Any) -> str:
    """Render inline coordinates as a readable triple-quoted atom block."""

    if isinstance(value, str) and "\n" in value and '"""' not in value:
        return '"""\n%s\n"""' % value.strip("\n")
    return _render_value(value)


def _render_call(call: CallSpec) -> str:
    args = [_render_value(value) for value in call.args]
    args.extend("%s=%s" % (key, _render_value(value)) for key, value in call.kwargs.items())
    display_name = {"optimize": "opt"}.get(call.name, call.name)
    if not args and (
        call.name in set(PRIMARY_ALIASES.values()) or call.name in BARE_MODIFIER_CALLS
    ):
        return display_name
    return "%s(%s)" % (display_name, ",".join(args))


def render_canonical_oqp(spec: CalculationSpec) -> str:
    """Render stable, readable canonical input that reparses identically.

    The route, options, driver, and modifiers are written as separate logical
    lines. Geometry is deliberately last; inline coordinates use a
    triple-quoted block with one atom per source line. The parser also accepts
    the equivalent single-line spelling.
    """

    model_args = ",".join(
        "%s=%s" % (key, _render_value(value)) for key, value in spec.model_options.items()
    )
    display_model = {
        "tda-dftb": "tda-tddftb", "sf-dftb": "sf-tddftb",
        "mrsf-dftb": "mrsf-tddftb", "mrsf-hf": "mrsf-tdhf",
        "umrsf-hf": "umrsf-tdhf", "sf-hf": "sf-tdhf", "tda-hf": "tda-tdhf",
    }.get(spec.model, spec.model)
    route = display_model + (("(" + model_args + ")") if model_args else "")
    if spec.functional:
        route += "/" + spec.functional
        if spec.basis:
            route += "/" + spec.basis
    elif spec.basis:
        route += "/" + spec.basis
    option_order = (
        "charge", "mult", "library", "ispher", "perf", "d4", "qmmm_flag",
        "omp_threads",
    )
    option_parts = [
        "%s=%s" % (key, _render_value(spec.options[key]))
        for key in option_order
        if key in spec.options
        and not (
            key == "charge"
            and _is_integer(spec.options[key])
            and spec.options[key] == 0
        )
        and not (
            key == "mult"
            and spec.model in DEFAULT_SINGLET_MODELS
            and _is_integer(spec.options[key])
            and spec.options[key] == 1
        )
    ]
    extra = sorted(set(spec.options) - set(option_order) - {"geom", "geom2"})
    option_parts.extend("%s=%s" % (key, _render_value(spec.options[key])) for key in extra)
    driver = _normalize_driver_defaults(spec.model, spec.driver)
    normalized_modifiers = [
        _normalize_geometry_owned_section_call(
            normalized,
            spec.options.get("geom"),
        )
        for call in spec.modifiers
        if (normalized := _normalize_default_section_call(call)) is not None
    ]
    # A bare energy calculation is the language default.  Keep an explicit
    # energy state selector (for example energy(T0)), but do not make generated
    # inputs repeat an otherwise empty ``energy`` driver line.
    calls = []
    if driver.name != "energy" or _driver_states(driver) or _driver_options(driver):
        calls.append(_render_call(driver))
    calls.extend(_render_call(call) for call in normalized_modifiers)
    geometry_parts = [
        "%s=%s" % (key, _render_geometry_value(spec.options[key]))
        for key in ("geom", "geom2")
        if key in spec.options
    ]
    return "\n".join([route] + option_parts + calls + geometry_parts) + "\n"


def _natural_model(text: str) -> Optional[str]:
    checks = (
        (r"\bumrsf(?:[- ]?tddft)?\b", "umrsf"),
        (r"\bmrsf(?:[- ]?tddftb)\b", "mrsf-dftb"),
        (r"\bsf(?:[- ]?tddftb)\b|\bsftddftb\b", "sf-dftb"),
        (r"\btda(?:[- ]?(?:td)?dftb)\b|\btddftb[- ]?tda\b", "tda-dftb"),
        (r"\bdftb0\b|\bdftb[- ]?(?:no|non)scc\b", "dftb0"),
        (r"\bmrsf(?:[- ]?tddft)?\b|다중\s*스핀\s*플립", "mrsf"),
        (r"\bsf[- ]?tddft\b|스핀\s*플립", "sf"),
        (r"\btd[- ]?dftb\b", "tddftb"),
        (r"\bdftb\b", "dftb"),
        (r"\btda\b", "tda"),
        (r"\btddft\b|\btd[- ]?dft\b", "tddft"),
        (r"\btdhf\b", "tdhf"),
        (r"\bmp2\b", "mp2"),
        (r"\bdft\b|밀도\s*범함수", "dft"),
        (r"\bhf\b|hartree[- ]?fock|하트리", "hf"),
    )
    lower = text.lower()
    for pattern, model in checks:
        if re.search(pattern, lower, re.I):
            return model
    return None


def _natural_route_components(text: str, model: str) -> Tuple[str, str]:
    functional = ""
    basis = ""
    func_match = re.search(r"\bfunctional\s*=\s*([\w+.-]+)", text, re.I)
    basis_match = re.search(r"\bbasis\s*=\s*([\w+*().,-]+)", text, re.I)
    if func_match:
        functional = func_match.group(1).lower()
    if basis_match:
        basis = basis_match.group(1).lower()

    # Natural requests commonly retain the familiar method/functional/basis
    # slash triplet even when the rest is Korean or English prose.
    slash = re.search(
        r"(?:umrsf(?:[- ]?tddft)?|mrsf(?:[- ]?tddftb?)?|sf[- ]?tddft|tddft|tdhf|tda|dft|hf|mp2)"
        r"\s*/\s*([A-Za-z0-9+*()._-]+)(?:\s*/\s*([A-Za-z0-9+*()._-]+))?",
        text, re.I,
    )
    if slash:
        first, second = slash.group(1), slash.group(2)
        if second:
            functional = first.lower()
            basis = second.lower().rstrip(".")
        elif model in {"hf", "mp2", "tdhf", "dftb", "tddftb", "mrsf-dftb"}:
            basis = first.lower().rstrip(".")
        else:
            functional = first.lower()

    if not functional and model in {"dft", "tddft", "tda", "mrsf", "umrsf", "sf"}:
        known = re.search(
            r"\b(bhhlyp|pbe0|pbe|b3lyp5?|cam[- ]?b3lyp|wb97x(?:[- ]?d)?|m06[- ]?2x)\b",
            text, re.I,
        )
        if known:
            functional = known.group(1).replace(" ", "-").lower()
    if not basis:
        known_basis = re.search(
            r"\b(6-31(?:\+\+)?g(?:\([^)]*\)|\*{0,2})|sto-3g|def2-(?:svp|tzvp|qzvp)|cc-pv[dtq]z)\b",
            text, re.I,
        )
        if known_basis:
            basis = known_basis.group(1).lower()
    return functional, basis


def _find_natural_drivers(text: str) -> List[str]:
    lower = text.lower()
    patterns = (
        ("meci", r"\bmeci\b|원뿔\s*교차"),
        ("mecp", r"\bmecp\b|최소\s*에너지\s*교차점"),
        ("tci", r"\btci\b|삼중\s*교차"),
        ("neb", r"\bneb\b|nudged elastic band"),
        ("nac", r"\bnac\b|nonadiabatic coupling|비단열\s*결합"),
        ("hess", r"\bhess(?:ian)?\b|헤시안|진동수|frequency"),
        ("grad", r"\bgrad(?:ient)?(?:\b|(?=[가-힣]))|구배"),
        ("optimize", r"\b(?:geometry\s+)?optimi[sz](?:e|ation)?\b|구조\s*최적화|최적화"),
        ("mep", r"\bmep\b|minimum energy path"),
        ("ts", r"\btransition state\b|\bts\b|전이\s*상태"),
        ("irc", r"\birc\b|intrinsic reaction coordinate"),
        ("namd", r"\bnamd\b|nonadiabatic molecular dynamics"),
        ("md", r"\bmolecular dynamics\b|분자\s*동역학"),
        ("energy", r"\bsingle[- ]?point\b|\benergy\s+calculation\b|calculate\s+(?:the\s+)?energy|단일점|에너지\s*계산"),
    )
    found = [name for name, pattern in patterns if re.search(pattern, lower, re.I)]
    # "MECI optimization" and analogous conventional phrases describe one
    # crossing-search driver, not a sequence.  Other combinations remain fatal.
    if "optimize" in found and any(name in found for name in ("meci", "mecp", "tci")):
        found.remove("optimize")
    return found


def compile_natural_request(text: str) -> CalculationSpec:
    """Compile a conservative Korean/English request without an LLM.

    Ambiguity is an error.  The caller can show the generated canonical line,
    and execution uses a second canonical parse rather than trusting this pass.
    """

    if not text.strip():
        raise OQPInputError("Empty natural-language request")
    model = _natural_model(text)
    if model is None:
        raise OQPInputError("Natural request does not identify an electronic-structure method")
    functional, basis = _natural_route_components(text, model)
    if model in {"dft", "tddft", "tda", "mrsf", "umrsf", "sf"} and not functional:
        raise OQPInputError("Natural request does not identify a density functional")
    if model not in {
        "dftb", "dftb0", "tddftb", "tda-dftb", "sf-dftb", "mrsf-dftb"
    } and not basis:
        raise OQPInputError("Natural request does not identify a basis set")

    geom_match = re.search(
        r"(?:\bgeom(?:etry)?\s*=\s*)?[\"']([^\"']+\.(?:xyz|pdb))[\"']",
        text, re.I,
    )
    if not geom_match:
        geom_match = re.search(r"\b([A-Za-z0-9_./~+-]+\.(?:xyz|pdb))", text, re.I)
    if not geom_match:
        raise OQPInputError("Natural request does not identify a geometry .xyz/.pdb file")
    options: Dict[str, Any] = {"geom": geom_match.group(1)}
    charge = re.search(r"\bcharge\s*(?:=|is)?\s*([+-]?\d+)\b|전하\s*(?:=|는|은)?\s*([+-]?\d+)", text, re.I)
    if charge:
        options["charge"] = int(charge.group(1) or charge.group(2))

    nstate_match = re.search(
        r"\bnstate\s*=\s*(\d+)\b|\b(\d+)\s+(?:singlet\s+)?roots?(?![A-Za-z0-9_])|상태\s*(\d+)\s*개",
        text, re.I,
    )
    model_options: Dict[str, Any] = {}
    if nstate_match:
        model_options["nstate"] = int(next(group for group in nstate_match.groups() if group))

    drivers = _find_natural_drivers(text)
    if len(set(drivers)) > 1:
        raise OQPInputError(
            "Natural request contains multiple primary drivers (%s); use one .oqp file per workflow"
            % ", ".join(drivers)
        )
    driver_name = drivers[0] if drivers else "energy"
    state_labels = [
        StateRef(label=value.upper())
        for value in re.findall(r"(?<![A-Za-z0-9_])([STQ]\d+)(?![A-Za-z0-9_])", text, re.I)
    ]
    required = {"meci": 2, "mecp": 2, "tci": 3, "nac": 2}.get(driver_name, 1)
    if driver_name == "energy":
        if not state_labels:
            lower = text.lower()
            if re.search(r"\btriplet\b|삼중항", lower):
                state_labels = [StateRef(label="T0")]
            elif re.search(r"\bquintet\b|오중항", lower):
                state_labels = [StateRef(label="Q0")]
            elif re.search(r"\bsinglet\b|단일항", lower):
                state_labels = [StateRef(label="S0")]
        manifolds = {state.label[0] for state in state_labels if state.label}
        if len(manifolds) > 1:
            raise OQPInputError("Energy correction request names multiple spin manifolds")
        state_labels = [StateRef(label=next(iter(manifolds)) + "0")] if manifolds else []
    elif driver_name in {"soc", "thermo", "prop"}:
        state_labels = []
    elif len(state_labels) > required:
        # Repeated mentions of the same target are harmless; distinct excess
        # labels are ambiguous.
        unique: List[StateRef] = []
        for state in state_labels:
            if state not in unique:
                unique.append(state)
        state_labels = unique
    args: List[Any] = state_labels[:required]

    driver_kwargs: Dict[str, Any] = {}
    maxit = re.search(r"\bmaxit\s*=\s*(\d+)\b|최대\s*(?:반복|iteration)?\s*(\d+)\s*회", text, re.I)
    if maxit and driver_name in {"optimize", "meci", "mecp", "tci", "mep", "ts", "irc", "neb"}:
        driver_kwargs["maxit"] = int(maxit.group(1) or maxit.group(2))
    if driver_name == "neb":
        product = re.search(r"\bproduct\s*=\s*[\"']([^\"']+)[\"']", text, re.I)
        images = re.search(r"\b(?:images|nimage)\s*=\s*(\d+)\b", text, re.I)
        if product:
            driver_kwargs["product"] = product.group(1)
        if images:
            driver_kwargs["images"] = int(images.group(1))
    driver = CallSpec(driver_name, tuple(args), driver_kwargs)
    driver = _normalize_driver_defaults(model, driver)

    modifiers: List[CallSpec] = []
    conv = re.search(
        r"\bscf\b[^\n;.]*?\bconv(?:ergence)?\s*=\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)",
        text, re.I,
    )
    if conv:
        modifiers.append(CallSpec("scf", kwargs={"conv": _parse_value(conv.group(1))}))
    pcm = re.search(
        r"\bpcm(?:\s*\(\s*([\w-]+)\s*\)|\s+solvent\s*=\s*([\w-]+))?"
        r"|용매\s*([\w-]+)",
        text, re.I,
    )
    if pcm:
        solvent = next((group for group in pcm.groups() if group), None)
        modifiers.append(CallSpec("pcm", args=((solvent,) if solvent else ())))
    nmr = re.search(r"\bnmr\b|핵자기\s*공명|차폐", text, re.I)
    if nmr:
        gauge = re.search(r"\b(giao|cgo)\b", text, re.I)
        modifiers.append(CallSpec("nmr", kwargs={"gauge": gauge.group(1).lower()} if gauge else {}))

    spec = CalculationSpec(
        model=model,
        functional=functional,
        basis=basis,
        model_options=model_options,
        options=options,
        driver=driver,
        modifiers=tuple(modifiers),
        source_text=text,
    )
    _validate_semantics(spec)
    return spec


def resolve_oqp_text(
    text: str,
    *,
    source_path: Optional[Path] = None,
    write_resolved: bool = False,
) -> OQPResolution:
    """Resolve canonical or prose text and optionally write ``.resolved.oqp``.

    Natural input is compiled, rendered, and then *reparsed* as canonical input
    before it is eligible for execution.  Canonical-looking syntax errors never
    fall back to prose interpretation.
    """

    source_path = Path(source_path).resolve() if source_path is not None else None
    was_natural = not looks_canonical(text)
    if was_natural:
        first_spec = compile_natural_request(text)
        canonical = render_canonical_oqp(first_spec)
        spec = parse_canonical_oqp(canonical)
    else:
        spec = parse_canonical_oqp(text)
        canonical = render_canonical_oqp(spec)

    resolved_path: Optional[Path] = None
    if was_natural and write_resolved:
        if source_path is None:
            raise OQPInputError("source_path is required to write a resolved .oqp file")
        resolved_path = source_path.with_name(source_path.stem + ".resolved.oqp")
        resolved_path.write_text(canonical, encoding="utf-8")
        # The on-disk artifact, not an in-memory intermediate, is what executes.
        spec = parse_canonical_oqp(resolved_path.read_text(encoding="utf-8"))
        canonical = render_canonical_oqp(spec)

    source_dir = source_path.parent if source_path is not None else None
    legacy = lower_to_legacy(spec, source_dir=source_dir)
    return OQPResolution(
        spec=spec,
        canonical_text=canonical,
        legacy_config=legacy,
        was_natural=was_natural,
        source_path=source_path,
        resolved_path=resolved_path,
    )


def resolve_oqp_file(path: Path, *, write_resolved: bool = True) -> OQPResolution:
    """Read and resolve a ``.oqp`` file without changing legacy ``.inp`` files."""

    path = Path(path).resolve()
    if path.suffix.lower() != ".oqp":
        raise OQPInputError("Semantic input files must use the .oqp extension")
    return resolve_oqp_text(
        path.read_text(encoding="utf-8"), source_path=path, write_resolved=write_resolved
    )


__all__ = [
    "CalculationSpec",
    "CallSpec",
    "DRIVER_OPTIONS",
    "GENERIC_SCHEMA_KEYS",
    "INTENTIONALLY_FORBIDDEN_SCHEMA_KEYS",
    "LEGACY_ONLY_SCHEMA_KEYS",
    "OQP_SCHEMA_KEYS",
    "OQPInputError",
    "OQPResolution",
    "ROUTE_DRIVER_SCHEMA_KEYS",
    "SCHEMA_KEY_OWNERS",
    "StateRef",
    "compile_natural_request",
    "looks_canonical",
    "lower_to_legacy",
    "parse_canonical_oqp",
    "render_canonical_oqp",
    "resolve_oqp_file",
    "resolve_oqp_text",
]
