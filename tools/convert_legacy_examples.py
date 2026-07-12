#!/usr/bin/env python3
"""Generate concise, one-line ``.oqp`` companions for legacy examples.

This is deliberately an example migration tool, not a second runtime parser.
It consumes every explicitly written legacy keyword, builds the public
``CalculationSpec``, renders it with :mod:`oqp_input`, and reparses the result.
Inputs which depend on behavior that the concise language cannot preserve are
reported as blockers and are never given a misleading ``.oqp`` companion.
"""

from __future__ import annotations

import argparse
import configparser
import hashlib
import importlib.util
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXAMPLES_ROOT = REPOSITORY_ROOT / "examples"


def _load_oqp_input():
    """Load the stdlib-only input module without initializing liboqp."""

    module_path = REPOSITORY_ROOT / "pyoqp" / "oqp" / "utils" / "oqp_input.py"
    module_name = "_openqp_example_oqp_input"
    if module_name in sys.modules:
        return sys.modules[module_name]
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:  # pragma: no cover - installation damage
        raise RuntimeError("cannot load %s" % module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


oqp_input = _load_oqp_input()
CalculationSpec = oqp_input.CalculationSpec
CallSpec = oqp_input.CallSpec
StateRef = oqp_input.StateRef


class ConversionError(ValueError):
    """A legacy example cannot be represented faithfully by concise input."""


@dataclass(frozen=True)
class GeometryAsset:
    name: str
    text: str


@dataclass
class GeometryCatalog:
    """Deduplicate real Cartesian geometries across all example inputs."""

    assets: Dict[Tuple[Tuple[str, Decimal, Decimal, Decimal], ...], GeometryAsset] = field(
        default_factory=dict
    )

    def add(self, value: str) -> GeometryAsset:
        atoms = _parse_cartesian_geometry(value)
        if atoms is None:
            raise ConversionError("geometry is not a Cartesian atom table")
        key = tuple(atoms)
        if key not in self.assets:
            canonical = "\n".join(
                "%s %s %s %s"
                % (symbol, _decimal_text(x), _decimal_text(y), _decimal_text(z))
                for symbol, x, y, z in key
            )
            digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:12]
            formula = _formula(symbol for symbol, *_ in key)
            name = "%s-%s.xyz" % (formula, digest)
            text = "%d\nOpenQP shared example geometry %s\n%s\n" % (
                len(key), formula, canonical
            )
            self.assets[key] = GeometryAsset(name=name, text=text)
        return self.assets[key]

    def ordered_assets(self) -> Tuple[GeometryAsset, ...]:
        return tuple(sorted(self.assets.values(), key=lambda asset: asset.name))


@dataclass(frozen=True)
class ConversionResult:
    source: Path
    target: Path
    text: str
    warnings: Tuple[str, ...] = ()


@dataclass(frozen=True)
class BlockedConversion:
    source: Path
    reason: str


@dataclass(frozen=True)
class ConversionRun:
    results: Tuple[ConversionResult, ...]
    blockers: Tuple[BlockedConversion, ...]
    catalog: GeometryCatalog


_ELEMENTS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe "
    "Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn "
    "Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re "
    "Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm "
    "Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
_ELEMENT_BY_NUMBER = {str(index): symbol for index, symbol in enumerate(_ELEMENTS) if index}
_ELEMENT_BY_NAME = {symbol.lower(): symbol for symbol in _ELEMENTS[1:]}
_INTEGER_RE = re.compile(r"^[+-]?\d+$")
_FLOAT_RE = re.compile(
    r"^[+-]?(?:\d+\.?(?:\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?$"
)
_MISSING = object()


def _symbol(value: str) -> Optional[str]:
    token = value.strip()
    if token in _ELEMENT_BY_NUMBER:
        return _ELEMENT_BY_NUMBER[token]
    return _ELEMENT_BY_NAME.get(token.lower())


def _parse_cartesian_geometry(
    value: str,
) -> Optional[Tuple[Tuple[str, Decimal, Decimal, Decimal], ...]]:
    lines = [line.split() for line in value.strip().splitlines() if line.strip()]
    if not lines or any(len(parts) != 4 for parts in lines):
        return None
    atoms = []
    for label, x_text, y_text, z_text in lines:
        symbol = _symbol(label)
        if symbol is None:
            return None
        try:
            coords = tuple(
                Decimal(item.replace("D", "E").replace("d", "e"))
                for item in (x_text, y_text, z_text)
            )
        except InvalidOperation:
            return None
        if any(not coordinate.is_finite() for coordinate in coords):
            return None
        atoms.append((symbol, *coords))
    return tuple(atoms)


def _decimal_text(value: Decimal) -> str:
    if value == 0:
        return "0.0"
    text = format(value.normalize(), "f")
    return text if "." in text else text + ".0"


def _formula(symbols: Iterable[str]) -> str:
    counts = Counter(symbols)
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else [])
        order.extend(sorted(set(counts) - set(order)))
    else:
        order = sorted(counts)
    return "".join(symbol + (str(counts[symbol]) if counts[symbol] != 1 else "") for symbol in order)


def _coerce(value: str) -> Any:
    text = value.strip()
    lowered = text.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if _INTEGER_RE.fullmatch(text):
        return int(text)
    if _FLOAT_RE.fullmatch(text):
        return float(text.replace("D", "E").replace("d", "e"))
    return text


def _read_legacy(path: Path) -> Dict[str, Dict[str, str]]:
    parser = configparser.ConfigParser(
        interpolation=None, allow_no_value=True, strict=False
    )
    parser.optionxform = str.lower
    with path.open(encoding="utf-8") as handle:
        parser.read_file(handle)
    return {
        section.lower(): {
            key.lower(): ("" if value is None else value.strip())
            for key, value in parser.items(section, raw=True)
        }
        for section in parser.sections()
    }


class _Consumer:
    def __init__(self, sections: Mapping[str, Mapping[str, str]]):
        self.sections = {
            section: dict(values) for section, values in sections.items()
        }

    def pop(self, section: str, key: str, default: Any = _MISSING) -> Any:
        values = self.sections.setdefault(section, {})
        if key in values:
            return values.pop(key)
        if default is _MISSING:
            raise ConversionError("missing required legacy keyword %s.%s" % (section, key))
        return default

    def has(self, section: str, key: str) -> bool:
        return key in self.sections.get(section, {})

    def take(self, section: str) -> Dict[str, str]:
        values = dict(self.sections.get(section, {}))
        self.sections[section] = {}
        return values

    def leftovers(self) -> Dict[str, Dict[str, str]]:
        return {
            section: values for section, values in self.sections.items() if values
        }


def _geometry_option(
    value: str, source: Path, examples_root: Path, catalog: GeometryCatalog
) -> str:
    if _parse_cartesian_geometry(value) is None:
        return value.strip()
    asset = catalog.add(value)
    geometry_path = examples_root / "geometries" / asset.name
    return Path(_relative_path(geometry_path, source.parent)).as_posix()


def _relative_path(target: Path, start: Path) -> str:
    import os

    return os.path.relpath(target, start)


def _state(model: str, root: int, multiplicity: int) -> StateRef:
    if root < 0:
        raise ConversionError("response root must be non-negative")
    if model in {"sf", "sf-hf", "sf-dftb"}:
        if root < 1:
            raise ConversionError("SF-TDDFT requires a positive implementation root")
        return StateRef(root=root)
    letter = {1: "S", 3: "T", 5: "Q"}.get(multiplicity)
    if letter is None:
        raise ConversionError("unsupported target multiplicity %s" % multiplicity)
    if model in {"mrsf", "mrsf-hf", "mrsf-dftb", "umrsf", "umrsf-hf"}:
        index = root - 1
    elif model in {"tddft", "tda", "tda-hf", "tdhf", "tddftb", "tda-dftb"}:
        index = root if multiplicity == 1 else root - 1
    else:
        if root != 0:
            raise ConversionError("ground-state method has nonzero state %s" % root)
        index = 0
    if index < 0:
        raise ConversionError(
            "legacy root %s is not valid for multiplicity %s" % (root, multiplicity)
        )
    return StateRef(label="%s%d" % (letter, index))


def _infer_model(
    method: str, functional: str, scf_type: str, tdhf_type: str, dftb_type: str
) -> str:
    if method == "hf":
        mapping = (
            {"rhf": "rks", "uhf": "uks", "rohf": "roks"}
            if functional
            else {"rhf": "rhf", "uhf": "uhf", "rohf": "rohf"}
        )
        if scf_type not in mapping:
            raise ConversionError("unsupported SCF type %r" % scf_type)
        return mapping[scf_type]
    if method == "mp2":
        return "mp2"
    if method == "dftb":
        mapping = {
            "auto": "dftb", "ground": "dftb", "ground_noscc": "dftb0",
            "tddftb": "tddftb", "rpa": "tddftb", "tda": "tda-dftb",
            "sf": "sf-dftb", "mrsf": "mrsf-dftb",
        }
        if dftb_type not in mapping:
            raise ConversionError("unsupported DFTB type %r" % dftb_type)
        return mapping[dftb_type]
    if method == "tdhf":
        if tdhf_type == "mrsf":
            return "mrsf" if functional else "mrsf-hf"
        if tdhf_type == "umrsf":
            return "umrsf" if functional else "umrsf-hf"
        if tdhf_type == "sf":
            return "sf" if functional else "sf-hf"
        if tdhf_type == "tda":
            return "tda" if functional else "tda-hf"
        if tdhf_type == "rpa":
            return "tddft" if functional else "tdhf"
        raise ConversionError("unsupported TDHF type %r" % tdhf_type)
    raise ConversionError("unsupported input.method %r" % method)


def _typed(values: Mapping[str, str]) -> Dict[str, Any]:
    return {key: _coerce(value) for key, value in values.items()}


def _pop_int(consumer: _Consumer, section: str, key: str, default: int) -> int:
    value = consumer.pop(section, key, str(default))
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ConversionError("%s.%s must be an integer" % (section, key)) from exc


def _consume_backend_options(
    consumer: _Consumer, runtime: str, warnings: list[str]
) -> None:
    if consumer.has("geometric", "coordsys") or consumer.sections.get("geometric"):
        raise ConversionError(
            "geomeTRIC settings/constraints have no exact concise native equivalent"
        )
    lib = consumer.pop("optimize", "lib", "oqp").strip().lower()
    if lib != "oqp":
        backend = "geomeTRIC" if lib == "geometric" else lib
        raise ConversionError(
            "%s optimizer backend is legacy-only; concise .oqp uses native OQP" % backend
        )
    optimizer = consumer.pop("optimize", "optimizer", "").strip().lower()
    if optimizer:
        if runtime == "meci" and optimizer == "bfgs":
            warnings.append("dropped optimize.optimizer=bfgs (legacy MECI default/unused)")
        else:
            raise ConversionError(
                "optimize.optimizer=%s is a legacy SciPy-backend option" % optimizer
            )
    for key in ("step_size", "step_tol", "mep_maxit"):
        if consumer.has("optimize", key):
            raise ConversionError(
                "optimize.%s is a legacy SciPy-backend option" % key
            )


def _driver_and_sections(
    consumer: _Consumer,
    runtime: str,
    model: str,
    target_mult: int,
    input_soc_2e: Optional[str],
    warnings: list[str],
) -> Tuple[CallSpec, Tuple[CallSpec, ...]]:
    driver_kwargs: Dict[str, Any] = {}
    states: list[StateRef] = []

    if runtime in {"optimize", "mep", "ts", "irc", "neb", "meci", "mecp", "tci"}:
        _consume_backend_options(consumer, runtime, warnings)

    if runtime == "grad":
        root = _pop_int(consumer, "properties", "grad", 0)
        states = [_state(model, root, target_mult)]
        driver_name = "grad"
    elif runtime in {"optimize", "mep", "ts", "irc", "neb"}:
        root = _pop_int(consumer, "optimize", "istate", 0)
        states = [_state(model, root, target_mult)]
        driver_name = runtime
        driver_kwargs.update(_typed(consumer.take("optimize")))
        if runtime == "mep":
            # OQPMEPOpt uses only istate/maxit from [optimize]; convergence is
            # controlled by oqp.path_gtol and step length by oqp.mep_step.
            for ignored in (
                "energy_shift", "energy_gap", "rmsd_grad", "rmsd_step",
                "max_grad", "max_step", "init_scf",
            ):
                if ignored in driver_kwargs:
                    driver_kwargs.pop(ignored)
                    warnings.append(
                        "dropped optimize.%s (native MEP path does not consume it)" % ignored
                    )
    elif runtime in {"meci", "tci"}:
        driver_name = runtime
        root_keys = ("istate", "jstate", "kstate") if runtime == "tci" else (
            "istate", "jstate", "kstate"
        )
        roots = []
        for key in root_keys:
            if consumer.has("optimize", key):
                roots.append(_pop_int(consumer, "optimize", key, 0))
        listed = consumer.pop("optimize", "states", "").strip()
        if listed:
            try:
                listed_roots = [int(item.strip()) for item in listed.split(",") if item.strip()]
            except ValueError as exc:
                raise ConversionError("optimize.states must be comma-separated roots") from exc
            if roots and roots != listed_roots:
                raise ConversionError("optimize.states disagrees with istate/jstate/kstate")
            roots = listed_roots
        expected = 3 if runtime == "tci" else 2
        if len(roots) < expected:
            raise ConversionError("%s needs at least %d explicit roots" % (runtime, expected))
        states = [_state(model, root, target_mult) for root in roots]
        driver_kwargs.update(_typed(consumer.take("optimize")))
        if runtime == "meci":
            public_names = {
                "meci_search": "algorithm",
                "pen_sigma": "sigma",
                "pen_alpha": "alpha",
                "pen_delta": "delta_beta",
                "pen_jump": "beta_schedule",
                "energy_gap": "gap",
            }
            driver_kwargs = {
                public_names.get(key, key): value
                for key, value in driver_kwargs.items()
            }
    elif runtime == "mecp":
        driver_name = "mecp"
        roots = [
            _pop_int(consumer, "optimize", "istate", 0),
            _pop_int(consumer, "optimize", "jstate", 0),
        ]
        mults = [
            _pop_int(consumer, "optimize", "imult", 1),
            _pop_int(consumer, "optimize", "jmult", 3),
        ]
        states = [_state(model, root, mult) for root, mult in zip(roots, mults)]
        driver_kwargs.update(_typed(consumer.take("optimize")))
    elif runtime == "hess":
        driver_name = "hess"
        root = _pop_int(consumer, "hess", "state", 0)
        states = [_state(model, root, target_mult)]
        driver_kwargs.update(_typed(consumer.take("hess")))
    elif runtime in {"nac", "nacme"}:
        driver_name = runtime
        raw_pairs = consumer.pop("nac", "states", "1 2")
        numbers = [int(item) for item in re.findall(r"\d+", raw_pairs)]
        if len(numbers) != 2:
            raise ConversionError("nac.states must contain exactly one state pair")
        states = [_state(model, root, target_mult) for root in numbers]
        driver_kwargs.update(_typed(consumer.take("nac")))
        if consumer.has("properties", "grad"):
            consumer.pop("properties", "grad")
            warnings.append("dropped properties.grad (NACME selects its pair through nac.states)")
    elif runtime == "prop":
        driver_name = "prop"
        root = _pop_int(consumer, "properties", "grad", 1)
        states = [_state(model, root, target_mult)]
        if consumer.sections.get("nac"):
            raise ConversionError("non-empty [nac] options alongside runtype=prop are ambiguous")
        consumer.take("nac")
    elif runtime == "ekt":
        driver_name = "ekt"
        driver_kwargs.update(_typed(consumer.take("ekt")))
        if consumer.has("tdhf", "target"):
            target = _pop_int(consumer, "tdhf", "target", 1)
            states = [_state(model, target, target_mult)]
    elif runtime == "soc":
        driver_name = "soc"
        if input_soc_2e is not None:
            driver_kwargs["soc_2e"] = _coerce(input_soc_2e)
        warnings.append("normalized legacy tdhf.multiplicity (SOC evaluates both spin manifolds)")
    elif runtime == "namd":
        driver_name = "namd"
        md_values = consumer.take("md")
        active = int(md_values.pop("active", md_values.pop("init_state", "1")))
        if bool(_coerce(md_values.get("soc", "false"))):
            # SOC-NAMD's active index addresses the spin-adiabatic propagation
            # manifold, not an MCH state label.  Preserve the numeric surface
            # exactly rather than inventing a Tn character selector.
            driver_kwargs["active"] = active
        else:
            states = [_state(model, active, target_mult)]
        driver_kwargs.update(_typed(md_values))
        if consumer.has("properties", "grad"):
            consumer.pop("properties", "grad")
            warnings.append("dropped properties.grad (NAMD derives gradients from active state)")
    elif runtime == "md":
        driver_name = "md"
        states = [StateRef(label="S0")]
    elif runtime in {"energy", "data"}:
        driver_name = runtime
        response_models = set(oqp_input.RESPONSE_MODELS)
        if runtime == "energy" and model in response_models and model not in {"sf", "sf-hf", "sf-dftb"}:
            states = [_state(model, 1, target_mult)]
            # energy accepts only the manifold selector, whose physical index is zero.
            states = [StateRef(label={1: "S0", 3: "T0", 5: "Q0"}[target_mult])]
        if consumer.has("properties", "grad"):
            grad = _pop_int(consumer, "properties", "grad", 0)
            if grad != 0:
                raise ConversionError("runtype=energy has non-default properties.grad=%d" % grad)
    else:
        raise ConversionError("unsupported input.runtype %r" % runtime)

    # Route-owned optimizer/Hessian/NEB/native-engine options must be merged
    # into the primary call because those section names are themselves drivers.
    oqp_values = consumer.take("oqp")
    if runtime in {"optimize", "meci", "mecp", "tci"}:
        driver_kwargs.update(_typed(oqp_values))
    elif runtime == "ts":
        aliases = {"init_hessian": "hessian"}
        driver_kwargs.update(_typed({aliases.get(k, k): v for k, v in oqp_values.items()}))
    elif runtime == "mep":
        for ignored in ("coordsys", "trust", "trust_max"):
            if ignored in oqp_values:
                oqp_values.pop(ignored)
                warnings.append("dropped oqp.%s (native MEP path does not consume it)" % ignored)
        aliases = {"mep_step": "step", "path_gtol": "gtol"}
        driver_kwargs.update(_typed({aliases.get(k, k): v for k, v in oqp_values.items()}))
    elif runtime == "irc":
        aliases = {"irc_step": "step", "irc_direction": "direction", "path_gtol": "gtol"}
        driver_kwargs.update(_typed({aliases.get(k, k): v for k, v in oqp_values.items()}))
        hess_values = consumer.take("hess")
        hess_state = int(hess_values.pop("state", "0"))
        if states and _state(model, hess_state, target_mult) != states[0]:
            raise ConversionError("hess.state disagrees with optimize.istate for IRC")
        if "type" in hess_values:
            driver_kwargs["hessian"] = _coerce(hess_values.pop("type"))
        if hess_values:
            raise ConversionError("unsupported IRC [hess] keys: %s" % ", ".join(sorted(hess_values)))
    elif runtime == "neb":
        neb_values = consumer.take("neb")
        driver_kwargs.update(_typed(neb_values))
        aliases = {"neb_dt": "dt", "neb_output": "output"}
        driver_kwargs.update(_typed({aliases.get(k, k): v for k, v in oqp_values.items()}))
    elif oqp_values:
        raise ConversionError("[oqp] options are unused by runtype=%s" % runtime)

    if runtime == "md" and consumer.sections.get("md"):
        raise ConversionError("ground-state QM/MM md options belong in [qmmm], not [md]")

    if model in {"sf", "sf-hf", "sf-dftb"} and states:
        if len(states) != 1 or states[0].root is None:
            raise ConversionError("SF workflow needs one explicit implementation root")
        driver_kwargs["root"] = states[0].root
        states = []

    driver = CallSpec(
        driver_name,
        args=tuple(states),
        kwargs=driver_kwargs,
        explicit=True,
    )

    modifiers = []
    generic_sections = (
        "mp2", "guess", "pcm", "dftb", "symmetry", "scf", "dftgrid",
        "tdhf", "properties", "qmmm", "json", "tests",
    )
    for section in generic_sections:
        values = consumer.take(section)
        if values:
            modifiers.append(CallSpec(section, kwargs=_typed(values)))
    # Empty compatibility sections are harmless, but non-empty sections must
    # never disappear silently.
    for section in ("input", "optimize", "geometric", "neb", "hess", "nac", "md", "ekt", "oqp"):
        values = consumer.take(section)
        if values:
            raise ConversionError(
                "unrepresented legacy keywords: %s" % ", ".join(
                    "%s.%s" % (section, key) for key in sorted(values)
                )
            )
    leftovers = consumer.leftovers()
    if leftovers:
        names = ["%s.%s" % (section, key) for section, values in leftovers.items() for key in values]
        raise ConversionError("unrepresented legacy keywords: %s" % ", ".join(sorted(names)))
    return driver, tuple(modifiers)


def convert_file(
    source: Path, examples_root: Path, catalog: GeometryCatalog
) -> ConversionResult:
    source = source.resolve()
    examples_root = examples_root.resolve()
    consumer = _Consumer(_read_legacy(source))
    warnings: list[str] = []

    raw_system = consumer.pop("input", "system", "")
    if not raw_system:
        raw_system = consumer.sections.get("qmmm", {}).get("pdb_file", "")
        if not raw_system:
            raise ConversionError("input.system is empty and no qmmm.pdb_file is available")
        warnings.append("synthesized geom from qmmm.pdb_file")
    geom = _geometry_option(raw_system, source, examples_root, catalog)
    raw_system2 = consumer.pop("input", "system2", "")
    geom2 = (
        _geometry_option(raw_system2, source, examples_root, catalog)
        if raw_system2 else None
    )
    charge = _coerce(consumer.pop("input", "charge", "0"))
    method = consumer.pop("input", "method", "hf").strip().lower()
    runtime = consumer.pop("input", "runtype", "energy").strip().lower()
    basis = consumer.pop("input", "basis", "6-31g*").strip().lower()
    functional = consumer.pop("input", "functional", "").strip().lower()
    scf_type = consumer.pop("scf", "type", "rhf").strip().lower()
    scf_mult = _pop_int(consumer, "scf", "multiplicity", 1)
    had_tdhf_type = consumer.has("tdhf", "type")
    tdhf_type = consumer.pop("tdhf", "type", "rpa").strip().lower()
    target_mult = _pop_int(consumer, "tdhf", "multiplicity", 1)
    dftb_type = consumer.pop("dftb", "type", "auto").strip().lower()
    model = _infer_model(method, functional, scf_type, tdhf_type, dftb_type)

    model_options: Dict[str, Any] = {}
    if model in set(oqp_input.RESPONSE_MODELS):
        model_options["nstate"] = _pop_int(consumer, "tdhf", "nstate", 1)
    elif had_tdhf_type or consumer.has("tdhf", "nstate"):
        # One historical ground-state TRAH deck carries an explicitly
        # commented-as-unused [tdhf] block.  It does not participate in the HF
        # route and must not turn the concise request into a response job.
        consumer.pop("tdhf", "nstate", "1")
        warnings.append("dropped unused [tdhf] route keys from a ground-state method")
    if model == "mp2":
        model_options["reference"] = scf_type
        for key in oqp_input.MP2_MODEL_OPTIONS - {"reference"}:
            if consumer.has("mp2", key):
                model_options[key] = _coerce(consumer.pop("mp2", key))

    if model in {"dftb", "dftb0", "tddftb", "tda-dftb", "sf-dftb", "mrsf-dftb"}:
        if basis:
            warnings.append("dropped input.basis (DFTB does not consume an AO basis)")
        if functional:
            warnings.append("dropped input.functional (DFTB route owns its parameterization)")
        basis = ""
        functional = ""

    options: Dict[str, Any] = {"geom": geom, "charge": charge}
    if geom2:
        options["geom2"] = geom2
    if model not in {"mrsf", "mrsf-hf", "mrsf-dftb", "umrsf", "umrsf-hf", "sf", "sf-hf", "sf-dftb"}:
        options["mult"] = scf_mult
    elif scf_mult != 3:
        raise ConversionError(
            "%s requires the automatic triplet reference, but legacy scf.multiplicity=%d"
            % (model, scf_mult)
        )
    for key in ("library", "ispher", "perf", "d4", "omp_threads"):
        if consumer.has("input", key):
            options[key] = _coerce(consumer.pop("input", key))
    qmmm_flag = consumer.pop("input", "qmmm_flag", "false")
    if consumer.sections.get("qmmm"):
        if str(qmmm_flag).strip().lower() not in {"true", "1", "yes", "on"}:
            warnings.append("qmmm(...) enables input.qmmm_flag automatically")
    elif str(qmmm_flag).strip().lower() in {"true", "1", "yes", "on"}:
        options["qmmm_flag"] = True
    input_soc_2e = consumer.pop("input", "soc_2e", None)

    driver, modifiers = _driver_and_sections(
        consumer, runtime, model, target_mult, input_soc_2e, warnings
    )
    spec = CalculationSpec(
        model=model,
        functional=functional,
        basis=basis,
        model_options=model_options,
        options=options,
        driver=driver,
        modifiers=modifiers,
        source_text=str(source),
    )
    text = oqp_input.render_canonical_oqp(spec)
    try:
        reparsed = oqp_input.parse_canonical_oqp(text)
    except Exception as exc:
        raise ConversionError("rendered .oqp does not reparse: %s" % exc) from exc
    if oqp_input.render_canonical_oqp(reparsed) != text:
        raise ConversionError("rendered .oqp is not canonical/idempotent")
    return ConversionResult(
        source=source,
        target=source.with_suffix(".oqp"),
        text=text,
        warnings=tuple(warnings),
    )


def convert_all(examples_root: Path = DEFAULT_EXAMPLES_ROOT) -> ConversionRun:
    examples_root = examples_root.resolve()
    catalog = GeometryCatalog()
    results = []
    blockers = []
    for source in sorted(examples_root.rglob("*.inp")):
        try:
            results.append(convert_file(source, examples_root, catalog))
        except (ConversionError, configparser.Error, ValueError) as exc:
            blockers.append(BlockedConversion(source=source.resolve(), reason=str(exc)))
    return ConversionRun(tuple(results), tuple(blockers), catalog)


def write_run(run: ConversionRun, examples_root: Path = DEFAULT_EXAMPLES_ROOT) -> None:
    geometry_root = examples_root.resolve() / "geometries"
    geometry_root.mkdir(parents=True, exist_ok=True)
    for asset in run.catalog.ordered_assets():
        (geometry_root / asset.name).write_text(asset.text, encoding="utf-8")
    for result in run.results:
        result.target.write_text(result.text, encoding="utf-8")


def _display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError:
        return str(path)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--examples", type=Path, default=DEFAULT_EXAMPLES_ROOT)
    parser.add_argument("--write", action="store_true", help="write .oqp and shared .xyz files")
    args = parser.parse_args(argv)

    run = convert_all(args.examples)
    if args.write:
        write_run(run, args.examples)
    verb = "generated" if args.write else "convertible"
    print("%s: %d" % (verb, len(run.results)))
    print("shared Cartesian geometries: %d" % len(run.catalog.assets))
    print("blocked: %d" % len(run.blockers))
    for blocker in run.blockers:
        print("BLOCKED %s: %s" % (_display_path(blocker.source), blocker.reason))
    warning_count = sum(len(result.warnings) for result in run.results)
    print("audited normalizations: %d" % warning_count)
    return 2 if run.blockers else 0


if __name__ == "__main__":
    raise SystemExit(main())
