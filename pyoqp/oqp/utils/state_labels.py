"""User-facing electronic-state labels for OpenQP logs.

The native MRSF drivers deliberately use a high-spin reference and one-based
response roots.  Those are useful implementation details, but they are not the
states a user asked OpenQP to calculate.  This module keeps that distinction in
one place so input summaries and progress logs can consistently say ``S0`` or
``S1`` while still identifying the internal ROHF/UHF reference when relevant.
"""


_SPIN_NAMES = {
    1: ("singlet", "S"),
    2: ("doublet", "D"),
    3: ("triplet", "T"),
    4: ("quartet", "Q"),
    5: ("quintet", "Q"),
}

_RUN_NAMES = {
    "energy": "Single-point energy",
    "grad": "Gradient",
    "optimize": "Geometry optimization",
    "meci": "Minimum-energy conical intersection",
    "mecp": "Minimum-energy crossing point",
    "tci": "Three-state conical intersection",
    "mep": "Minimum-energy path",
    "ts": "Transition-state search",
    "irc": "Intrinsic reaction coordinate",
    "neb": "Nudged elastic band",
    "hess": "Hessian / frequencies",
    "nac": "Nonadiabatic coupling vector",
    "nacme": "Nonadiabatic coupling matrix element",
    "soc": "Spin-orbit coupling",
    "ekt": "Extended Koopmans calculation",
    "namd": "Nonadiabatic molecular dynamics",
    "md": "Molecular dynamics",
    "prop": "Properties",
    "data": "Data export",
    "thermo": "Thermochemistry",
}


def _section(config, name):
    value = config.get(name, {}) if isinstance(config, dict) else {}
    return value if isinstance(value, dict) else {}


def _int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _int_list(value):
    if isinstance(value, (list, tuple)):
        return [_int(item) for item in value]
    if value in (None, ""):
        return []
    return [_int(item) for item in str(value).replace(",", " ").split()]


def response_type(config):
    """Return the normalized response method used by the calculation."""
    return str(_section(config, "tdhf").get("type", "")).strip().lower()


def is_mrsf(config):
    """True for MRSF/UMRSF response calculations, including TD-DFTB."""
    return response_type(config) in {"mrsf", "umrsf"}


def spin_name(multiplicity):
    """Return a readable spin name, retaining unusual multiplicities."""
    mult = _int(multiplicity, 1)
    return _SPIN_NAMES.get(mult, ("multiplicity %d" % mult, "M%d-" % mult))[0]


def public_method_name(config):
    """Return the method name a user recognizes, not its legacy dispatch key."""
    inp = _section(config, "input")
    method = str(inp.get("method", "hf")).strip().lower()
    td_type = response_type(config)

    if method == "dftb":
        dftb_type = str(_section(config, "dftb").get("type", td_type)).lower()
        if td_type == "mrsf" or dftb_type in {"mrsf", "mrsftddftb", "mrsf-tddftb"}:
            return "MRSF-TDDFTB"
        if td_type == "sf" or dftb_type in {"sf", "sftddftb", "sf-tddftb"}:
            return "SF-TDDFTB"
        return "DFTB"
    if method == "tdhf":
        dft = bool(str(inp.get("functional", "")).strip())
        return {
            "mrsf": "MRSF-TDDFT" if dft else "MRSF-TDHF",
            "umrsf": "UMRSF-TDDFT" if dft else "UMRSF-TDHF",
            "sf": "SF-TDDFT" if dft else "SF-TDHF",
            "tda": "TDA-TDDFT" if dft else "TDA-TDHF",
            "rpa": "TDDFT" if dft else "TDHF",
        }.get(td_type, "TDDFT" if dft else "TDHF")
    if method == "mp2":
        variant = str(_section(config, "mp2").get("variant", "mp2")).upper()
        return variant
    functional = str(inp.get("functional", "")).strip()
    return "DFT" if functional else "HF"


def public_state_label(config, root, multiplicity=None, *, show_reference=True):
    """Translate an engine state/root number to a physical state label.

    MRSF labels are zero-based inside each spin manifold: response root 1 is
    ``S0``, ``T0``, or ``Q0`` according to the requested target multiplicity.
    Root zero in an MRSF energy array is the high-spin SCF working reference,
    not a physical target state.
    Other methods retain the legacy numeric state label to avoid asserting a
    spin assignment that the input may not contain.
    """
    root = _int(root)
    if not is_mrsf(config):
        return "state %d" % root

    scf = _section(config, "scf")
    td = _section(config, "tdhf")
    mult = _int(td.get("multiplicity", 1) if multiplicity is None else multiplicity, 1)
    if root == 0 and show_reference:
        ref_type = str(scf.get("type", "rohf")).upper()
        ref_mult = _int(scf.get("multiplicity", 3), 3)
        return "%s %s reference (internal)" % (spin_name(ref_mult), ref_type)

    prefix = _SPIN_NAMES.get(mult, ("multiplicity %d" % mult, "M%d-" % mult))[1]
    public_index = root - 1
    return "%s%d" % (prefix, public_index)


def reference_description(config):
    """Describe the working SCF reference, explicitly marking it internal."""
    scf = _section(config, "scf")
    ref_type = str(scf.get("type", "rhf")).upper()
    ref_mult = _int(scf.get("multiplicity", 1), 1)
    return "%s %s (internal working reference)" % (spin_name(ref_mult), ref_type)


def _state_range(config):
    td = _section(config, "tdhf")
    mult = _int(td.get("multiplicity", 1), 1)
    nstate = max(1, _int(td.get("nstate", 1), 1))
    first = 0
    prefix = _SPIN_NAMES.get(mult, ("multiplicity %d" % mult, "M%d-" % mult))[1]
    last = nstate - 1
    if first == last:
        return "%s%d" % (prefix, first)
    return "%s%d-%s%d" % (prefix, first, prefix, last)


def requested_states(config):
    """Return physical target state(s) for the active workflow."""
    if not is_mrsf(config):
        return ""
    inp = _section(config, "input")
    opt = _section(config, "optimize")
    td = _section(config, "tdhf")
    runtype = str(inp.get("runtype", "energy")).strip().lower()
    td_mult = _int(td.get("multiplicity", 1), 1)

    if runtype == "soc":
        ns = max(1, _int(td.get("nstate_s", 0), 0) or _int(td.get("nstate", 1), 1))
        nt = max(1, _int(td.get("nstate_t", 0), 0) or _int(td.get("nstate", 1), 1))
        srange = "S0" if ns == 1 else "S0-S%d" % (ns - 1)
        trange = "T0" if nt == 1 else "T0-T%d" % (nt - 1)
        return "%s and %s" % (srange, trange)
    if runtype == "energy":
        return _state_range(config)
    if runtype == "ekt":
        return public_state_label(config, _section(config, "tdhf").get("target", 1))
    if runtype == "grad":
        roots = _int_list(_section(config, "properties").get("grad", []))
        return ", ".join(public_state_label(config, root) for root in roots)
    if runtype in {"prop", "data"}:
        roots = _int_list(_section(config, "properties").get("grad", []))
        return ", ".join(public_state_label(config, root) for root in roots)
    if runtype == "hess":
        return public_state_label(config, _section(config, "hess").get("state", 0))
    if runtype in {"optimize", "mep", "ts", "irc", "neb"}:
        return public_state_label(config, opt.get("istate", 1))
    if runtype == "meci":
        search = str(opt.get("meci_search", "penalty")).strip().lower()
        roots = (_int_list(opt.get("states", [])) if search == "baeka" else [])
        if not roots:
            roots = [opt.get("istate", 1), opt.get("jstate", 2)]
        return "/".join(public_state_label(config, root) for root in roots)
    if runtype == "tci":
        roots = [opt.get("istate", 1), opt.get("jstate", 2), opt.get("kstate", 3)]
        return "/".join(public_state_label(config, root) for root in roots)
    if runtype == "mecp":
        first = public_state_label(config, opt.get("istate", 1), opt.get("imult", td_mult))
        second = public_state_label(config, opt.get("jstate", 1), opt.get("jmult", 3))
        return "%s/%s" % (first, second)
    if runtype in {"nac", "nacme"}:
        pairs = _section(config, "nac").get("states", [])
        flat = []
        if isinstance(pairs, (list, tuple)):
            for item in pairs:
                flat.extend(item if isinstance(item, (list, tuple)) else [item])
        else:
            flat = _int_list(pairs)
        return ", ".join(public_state_label(config, root) for root in flat)
    if runtype == "namd":
        md = _section(config, "md")
        init_state = str(md.get("init_state", "")).strip()
        if bool(md.get("soc", False)):
            active = _int(md.get("active", 1), 1)
            active_label = "spin-adiabatic surface %d" % active
            if init_state:
                return "%s initialization; %s" % (
                    init_state.upper(), active_label
                )
            return active_label
        if init_state:
            return init_state.upper()
        return public_state_label(config, md.get("active", 1))
    return ""


def calculation_request_lines(config, source=None, resolved=None):
    """Build the concise, user-facing calculation summary shown at log start."""
    inp = _section(config, "input")
    runtype = str(inp.get("runtype", "energy")).lower()
    calculation = _RUN_NAMES.get(runtype, runtype)
    baeka = (
        runtype == "meci"
        and str(_section(config, "optimize").get("meci_search", "")).strip().lower()
        == "baeka"
    )
    if baeka:
        calculation = "BaekA multistate conical intersection"
    if runtype == "nac" and bool(_section(config, "nac").get("bp", False)):
        calculation = "Branching-plane analysis"
    lines = [
        ("Method", public_method_name(config)),
        ("Calculation", calculation),
    ]
    if baeka:
        lines.append(("Algorithm", "Baek adaptive penalty (BaekA)"))
    functional = str(inp.get("functional", "")).strip()
    basis = str(inp.get("basis", "")).strip()
    if functional and basis:
        lines.append(("Model chemistry", "%s/%s" % (functional, basis)))
    elif basis:
        lines.append(("Basis", basis))

    targets = requested_states(config)
    if targets:
        lines.append(("Physical target state(s)", targets))
    freeze = str(_section(config, "oqp").get("freeze", "")).strip()
    if freeze:
        lines.append(("Frozen distance(s)", freeze))
    if is_mrsf(config):
        target_mult = _int(_section(config, "tdhf").get("multiplicity", 1), 1)
        runtype = str(inp.get("runtype", "energy")).lower()
        multi_spin = runtype in {"soc", "mecp"} or (
            runtype == "namd" and bool(_section(config, "md").get("soc", False))
        )
        if not multi_spin:
            lines.append(("Target spin", spin_name(target_mult)))
        lines.append(("Reference", reference_description(config)))
        lines.append(("State labels", "physical labels shown; engine root numbers are internal"))
    if source:
        lines.append(("Input", str(source)))
    if resolved and str(resolved) != str(source):
        lines.append(("Resolved input", str(resolved)))
    return lines


def format_calculation_request(config, source=None, resolved=None):
    """Format :func:`calculation_request_lines` as a log-ready text block."""
    return "\n".join("   PyOQP %-28s %s" % (label + ":", value)
                     for label, value in calculation_request_lines(
                         config, source=source, resolved=resolved))
