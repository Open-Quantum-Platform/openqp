"""User-facing electronic-state labels for OpenQP logs.

The native MRSF drivers deliberately use a high-spin reference and one-based
response roots.  Those are useful implementation details, but they are not the
states a user asked OpenQP to calculate.  This module keeps that distinction in
one place so input summaries and progress logs can consistently say ``S0`` or
``S1`` while still identifying the internal ROHF/UHF reference when relevant.
"""

try:
    from oqp.utils.log_format import format_log_fields
except ModuleNotFoundError as exc:
    # Several source-only tests load this file directly, without importing the
    # OpenQP package or compiled library.  Load the pure formatting sibling in
    # that narrow case; do not hide unrelated missing dependencies.
    if exc.name not in {"oqp", "oqp.utils", "oqp.utils.log_format"}:
        raise
    import importlib.util
    from pathlib import Path

    _format_path = Path(__file__).with_name("log_format.py")
    _format_spec = importlib.util.spec_from_file_location(
        "oqp_log_format_for_state_labels", _format_path)
    if _format_spec is None or _format_spec.loader is None:
        raise
    _format_module = importlib.util.module_from_spec(_format_spec)
    _format_spec.loader.exec_module(_format_module)
    format_log_fields = _format_module.format_log_fields


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

_DFTB_TYPE_ALIASES = {
    "auto": "auto",
    "ground": "ground",
    "dftb": "ground",
    "dftb0": "ground_noscc",
    "noscc": "ground_noscc",
    "ground_noscc": "ground_noscc",
    "tddftb": "tddftb",
    "tda": "tddftb",
    "td-dftb": "tddftb",
    "rpa": "tddftb",
    "sf": "sf",
    "sftddftb": "sf",
    "sf-tddftb": "sf",
    "mrsf": "mrsf",
    "mrsftddftb": "mrsf",
    "mrsf-tddftb": "mrsf",
}

_DFTB_METHOD_NAMES = {
    "ground": "DFTB",
    "ground_noscc": "DFTB0",
    # The openqp-dftb conventional response kernel implements the TDA
    # eigenproblem for both accepted conventional route spellings.
    "tddftb": "TD-DFTB (TDA)",
    "sf": "SF-TDDFTB",
    "mrsf": "MRSF-TDDFTB",
}

# Additive openqp-dftb C capability-mask bits.  These assignments are part of
# the public C contract and deliberately do not follow the state-gradient ABI
# version: an optional feature may be added without changing that argument
# list.
DFTB_CAP_STRUCTURED_TRACE = 1
DFTB_CAP_STATE_SPECTRUM = 2
DFTB_CAP_SCC_FINAL_TRUST = 4

# Named operator presets published by openqp-dftb ([dftb] model=...).  The
# parameter vectors themselves live in openqp-dftb (src/openqp_dftb_presets.F90,
# the single source of truth) and are reported back through the
# openqp_dftb_resolved_options record; this tuple only feeds the settings log.
# Canonical spellings only -- the legacy aliases ("dtcam-tb", "dtcam-tb2",
# "dtcam-tb-erf", "dftb+") still parse but are no longer advertised.
DFTB_KNOWN_PRESETS = ("dtcam", "dtcam2", "dtcam-erf", "ob2")


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


def canonical_dftb_type(value):
    """Return the backend spelling shared by API, logs, and runtime dispatch."""
    normalized = str(value or "auto").strip().lower()
    return _DFTB_TYPE_ALIASES.get(normalized, normalized)


def resolved_dftb_type(config):
    """Resolve ``[dftb] type=auto`` using the same rules as the adapter."""
    explicit = canonical_dftb_type(_section(config, "dftb").get("type", "auto"))
    if explicit != "auto":
        return explicit

    inp = _section(config, "input")
    runtype = str(inp.get("runtype", "energy")).strip().lower()
    istate = _int(_section(config, "optimize").get("istate", 0))
    if runtype in {"optimize", "mep"} and istate == 0:
        return "ground"

    td_type = response_type(config)
    if runtype in {"energy", "grad"} and td_type in {"", "rpa", "tda"}:
        grad_states = _int_list(_section(config, "properties").get("grad", []))
        if not any(state > 0 for state in grad_states):
            return "ground"

    return {
        "rpa": "tddftb",
        "tda": "tddftb",
        "sf": "sf",
        "mrsf": "mrsf",
    }.get(td_type, td_type or "ground")


def dftb_method_name(config):
    """Return the scientifically accurate display name for a DFTB route."""
    return _DFTB_METHOD_NAMES.get(resolved_dftb_type(config), "DFTB")


def dftb_reference_description(config):
    """Describe the reference actually requested from openqp-dftb."""
    method = resolved_dftb_type(config)
    dftb = _section(config, "dftb")
    explicit = _int(dftb.get("reference_multiplicity", 0))
    if method in {"sf", "mrsf"}:
        # The native backend treats 0/1 as "use the required triplet
        # high-spin reference"; report that effective value in the log.
        multiplicity = explicit if explicit > 1 else 3
    else:
        multiplicity = explicit if explicit > 0 else 1
    spin = spin_name(multiplicity)
    if method in {"sf", "mrsf"}:
        return "%s ROKS (internal high-spin reference)" % spin
    if method == "ground_noscc":
        return "%s non-SCC DFTB0 reference" % spin
    if method == "tddftb":
        return "singlet restricted SCC-DFTB reference"
    return "%s SCC-DFTB reference" % spin


def dftb_target_description(config):
    """Describe the physical state manifold produced by the DFTB route."""
    method = resolved_dftb_type(config)
    if method in {"ground", "ground_noscc"}:
        dftb = _section(config, "dftb")
        explicit = _int(dftb.get("reference_multiplicity", 0))
        multiplicity = explicit if explicit > 0 else 1
        prefix = _SPIN_NAMES.get(
            multiplicity, ("multiplicity %d" % multiplicity, "M%d-" % multiplicity)
        )[1]
        return "%s0 ground state" % prefix

    nstate = max(1, _int(_section(config, "tdhf").get("nstate", 1), 1))
    runtype = str(_section(config, "input").get("runtype", "energy")).lower()
    selected_roots = []
    if runtype in {"grad", "data", "prop"}:
        selected_roots = _int_list(_section(config, "properties").get("grad", []))
    elif runtype in {"optimize", "mep"}:
        selected_roots = [_int(_section(config, "optimize").get("istate", 0))]
    elif runtype == "meci":
        optimize = _section(config, "optimize")
        selected_roots = [
            _int(optimize.get("istate", 1)),
            _int(optimize.get("jstate", 2)),
        ]

    if method == "tddftb":
        if selected_roots:
            return ", ".join("S%d" % root for root in selected_roots)
        return "singlet S1" if nstate == 1 else "singlet S1-S%d" % nstate
    if method == "sf":
        if selected_roots:
            roots = ", ".join("root %d" % root for root in selected_roots)
        else:
            roots = "root 1" if nstate == 1 else "roots 1-%d" % nstate
        return "spin-flip %s (character assigned after calculation)" % roots

    targets = requested_states(config)
    if targets:
        return targets
    multiplicity = _int(_section(config, "dftb").get("target_multiplicity", 1), 1)
    prefix = _SPIN_NAMES.get(
        multiplicity, ("multiplicity %d" % multiplicity, "M%d-" % multiplicity)
    )[1]
    return "%s0" % prefix if nstate == 1 else "%s0-%s%d" % (
        prefix, prefix, nstate - 1
    )


def _yes_no(value):
    return "yes" if bool(value) else "no"


def _display_number(value):
    try:
        return "%.8g" % float(value)
    except (TypeError, ValueError):
        return str(value)


def format_dftb_settings(
    config,
    *,
    backend=None,
    parameter_path=None,
    library_path=None,
    executable=None,
    abi_version=None,
    capabilities=None,
):
    """Format one concise, grouped summary of effective DFTB request settings.

    A named model preset is resolved inside openqp-dftb.  In that case this
    formatter deliberately does not present schema defaults as effective
    operator/SCC values.
    """
    dftb = _section(config, "dftb")
    td = _section(config, "tdhf")
    method = resolved_dftb_type(config)
    requested_backend = str(dftb.get("backend", "native")).strip().lower()
    actual_backend = str(backend or requested_backend).strip().lower()
    if actual_backend == "probe":
        # The command-line probe exposes energies/gradients only.  It has no
        # C capability symbol and does not forward the native diagnostic
        # environment hooks.
        effective_capabilities = 0
    elif capabilities is not None:
        effective_capabilities = int(capabilities)
    elif abi_version is not None:
        # Once a native library has been loaded, a missing optional symbol is
        # authoritatively the zero mask, including for ABI 3 builds predating
        # capability discovery.
        effective_capabilities = 0
    else:
        effective_capabilities = None

    def has_capability(bit):
        return (
            effective_capabilities is not None
            and bool(effective_capabilities & bit)
        )

    def unavailable_capability(label):
        if actual_backend == "probe":
            return "unavailable with probe backend"
        if effective_capabilities is None:
            return "availability checked after native library load"
        abi = " (C ABI %s)" % abi_version if abi_version is not None else ""
        return "unavailable; loaded library%s does not advertise %s" % (
            abi, label)

    backend_text = actual_backend
    if requested_backend != actual_backend:
        backend_text += " (requested: %s)" % requested_backend
    model = str(dftb.get("model", "")).strip()
    if model.lower() == "none":
        # Explicit opt-out of the SF/MRSF dtcam-tb default.
        model = ""

    groups = []
    method_rows = [
        ("Method", dftb_method_name(config)),
        ("Backend", backend_text),
        ("Parameters", parameter_path or dftb.get("parameter_path") or "bundled default"),
        ("Reference", dftb_reference_description(config)),
        ("Target", dftb_target_description(config)),
    ]
    if actual_backend == "native":
        level = dftb.get("print_level", 1)
        if effective_capabilities is None or has_capability(
                DFTB_CAP_STRUCTURED_TRACE):
            trace_text = "level %s" % level
        else:
            trace_text = "requested level %s; %s" % (
                level,
                unavailable_capability("structured progress tracing"),
            )
        method_rows.insert(2, ("Native progress trace", trace_text))
    if library_path:
        method_rows.insert(2, ("Shared library", library_path))
    if abi_version is not None:
        method_rows.insert(3, ("C API ABI", abi_version))
    if actual_backend == "native" and effective_capabilities is not None:
        names = []
        for bit, name in (
            (DFTB_CAP_STRUCTURED_TRACE, "structured trace"),
            (DFTB_CAP_STATE_SPECTRUM, "state spectrum"),
            (DFTB_CAP_SCC_FINAL_TRUST, "SCC final trust recovery"),
        ):
            if effective_capabilities & bit:
                names.append(name)
        capability_text = ", ".join(names) if names else "none advertised"
        method_rows.insert(4, ("C API capabilities", capability_text))
    if executable:
        method_rows.insert(2, ("Probe executable", executable))
    groups.append(("calculation", method_rows))

    if method == "ground_noscc":
        scc_rows = [("Enabled", "no (DFTB0)")]
    elif model:
        scc_rows = [
            ("Enabled", "yes"),
            ("Protocol", "%s preset" % model),
        ]
    else:
        mixer = str(dftb.get("scc_mixer", "auto")).strip().lower()
        scc_rows = [
            ("Enabled", "yes"),
            ("Mixer", "%s (requested)" % mixer),
            ("Tolerance", _display_number(dftb.get("scc_tolerance", 1.0e-8))),
            ("Maximum iterations", dftb.get("max_scc_iterations", 1200)),
            ("Mixing / history / max step", "%s / %s / %s" % (
                _display_number(dftb.get("scc_mixing", 0.35)),
                dftb.get("scc_history", 12),
                _display_number(dftb.get("scc_max_step", 0.5)),
            )),
        ]
    if method != "ground_noscc":
        if has_capability(DFTB_CAP_SCC_FINAL_TRUST):
            recovery = "final charge/spin trust-TRAH pass when eligible"
        else:
            recovery = unavailable_capability("final trust-region SCC recovery")
        scc_rows.append(("Failed-cycle recovery", recovery))
    groups.append(("SCC", scc_rows))

    if method not in {"ground", "ground_noscc"}:
        solver = (
            "%s preset" % model
            if model else
            "%s (requested)" % str(dftb.get("response_solver", "auto")).lower()
        )
        groups.append(("response", [
            ("Solver", solver),
            ("States", td.get("nstate", 1)),
            ("Tolerance", _display_number(dftb.get("response_tolerance", 1.0e-6))),
            ("Maximum iterations / subspace", "%s / %s" % (
                dftb.get("response_max_iterations", 50),
                dftb.get("response_max_subspace", 100),
            )),
        ]))
        if method in {"sf", "mrsf"}:
            groups[-1][1].append(
                ("Spin completion", _yes_no(dftb.get("spin_complete", True)))
            )
        else:
            groups[-1][1].append(("Response manifold", "closed-shell singlet TDA"))
        spectrum_requested = bool(dftb.get("state_to_state_spectrum", True))
        if not spectrum_requested:
            spectrum_text = "no"
        elif has_capability(DFTB_CAP_STATE_SPECTRUM):
            spectrum_text = "yes (unrelaxed TDA/state interaction)"
        else:
            spectrum_text = "requested; %s" % unavailable_capability(
                "all-pair state spectrum"
            )
        groups[-1][1].append(("State-to-state spectrum", spectrum_text))
        zvector = (
            "%s preset" % model
            if model else _yes_no(dftb.get("zvector", True))
        )
        groups.append(("gradient", [("Relaxed Z-vector", zvector)]))

    if model:
        operator_rows = [
            ("Model preset", "%s (operator resolved by openqp-dftb backend)" % model),
            ("Effective values", "reported under 'DFTB effective options'"),
        ]
    elif method == "ground_noscc":
        operator_rows = [("Hamiltonian", "non-SCC DFTB0")]
    else:
        operator_rows = [
            ("Model preset", "none (explicit [dftb] keys; schema defaults)"),
            ("LC ground state", _yes_no(dftb.get("lc_ground_state", False))),
        ]
        if method in {"sf", "mrsf"} or bool(dftb.get("lc_ground_state", False)):
            operator_rows.extend([
                ("Gamma kernel / omega", "%s / %s" % (
                    str(dftb.get("lc_gamma", "yukawa")).lower(),
                    _display_number(dftb.get("omega", 0.3)),
                )),
                ("CAM alpha / beta", "%s / %s" % (
                    _display_number(dftb.get("cam_alpha", 0.0)),
                    _display_number(dftb.get("cam_beta", 1.0)),
                )),
            ])
        if method in {"sf", "mrsf"}:
            spc = dftb.get("spc", 0.5)
            channels = []
            for key in ("spc_coco", "spc_ovov", "spc_coov"):
                value = dftb.get(key, -999.0)
                try:
                    value = spc if float(value) <= -900.0 else value
                except (TypeError, ValueError):
                    pass
                channels.append(_display_number(value))
            operator_rows.append(("SPC coco / ovov / coov", " / ".join(channels)))
            if method == "mrsf":
                operator_rows.append(("MRSF shifts oo/co/ov/cv", " / ".join(
                    _display_number(dftb.get(key, 0.0))
                    for key in (
                        "mrsf_shift_oo", "mrsf_shift_co",
                        "mrsf_shift_ov", "mrsf_shift_cv",
                    )
                )))
    if method != "ground_noscc":
        operator_rows.append(
            ("Available presets", ", ".join(DFTB_KNOWN_PRESETS))
        )
    groups.append(("operator", operator_rows))

    titles = {
        "calculation": "Calculation",
        "SCC": "SCC",
        "response": "Response",
        "gradient": "Gradient",
        "operator": "Operator",
    }
    blocks = []
    for group, rows in groups:
        lines = ["   %s" % titles.get(group, group.capitalize())]
        lines.append(format_log_fields(
            rows, prefix=None, label_width=29, indent=5))
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks)


def is_mrsf(config):
    """True for MRSF/UMRSF response calculations, including TD-DFTB."""
    if str(_section(config, "input").get("method", "")).strip().lower() == "dftb":
        return resolved_dftb_type(config) == "mrsf"
    return response_type(config) in {"mrsf", "umrsf"}


def spin_name(multiplicity):
    """Return a readable spin name, retaining unusual multiplicities."""
    mult = _int(multiplicity, 1)
    return _SPIN_NAMES.get(mult, ("multiplicity %d" % mult, "M%d-" % mult))[0]


def _target_multiplicity(config):
    """Return the multiplicity consumed by the active backend."""
    if (
        str(_section(config, "input").get("method", "")).strip().lower() == "dftb"
        and resolved_dftb_type(config) == "mrsf"
    ):
        return _int(_section(config, "dftb").get("target_multiplicity", 1), 1)
    return _int(_section(config, "tdhf").get("multiplicity", 1), 1)


def public_method_name(config):
    """Return the method name a user recognizes, not its legacy dispatch key."""
    inp = _section(config, "input")
    method = str(inp.get("method", "hf")).strip().lower()
    td_type = response_type(config)

    if method == "dftb":
        return dftb_method_name(config)
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
    if method in ("ccsd", "ccsd(t)"):
        # Without this these fall through to the empty-functional case below
        # and a correlated calculation reports itself as plain HF.
        return method.upper()
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
    mult = _target_multiplicity(config) if multiplicity is None else _int(multiplicity, 1)
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
    mult = _target_multiplicity(config)
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
    td_mult = _target_multiplicity(config)

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
        search = str(opt.get("meci_search", "auto")).strip().lower()
        roots = (
            _int_list(opt.get("states", []))
            if search in {"auto", "baeka"} else []
        )
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
    auto_meci = (
        runtype == "meci"
        and str(_section(config, "optimize").get("meci_search", "")).strip().lower()
        == "auto"
    )
    if baeka:
        calculation = "BaekA multistate conical intersection"
    elif auto_meci:
        calculation = "Automatic conical intersection search"
    if runtype == "nac" and bool(_section(config, "nac").get("bp", False)):
        calculation = "Branching-plane analysis"
    lines = [
        ("Method", public_method_name(config)),
        ("Calculation", calculation),
    ]
    if baeka:
        lines.append(("Algorithm", "Baek adaptive penalty (BaekA)"))
    elif auto_meci:
        lines.append(("Algorithm", "penalty → BaekA on nonconvergence"))
    is_dftb_request = str(inp.get("method", "")).strip().lower() == "dftb"
    if not is_dftb_request:
        functional = str(inp.get("functional", "")).strip()
        basis = str(inp.get("basis", "")).strip()
        if functional and basis:
            lines.append(("Model chemistry", "%s/%s" % (functional, basis)))
        elif basis:
            lines.append(("Basis", basis))

    targets = dftb_target_description(config) if is_dftb_request else requested_states(config)
    if targets:
        lines.append(("Physical target state(s)", targets))
    freeze = str(_section(config, "oqp").get("freeze", "")).strip()
    if freeze:
        lines.append(("Frozen distance(s)", freeze))
    if is_dftb_request:
        lines.append(("Reference", dftb_reference_description(config)))
    if is_mrsf(config):
        target_mult = _target_multiplicity(config)
        runtype = str(inp.get("runtype", "energy")).lower()
        multi_spin = runtype in {"soc", "mecp"} or (
            runtype == "namd" and bool(_section(config, "md").get("soc", False))
        )
        if not multi_spin:
            lines.append(("Target spin", spin_name(target_mult)))
        if not is_dftb_request:
            lines.append(("Reference", reference_description(config)))
        lines.append(("State labels", "physical labels shown; engine root numbers are internal"))
    if source:
        lines.append(("Input", str(source)))
    if resolved and str(resolved) != str(source):
        lines.append(("Resolved input", str(resolved)))
    return lines


def format_calculation_request(config, source=None, resolved=None):
    """Format :func:`calculation_request_lines` as a log-ready text block."""
    return format_log_fields(calculation_request_lines(
        config, source=source, resolved=resolved))
