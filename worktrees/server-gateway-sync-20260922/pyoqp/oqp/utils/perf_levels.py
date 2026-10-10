"""
Centralised resolver for OpenQP's performance options and the ``perf`` preset.

Each performance knob is an ordinary **input key** (e.g. ``[scf] xc_c2f``,
``[tdhf] resp_cutoff``) backed by a field on the shared ``control`` struct -- there are
no environment variables involved. ``[input] perf = 0|1|2|3`` is a preset that fills in
those input keys (the ones the user left at the sentinel ``auto``) before the config is
pushed to the native library. An explicit input key always overrides the preset.

  perf = 0  strict reference / reproducible  -- every accelerator off, tightest cutoffs.
  perf = 1  recommended production (exact)   -- only exact, proven-helpful knobs:
                                                MRSF response cutoff 1e-8, z-vector
                                                warm-start (+ the always-on Fock digest).
  perf = 2  faster, tiny degradation         -- + coarse-to-fine XC grid + grad Schwarz
                                                cutoff 1e-8 (<=5e-7 a.u. on gradients).
  perf = 3  aggressive, degradation allowed  -- looser cutoffs traded for speed: grad
                                                1e-7 (~1e-5 a.u.), response 1e-6 (~few ueV).

The accelerators that the CPU benchmark found net-negative or experimental -- Phi-cache,
IncDFT, FP32 and progressive screening (pscreen) -- are NOT enabled by any preset (they
stay available as explicit input keys, e.g. for GPU XC or regimes not covered here).

Precedence (low -> high): control-struct default  <  perf preset  <  explicit input key.
An input key set to ``auto`` (its default) defers to the preset; any other value overrides
it. ``perf`` unset (-1) leaves every key at its default -> identical to legacy behaviour.
"""

UNSET = -1

# Per-knob preset table.  Columns: section, input-key, [value@0, @1, @2, @3].
# Every key here is a sentinel-'auto' string input key: 'auto'/'' means "defer to the
# preset"; any other value is an explicit override (parsed by the setter in oqpdata.py).
KNOBS = [
    ("scf",  "xc_c2f",        ["off", "off", "on",  "on" ]),
    ("scf",  "xc_phi_cache",  ["off", "off", "off", "off"]),
    ("scf",  "xc_incdft",     ["off", "off", "off", "off"]),
    ("scf",  "grad_cutoff",   ["1.0d-10", "1.0d-10", "1.0d-8", "1.0d-7"]),
    ("tdhf", "resp_cutoff",   ["5e-11", "1e-8", "1e-8", "1e-6"]),
    ("tdhf", "fp32",          ["off", "off", "off", "off"]),
    ("tdhf", "zv_warmstart",  ["off", "on",  "on",  "on" ]),
]

# Input keys this module introduces (added to OQP_CONFIG_SCHEMA as sentinel 'auto' strings).
NEW_INPUT_KEYS = {
    "scf":  ["xc_c2f", "xc_phi_cache", "xc_incdft", "grad_cutoff"],
    "tdhf": ["resp_cutoff", "fp32", "zv_warmstart"],
}


def _truthy(v):
    return str(v).strip().lower() in ("1", "y", "yes", "t", "true", "on")


def _is_auto(v):
    return v is None or str(v).strip().lower() in ("", "auto")


def resolve(config, perf):
    """Fill the perf input keys in ``config`` from the preset (in place).

    A key is filled from the preset only when its value is the sentinel ``auto`` -- so an
    explicit value (from the input file or a caller-supplied dict, even a fully
    default-expanded one) always wins. Returns a list of (label, shown_value, source) with
    source in {"input", "perf"}. ``perf`` unset -> no change, empty report.
    """
    if perf is None or perf == UNSET:
        return []
    if perf not in (0, 1, 2, 3):
        raise ValueError("perf=%r is invalid; choose 0, 1, 2 or 3" % perf)
    report = []
    for sec, key, vals in KNOBS:
        cfg_val = config.get(sec, {}).get(key)
        if not _is_auto(cfg_val):                         # explicit value wins
            report.append(("%s.%s" % (sec, key), cfg_val, "input"))
            continue
        config.setdefault(sec, {})[key] = vals[perf]      # else apply the preset
        report.append(("%s.%s" % (sec, key), vals[perf], "perf"))
    return report


def validate(config, scf_conv=None, zv_conv=None):
    """Warnings for risky resolved settings (reads the post-resolution config)."""
    warns = []
    scf = config.get("scf", {})
    tdhf = config.get("tdhf", {})
    if _truthy(scf.get("xc_incdft")):
        warns.append("xc_incdft is experimental and usually *slows* SCF convergence.")
    if _truthy(scf.get("pscreen")):
        cap = scf.get("pscreen_cap")
        try:
            if cap is not None and float(str(cap).replace("d", "e")) > 1e-8:
                warns.append("pscreen_cap looser than 1e-8 risks a spurious SCF state.")
        except ValueError:
            pass
    if _truthy(tdhf.get("fp32")):
        warns.append("fp32 is non-reproducible, can flip near-degenerate states, and was "
                     "net-slower than fp64 on CPU in the perf benchmark.")
        for label, c in (("scf.conv", scf_conv), ("tdhf.zvconv", zv_conv)):
            if c is not None and c < 1e-9:
                warns.append("fp32 with a tight %s (%g) may stall before convergence." % (label, c))
    return warns


def format_report(perf, report, warns):
    if not report:
        return ""
    head = "perf = %s" % ("unset" if perf in (None, UNSET) else perf)
    lines = ["", "   Performance settings (%s)" % head]
    w = max(len(lbl) for lbl, _, _ in report)
    src = {"input": "(input)", "perf": "(preset)"}
    for lbl, val, s in report:
        lines.append("     %-*s = %-9s %s" % (w, lbl, val, src.get(s, "")))
    for warn in warns:
        lines.append("     WARNING: " + warn)
    lines.append("")
    return "\n".join(lines)


def apply(config, perf, scf_conv=None, zv_conv=None):
    """Resolve into ``config`` (in place) + validate; return (report, warns)."""
    report = resolve(config, perf)
    warns = validate(config, scf_conv=scf_conv, zv_conv=zv_conv) if report else []
    return report, warns
