"""Parse and pretty-print openqp-dftb structured native diagnostics.

The openqp-dftb shared library cannot write into the OpenQP log directly (the
log unit belongs to liboqp's Fortran runtime), so the adapter captures the
library's flushed stdout.  At trace level 2 that capture carries one
machine-readable record per SCC/Davidson/Z-vector event
(``openqp_dftb_scc_iter ...``).  This module turns those records into the same
kind of human-readable iteration tables OpenQP prints for its native SCF,
Davidson, and Z-vector solvers, so a DFTB log reads like any other OpenQP log.
"""

from __future__ import annotations

import re


HARTREE_TO_EV = 27.211386245988

# Update-kind display names for the per-iteration SCC labels.  The record label
# is "<pass label>[_<update kind>]" (e.g. rohf_trah is a TRAH step inside the
# rohf pass).
_SCC_ITER_MODES = {
    "": "",
    "trah": "TRAH",
    "trust": "trust",
    "ls": "level-shift",
    "ls_trust": "level-shift+trust",
}

_SCC_PASS_TITLES = {
    "rhf_seed": "restricted charge seed",
    "rohf": "open-shell (ROKS) reference",
    "restricted": "restricted reference",
}

_SPECTRUM_HEADER = "DFTB state-to-state oscillator strengths"
# Old printers wrote 3 numbers per row (dE, |mu|, f); current ones APPEND the
# dipole components (dE, |mu|, f, mu_x, mu_y, mu_z) so positional parsers of
# the historical layout keep reading the same fields.
_SPECTRUM_ROW = re.compile(
    r"^\s*(S\d+|root \d+)\s+(S\d+|root \d+)((?:\s+[-+0-9.EeDd]+){3,6})\s*$"
)
_LINEAR_RECORD = re.compile(
    r"^openqp_dftb_(gmres|minres|bicgstab|cg)_(start|iter|end)$"
)


def _float(token):
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except (AttributeError, ValueError):
        return None


def _int(token):
    try:
        return int(token)
    except (TypeError, ValueError):
        return None


def _bool(token):
    return str(token).strip() in {"T", "t", "true", "True", "1"}


def _pairs(tokens):
    """Map ``key value key value ...`` tokens; last value wins on repeats."""
    mapping = {}
    for key, value in zip(tokens[0::2], tokens[1::2]):
        mapping[key] = value
    return mapping


def _typed_value(token):
    """Best-effort typed record value: bool ('T'/'F'), int, float, or str."""
    if token in ("T", "F"):
        return token == "T"
    if re.fullmatch(r"[-+]?\d+", token or ""):
        return int(token)
    number = _float(token)
    return number if number is not None else token


def parse_native_trace(text):
    """Split a captured native stdout into structured events.

    Returns a dict with ``scc_passes``, ``post_updates``, ``recoveries``,
    ``davidson``, ``zvector``, ``zvector_dense``, ``spectrum``, ``warnings``,
    and ``other`` (verbatim unrecognized lines, e.g. ``[zvec-stage]`` timing
    breakdowns).  The parser is tolerant: an unknown or malformed record line
    degrades to ``other`` instead of failing the calculation log.
    """
    trace = {
        "scc_passes": [],
        "post_updates": [],
        "recoveries": [],
        "davidson": [],
        "zvector": [],
        "zvector_dense": [],
        "energy_components": None,
        "resolved_options": None,
        "spectrum": None,
        "linear": [],
        "warnings": [],
        "other": [],
    }
    open_scc = {}
    open_davidson = None
    open_zvector = None
    open_dense = None
    open_linear = None
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        line = lines[index]
        index += 1
        stripped = line.strip()
        if not stripped:
            continue
        if stripped == _SPECTRUM_HEADER:
            spectrum, index = _parse_spectrum(lines, index)
            trace["spectrum"] = spectrum
            continue
        tokens = stripped.split()
        record = tokens[0]
        try:
            if record == "openqp_dftb_scc_pass":
                label = tokens[1]
                info = _pairs(tokens[2:])
                entry = {
                    "label": label,
                    "pass": _int(info.get("pass")),
                    "stage": info.get("stage", ""),
                    "mixer": info.get("mixer", ""),
                    "iterations": [],
                    "trah_steps": 0,
                    "outcome": None,
                    "iteration_count": None,
                    "elapsed": None,
                }
                trace["scc_passes"].append(entry)
                open_scc[label] = entry
            elif record == "openqp_dftb_scc_iter":
                label = tokens[1]
                iteration = _int(tokens[2])
                info = _pairs(tokens[3:])
                base, mode = _split_scc_iter_label(label)
                entry = open_scc.get(base)
                if entry is None:
                    entry = {
                        "label": base, "pass": _int(info.get("pass")),
                        "stage": "", "mixer": "", "iterations": [],
                        "trah_steps": 0, "outcome": None,
                        "iteration_count": None, "elapsed": None,
                    }
                    trace["scc_passes"].append(entry)
                    open_scc[base] = entry
                entry["iterations"].append({
                    "iteration": iteration,
                    "residual": _float(info.get("residual")),
                    "charge": _float(info.get("charge")),
                    "spin": _float(info.get("spin")),
                    "orbital": _float(info.get("orbital")),
                    "tolerance": _float(info.get("tol")),
                    "mode": mode,
                })
            elif record == "openqp_dftb_scc_trah_step":
                base, _ = _split_scc_iter_label(tokens[1])
                entry = open_scc.get(base)
                if entry is not None:
                    entry["trah_steps"] += 1
            elif record == "openqp_dftb_scc_pass_end":
                label = tokens[1]
                info = _pairs(tokens[2:])
                entry = open_scc.pop(label, None)
                if entry is None:
                    entry = {
                        "label": label, "pass": _int(info.get("pass")),
                        "stage": "", "mixer": "", "iterations": [],
                        "trah_steps": 0,
                    }
                    trace["scc_passes"].append(entry)
                entry["outcome"] = info.get("outcome")
                entry["iteration_count"] = _int(info.get("iterations"))
                entry["elapsed"] = _float(info.get("elapsed_seconds"))
            elif record == "openqp_dftb_scc_post_update":
                info = _pairs(tokens[2:])
                trace["post_updates"].append({
                    "label": tokens[1],
                    "residual": _float(info.get("residual")),
                    "charge": _float(info.get("charge")),
                    "spin": _float(info.get("spin")),
                    "density": _float(info.get("density")),
                    "tolerance": _float(info.get("tol")),
                    "converged": _bool(info.get("converged")),
                    "elapsed": _float(info.get("elapsed_seconds")),
                })
            elif record == "openqp_dftb_scc_recovery":
                info = _pairs(tokens[2:])
                trace["recoveries"].append({
                    "label": tokens[1],
                    "kind": info.get("kind", ""),
                    "attempted": _bool(info.get("attempted")),
                    "succeeded": _bool(info.get("succeeded")),
                })
            elif record == "openqp_dftb_davidson_start":
                info = _pairs(tokens[1:])
                open_davidson = {
                    "dimension": _int(info.get("dimension")),
                    "roots": _int(info.get("roots")),
                    "max_subspace": _int(info.get("max_subspace")),
                    "max_iterations": _int(info.get("max_iterations")),
                    "tolerance": _float(info.get("tolerance")),
                    "iterations": [],
                    "outcome": None,
                }
                trace["davidson"].append(open_davidson)
            elif record == "openqp_dftb_davidson_iter":
                info = _pairs(tokens[1:])
                if open_davidson is None:
                    open_davidson = {
                        "dimension": None, "roots": None, "max_subspace": None,
                        "max_iterations": None, "tolerance": None,
                        "iterations": [], "outcome": None,
                    }
                    trace["davidson"].append(open_davidson)
                open_davidson["iterations"].append({
                    "iteration": _int(info.get("iteration")),
                    "basis": _int(info.get("basis")),
                    "converged_roots": _int(info.get("converged_roots")),
                    "residual": _float(info.get("residual")),
                    "matvecs": _int(info.get("matvecs")),
                    "restarts": _int(info.get("restarts")),
                })
            elif record == "openqp_dftb_davidson_end":
                info = _pairs(tokens[1:])
                if open_davidson is None:
                    open_davidson = {
                        "dimension": None, "roots": None, "max_subspace": None,
                        "max_iterations": None, "tolerance": None,
                        "iterations": [],
                    }
                    trace["davidson"].append(open_davidson)
                open_davidson.update({
                    "outcome": info.get("outcome"),
                    "iteration_count": _int(info.get("iterations")),
                    "residual": _float(info.get("residual")),
                    "matvecs": _int(info.get("matvecs")),
                    "restarts": _int(info.get("restarts")),
                    "elapsed": _float(info.get("elapsed_seconds")),
                })
                open_davidson = None
            elif record == "openqp_dftb_zvector_start":
                info = _pairs(tokens[1:])
                open_zvector = {
                    "dimension": _int(info.get("dimension")),
                    "max_iterations": _int(info.get("max_iterations")),
                    "tolerance": _float(info.get("tolerance")),
                    "iterations": [],
                    "outcome": None,
                }
                trace["zvector"].append(open_zvector)
            elif record == "openqp_dftb_zvector_iter":
                info = _pairs(tokens[1:])
                if open_zvector is None:
                    open_zvector = {
                        "dimension": None, "max_iterations": None,
                        "tolerance": None, "iterations": [], "outcome": None,
                    }
                    trace["zvector"].append(open_zvector)
                open_zvector["iterations"].append({
                    "iteration": _int(info.get("iteration")),
                    "residual": _float(info.get("residual")),
                    "inner_solves": _int(info.get("inner_solves")),
                    "inner_matvecs": _int(info.get("inner_matvecs")),
                })
            elif record == "openqp_dftb_zvector_end":
                info = _pairs(tokens[1:])
                if open_zvector is None:
                    open_zvector = {
                        "dimension": None, "max_iterations": None,
                        "tolerance": None, "iterations": [],
                    }
                    trace["zvector"].append(open_zvector)
                open_zvector.update({
                    "outcome": info.get("outcome"),
                    "pair_iterations": _int(info.get("pair_iterations")),
                    "anderson_iterations": _int(info.get("anderson_iterations")),
                    "fallback_gmres_iterations": _int(
                        info.get("fallback_gmres_iterations")),
                    "residual": _float(info.get("residual")),
                    "elapsed": _float(info.get("elapsed_seconds")),
                })
                open_zvector = None
            elif record == "openqp_dftb_zvector_dense_start":
                info = _pairs(tokens[1:])
                open_dense = {
                    "dimension": _int(info.get("dimension")),
                    "solver": info.get("solver", "dense_direct"),
                    "stages": [],
                    "outcome": None,
                }
                trace["zvector_dense"].append(open_dense)
            elif record == "openqp_dftb_zvector_dense_stage":
                info = _pairs(tokens[1:])
                if open_dense is not None:
                    open_dense["stages"].append(
                        (info.get("stage", ""), _float(info.get("elapsed_seconds")))
                    )
            elif record == "openqp_dftb_zvector_dense_end":
                info = _pairs(tokens[1:])
                if open_dense is None:
                    open_dense = {
                        "dimension": None, "solver": "dense_direct",
                        "stages": [],
                    }
                    trace["zvector_dense"].append(open_dense)
                open_dense.update({
                    "outcome": info.get("outcome"),
                    "solve_succeeded": _bool(info.get("solve_succeeded")),
                    "residual": _float(info.get("residual")),
                    "elapsed": _float(info.get("elapsed_seconds")),
                })
                open_dense = None
            elif record == "openqp_dftb_energy_components":
                info = _pairs(tokens[1:])
                trace["energy_components"] = {
                    key: _float(value) for key, value in info.items()
                }
            elif record == "openqp_dftb_resolved_options":
                info = _pairs(tokens[1:])
                trace["resolved_options"] = {
                    key: _typed_value(value) for key, value in info.items()
                }
            elif record == "openqp_dftb_state_spectrum_warning":
                trace["warnings"].append(stripped)
            elif _LINEAR_RECORD.match(record):
                match = _LINEAR_RECORD.match(record)
                solver, event = match.group(1), match.group(2)
                info = _pairs(tokens[1:])
                if event == "start":
                    open_linear = {
                        "solver": solver,
                        "dimension": _int(info.get("dimension")),
                        "max_iterations": _int(info.get("max_iterations")),
                        "max_subspace": _int(info.get("max_subspace")),
                        "tolerance": _float(info.get("tolerance")),
                        "iterations": [],
                        "outcome": None,
                    }
                    # A Krylov solve running inside the Z-vector window IS the
                    # Z-vector solve (pair-space single solve / full-space
                    # fallback); attach it there so its iterations can be
                    # rendered as the Z-vector iteration table.
                    if open_zvector is not None:
                        open_zvector.setdefault("krylov_solves", []).append(
                            open_linear)
                    else:
                        trace["linear"].append(open_linear)
                elif event == "iter" and open_linear is not None:
                    open_linear["iterations"].append({
                        "iteration": _int(info.get("iteration")),
                        "cycle": _int(info.get("cycle")),
                        "residual": _float(info.get("residual")),
                        "matvecs": _int(info.get("matvecs")),
                    })
                elif event == "end" and open_linear is not None:
                    open_linear.update({
                        "outcome": info.get("outcome"),
                        "iteration_count": _int(info.get("iterations")),
                        "residual": _float(info.get("residual")),
                        "matvecs": _int(info.get("matvecs")),
                        "elapsed": _float(info.get("elapsed_seconds")),
                    })
                    open_linear = None
            else:
                trace["other"].append(line.rstrip())
        except (IndexError, KeyError, TypeError, ValueError):
            trace["other"].append(line.rstrip())
    return trace


def _split_scc_iter_label(label):
    """Return (pass label, update-kind display) for a per-iteration label."""
    for suffix, display in sorted(
            _SCC_ITER_MODES.items(), key=lambda item: -len(item[0])):
        if suffix and label.endswith("_" + suffix):
            return label[: -len(suffix) - 1], display
    return label, ""


def _parse_spectrum(lines, index):
    """Consume the printed state-to-state block; return (spectrum, new index)."""
    spectrum = {"approximation": "", "rows": []}
    while index < len(lines):
        stripped = lines[index].strip()
        if stripped.startswith("Approximation:"):
            spectrum["approximation"] = stripped[len("Approximation:"):].strip()
            index += 1
            continue
        if stripped.startswith("Initial"):
            index += 1
            continue
        match = _SPECTRUM_ROW.match(lines[index])
        if match:
            numbers = [_float(token) for token in match.group(3).split()]
            row = {
                "initial": match.group(1),
                "final": match.group(2),
                "de_ev": numbers[0],
                "dipole": numbers[1] if len(numbers) > 1 else None,
                "oscillator": numbers[2] if len(numbers) > 2 else None,
                "dipole_xyz": numbers[3:6] if len(numbers) >= 6 else None,
            }
            spectrum["rows"].append(row)
            index += 1
            continue
        if not stripped:
            index += 1
            # A blank line ends the block unless the next line is still a row.
            if index < len(lines) and _SPECTRUM_ROW.match(lines[index]):
                continue
            break
        break
    return spectrum, index


def _fmt_seconds(value):
    if value is None:
        return ""
    return " (%.3f s)" % value


def _fmt_exp(value, width=16):
    if value is None:
        return " " * width
    return ("%" + str(width) + ".6E") % value


def format_scc_pass(entry):
    """One OpenQP-style iteration table for a single SCC pass."""
    title = _SCC_PASS_TITLES.get(entry["label"], entry["label"])
    mixer = entry.get("mixer") or ""
    qualifier = []
    if mixer:
        qualifier.append("%s mixer" % mixer)
    if entry.get("pass") and entry["pass"] > 1:
        qualifier.append("pass %d" % entry["pass"])
    if entry.get("stage") and entry["stage"] not in {"initial", "seed"}:
        qualifier.append("stage %s" % entry["stage"])
    header = "   SCC iterations: %s%s" % (
        title, " (%s)" % ", ".join(qualifier) if qualifier else "")

    lines = [header]
    iterations = entry.get("iterations") or []
    if iterations:
        has_components = any(
            it.get("charge") is not None for it in iterations)
        rule = "   " + "-" * (76 if has_components else 42)
        lines.append(rule)
        if has_components:
            lines.append("    Iter        Residual          Charge"
                         "            Spin          Step")
        else:
            lines.append("    Iter        Residual")
        lines.append(rule)
        for it in iterations:
            row = "   %5s" % (it.get("iteration") if it.get("iteration")
                              is not None else "?")
            row += _fmt_exp(it.get("residual"))
            if has_components:
                row += _fmt_exp(it.get("charge"))
                row += _fmt_exp(it.get("spin"))
                mode = it.get("mode") or ""
                row += ("    %s" % mode if mode else "")
            lines.append(row.rstrip())
        lines.append(rule)
    outcome = entry.get("outcome")
    count = entry.get("iteration_count")
    tolerance = None
    if iterations:
        tolerance = iterations[-1].get("tolerance")
    if outcome == "converged":
        summary = "   SCC converged in %s iterations" % (
            count if count is not None else len(iterations))
    elif outcome == "exhausted":
        summary = "   SCC did NOT converge within %s iterations" % (
            count if count is not None else len(iterations))
    elif outcome:
        summary = "   SCC ended (%s) after %s iterations" % (
            outcome, count if count is not None else len(iterations))
    else:
        summary = "   SCC pass recorded %d iterations" % len(iterations)
    if tolerance is not None:
        summary += ", tolerance %.2E" % tolerance
    summary += _fmt_seconds(entry.get("elapsed"))
    lines.append(summary)
    if entry.get("trah_steps"):
        lines.append("   TRAH trust-region steps accepted: %d"
                     % entry["trah_steps"])
    return "\n".join(lines)


def format_scc_post_update(entry):
    parts = ["   Final self-consistency: residual %.2E" % entry["residual"]
             if entry.get("residual") is not None else
             "   Final self-consistency:"]
    components = []
    for key in ("charge", "spin", "density"):
        if entry.get(key) is not None:
            components.append("%s %.2E" % (key, entry[key]))
    if components:
        parts.append("(%s)" % ", ".join(components))
    if entry.get("tolerance") is not None:
        parts.append("tolerance %.2E" % entry["tolerance"])
    parts.append("-- converged" if entry.get("converged")
                 else "-- NOT converged")
    return " ".join(parts)


def format_scc_recoveries(recoveries, verbose=False):
    lines = []
    for entry in recoveries:
        if not entry.get("attempted") and not verbose:
            continue
        outcome = ("succeeded" if entry.get("succeeded") else
                   "did not succeed") if entry.get("attempted") else "not needed"
        lines.append("   SCC recovery (%s): %s" % (entry.get("kind"), outcome))
    return "\n".join(lines)


def format_scc_block(trace, verbose=False):
    """All SCC passes + final residual check + recoveries as one text block."""
    blocks = [format_scc_pass(entry) for entry in trace.get("scc_passes", [])]
    for entry in trace.get("post_updates", []):
        blocks.append(format_scc_post_update(entry))
    recovery_text = format_scc_recoveries(trace.get("recoveries", []),
                                          verbose=verbose)
    if recovery_text:
        blocks.append(recovery_text)
    return "\n\n".join(block for block in blocks if block)


def format_scc_summary_lines(trace):
    """One line per SCC pass -- used when the reference is merely re-solved."""
    lines = []
    for entry in trace.get("scc_passes", []):
        title = _SCC_PASS_TITLES.get(entry["label"], entry["label"])
        outcome = entry.get("outcome") or "recorded"
        count = entry.get("iteration_count")
        if count is None:
            count = len(entry.get("iterations") or [])
        lines.append("   SCC %s: %s in %s iterations%s" % (
            title, outcome, count, _fmt_seconds(entry.get("elapsed"))))
    return "\n".join(lines)


def format_davidson_block(solve, method_display="DFTB"):
    lines = ["   Davidson algorithm for %s response states" % method_display]
    params = []
    for key, label in (("dimension", "dimension"), ("roots", "roots"),
                       ("max_subspace", "max subspace"),
                       ("max_iterations", "max iterations")):
        if solve.get(key) is not None:
            params.append("%s = %s" % (label, solve[key]))
    if solve.get("tolerance") is not None:
        params.append("tolerance = %.1E" % solve["tolerance"])
    if params:
        lines.append("   " + "    ".join(params))
    iterations = solve.get("iterations") or []
    if iterations:
        rule = "   " + "-" * 58
        lines.append("")
        lines.append("    Iter    Subspace    Converged roots      Max residual")
        lines.append(rule)
        for it in iterations:
            lines.append("   %5s   %7s    %11s     %s" % (
                it.get("iteration"), it.get("basis"),
                it.get("converged_roots"),
                _fmt_exp(it.get("residual"))))
        lines.append(rule)
    outcome = solve.get("outcome")
    count = solve.get("iteration_count")
    if outcome == "converged":
        summary = "   Davidson converged in %s iterations" % count
    elif outcome:
        summary = "   Davidson %s after %s iterations" % (outcome, count)
    else:
        summary = "   Davidson recorded %d iterations" % len(iterations)
    details = []
    if solve.get("residual") is not None:
        details.append("residual %.2E" % solve["residual"])
    if solve.get("matvecs") is not None:
        details.append("%d matvecs" % solve["matvecs"])
    if solve.get("restarts"):
        details.append("%d restarts" % solve["restarts"])
    if details:
        summary += ": " + ", ".join(details)
    summary += _fmt_seconds(solve.get("elapsed"))
    lines.append(summary)
    return "\n".join(lines)


def format_davidson_summary_line(solve):
    outcome = solve.get("outcome") or "recorded"
    count = solve.get("iteration_count",
                      len(solve.get("iterations") or []))
    extra = ""
    if solve.get("residual") is not None:
        extra = " (residual %.2E)" % solve["residual"]
    return "   Response Davidson: %s in %s iterations%s" % (
        outcome, count, extra)


def collapse_repeats(lines):
    """Merge consecutive identical lines into one line with a repeat count."""
    merged = []
    for line in lines:
        if merged and merged[-1][0] == line:
            merged[-1][1] += 1
        else:
            merged.append([line, 1])
    return [line if count == 1 else "%s  [x%d]" % (line, count)
            for line, count in merged]


_ZVECTOR_OUTCOMES = {
    "pair_space_converged": "converged (pair-space single solve)",
    "anderson_converged": "converged (Anderson fixed point)",
    "gmres_converged": "converged (GMRES)",
    "fallback_per_perturbation":
        "fell back to the per-perturbation gradient path",
    "converged": "converged",
}


def format_zvector_block(solve):
    lines = ["   Z-vector relaxed-density solve"]
    params = []
    if solve.get("dimension") is not None:
        params.append("dimension = %s" % solve["dimension"])
    if solve.get("max_iterations") is not None:
        params.append("max iterations = %s" % solve["max_iterations"])
    if solve.get("tolerance") is not None:
        params.append("tolerance = %.1E" % solve["tolerance"])
    if params:
        lines.append("   " + "    ".join(params))
    iterations = solve.get("iterations") or []
    if iterations:
        rule = "   " + "-" * 58
        lines.append("")
        lines.append("    Iter        Residual        Inner solves"
                     "    Inner matvecs")
        lines.append(rule)
        for it in iterations:
            lines.append("   %5s%s   %9s        %9s" % (
                it.get("iteration"), _fmt_exp(it.get("residual")),
                it.get("inner_solves") if it.get("inner_solves") is not None
                else "", it.get("inner_matvecs")
                if it.get("inner_matvecs") is not None else ""))
        lines.append(rule)
    else:
        # Production path: one pair-space Krylov solve, no outer fixed-point
        # loop.  Its per-iteration residuals ARE the Z-vector iterations; any
        # further Krylov solves inside the window (auxiliary contractions) are
        # summarized in one line each.
        krylov_solves = solve.get("krylov_solves") or []
        for position, krylov in enumerate(krylov_solves):
            krylov_iterations = _dedupe_iterations(
                krylov.get("iterations") or [])
            if position == 0 and krylov_iterations:
                rule = "   " + "-" * 40
                lines.append("")
                lines.append("   Pair-space %s iterations"
                             % krylov.get("solver", "gmres").upper())
                lines.append("    Iter        Residual")
                lines.append(rule)
                for it in krylov_iterations:
                    lines.append("   %5s%s" % (
                        it.get("iteration"), _fmt_exp(it.get("residual"))))
                lines.append(rule)
            else:
                summary = ("   Auxiliary %s solve: %s in %s iterations"
                           % (krylov.get("solver", "gmres").upper(),
                              krylov.get("outcome") or "recorded",
                              krylov.get("iteration_count",
                                         len(krylov_iterations))))
                if krylov.get("residual") is not None:
                    summary += " (residual %.2E)" % krylov["residual"]
                if position == 1:
                    lines.append("")
                lines.append(summary)
    outcome = _ZVECTOR_OUTCOMES.get(solve.get("outcome"),
                                    solve.get("outcome") or "recorded")
    summary = "   Z-vector %s" % outcome
    details = []
    if solve.get("pair_iterations"):
        details.append("%d pair-space Krylov iterations"
                       % solve["pair_iterations"])
    if solve.get("anderson_iterations"):
        details.append("%d fixed-point iterations"
                       % solve["anderson_iterations"])
    if solve.get("fallback_gmres_iterations"):
        details.append("%d fallback GMRES iterations"
                       % solve["fallback_gmres_iterations"])
    if solve.get("residual") is not None:
        details.append("residual %.2E" % solve["residual"])
    if details:
        summary += ": " + ", ".join(details)
    summary += _fmt_seconds(solve.get("elapsed"))
    lines.append(summary)
    return "\n".join(lines)


def _dedupe_iterations(iterations):
    """Keep the last record per iteration number (restart edges repeat one)."""
    merged = {}
    order = []
    for it in iterations:
        number = it.get("iteration")
        if number not in merged:
            order.append(number)
        merged[number] = it
    return [merged[number] for number in order]


def format_zvector_dense_block(solve):
    lines = ["   Z-vector dense direct solve"]
    if solve.get("dimension") is not None:
        lines.append("   dimension = %s" % solve["dimension"])
    for stage, elapsed in solve.get("stages", []):
        label = stage.replace("_", " ")
        lines.append("     %-24s %s" % (
            label, ("%.3f s" % elapsed) if elapsed is not None else ""))
    outcome = _ZVECTOR_OUTCOMES.get(solve.get("outcome"),
                                    solve.get("outcome") or "recorded")
    summary = "   Z-vector dense solve %s" % outcome
    if solve.get("residual") is not None:
        summary += ": residual %.2E" % solve["residual"]
    summary += _fmt_seconds(solve.get("elapsed"))
    lines.append(summary)
    return "\n".join(lines)


AU_TO_DEBYE = 2.541746451895

_SCC_MIXER_NAMES = {0: "linear", 1: "anderson", 2: "broyden", 3: "auto",
                    4: "diis", 5: "trust"}
_RESPONSE_SOLVER_NAMES = {0: "dense", 1: "davidson", 2: "auto"}


def format_resolved_options(options):
    """Grouped table of the EFFECTIVE openqp-dftb options after resolution.

    The values come from the native openqp_dftb_resolved_options record, so a
    named preset (e.g. dtcam-tb) is expanded to its actual published parameter
    vector without duplicating those numbers in Python.
    """
    if not options:
        return ""
    preset = str(options.get("preset", "none"))

    def number(key, digits=4):
        value = options.get(key)
        if value is None:
            return "?"
        return "%.*f" % (digits, float(value))

    def selector(key):
        value = options.get(key)
        if value is None:
            return "?"
        return "off" if float(value) < 0.0 else number(key)

    def yes_no(key):
        return "yes" if options.get(key) else "no"

    def response_lc(key):
        value = options.get(key)
        if value is None:
            return "?"
        if float(value) < 0.0:
            return "inherit reference"
        return number(key)

    if preset != "none":
        header = ("   Effective operator and numerical controls "
                  "(preset: %s)" % preset)
    else:
        header = ("   Effective operator and numerical controls "
                  "(no preset: explicit [dftb] keys / schema defaults)")
    lines = [header, ""]

    lines.append("   Reference operator")
    lines.append("     %-31s %s" % (
        "LC ground state:", "yes" if options.get("lc_ground_state") else "no"))
    # The LC kernel parameters below also drive the SF/MRSF response coupling
    # (gamma^lr) whenever the response LC inherits the reference, so they are
    # reported even for a non-LC ground state.
    lines.append("     %-31s %s kernel" % (
        "LC gamma:", "erf" if options.get("lc_gamma_erf") else "yukawa"))
    lines.append("     %-31s %s / %s / %s" % (
        "omega / cam_alpha / cam_beta:", number("omega"),
        number("cam_alpha"), number("cam_beta")))
    lines.append("     %-31s %s (W scale %s)" % (
        "Spin constants:",
        "built-in" if options.get("builtin_spin_constants")
        else "external spinw.txt", number("w_scale")))

    lines.append("   Response operator")
    lines.append("     %-31s %s / %s / %s" % (
        "omega / cam_alpha / cam_beta:", response_lc("response_omega"),
        response_lc("response_cam_alpha"), response_lc("response_cam_beta")))
    lines.append("     %-31s %s" % (
        "Response spin W scale:", response_lc("response_w_scale")))
    lines.append("     %-31s %s / %s / %s" % (
        "SPC coco / ovov / coov:", number("spc_coco"), number("spc_ovov"),
        number("spc_coov")))
    lines.append("     %-31s %s / %s / %s" % (
        "Onsite ss / sp / pp:", number("onsite_ss"), number("onsite_sp"),
        number("onsite_pp")))
    if options.get("onsite_exchange_scale"):
        lines.append("     %-31s %s" % (
            "Onsite exchange scale:", number("onsite_exchange_scale")))
    lines.append("     %-31s %s / %s" % (
        "Legacy c_mrsf / c_mrsf_oo:", selector("c_mrsf"),
        selector("c_mrsf_oo")))
    if options.get("response_global_hybrid") or options.get(
            "global_hybrid_exchange"):
        lines.append("     %-31s response %s / reference %s" % (
            "Global hybrid exchange:", yes_no("response_global_hybrid"),
            yes_no("global_hybrid_exchange")))

    lines.append("   Numerical controls")
    mixer = _SCC_MIXER_NAMES.get(options.get("scc_mixer"),
                                 str(options.get("scc_mixer")))
    max_step = options.get("scc_max_step")
    max_step_text = ("no clamp" if max_step is not None
                     and float(max_step) <= 0.0
                     else number("scc_max_step", 2))
    lines.append("     %-31s %s (mixing %s, history %s, max step %s)" % (
        "SCC mixer:", mixer, number("scc_mixing"),
        options.get("scc_history", "?"), max_step_text))
    lines.append("     %-31s %.1E / %s" % (
        "SCC tolerance / max iterations:",
        float(options.get("scc_tolerance") or 0.0),
        options.get("max_scc_iterations", "?")))
    solver = _RESPONSE_SOLVER_NAMES.get(options.get("response_solver"),
                                        str(options.get("response_solver")))
    lines.append("     %-31s %s (tolerance %.1E, max iterations %s, "
                 "subspace %s)" % (
                     "Response solver:", solver,
                     float(options.get("response_tolerance") or 0.0),
                     options.get("response_max_iterations", "?"),
                     options.get("response_max_subspace", "?")))
    lines.append("     %-31s %s" % (
        "Analytic Z-vector gradient:", yes_no("zvector")))
    return "\n".join(lines)


def format_step_timing(step_cpu, step_wall, total_cpu, total_wall):
    """The native-OpenQP style per-step + cumulative timing footer."""
    return (
        "   Step  CPU time (seconds): %10.3f   Wall time (seconds): %10.3f\n"
        "   Total CPU time (seconds): %10.3f   Wall time (seconds): %10.3f"
        % (step_cpu, step_wall, total_cpu, total_wall))


def format_energy_components(components, reference_label="SCC reference"):
    """OpenQP-style 'Energy components' block for the converged reference.

    Semantics (from the native emitter): electronic = h0 + scc + spin +
    exchange; total = electronic + repulsive; band = sum of occupied
    eigenvalues (informational).
    """
    if not components or components.get("total") is None:
        return ""

    def row(label, key, always=False):
        value = components.get(key)
        if value is None:
            return None
        if not always and abs(value) < 5.0e-13:
            return None
        return "   %-33s = %18.10f" % (label, value)

    rule = "   %s   ------------------" % (" " * 33)
    lines = ["   Energy components (%s)" % reference_label, ""]
    for entry in (
        row("Core Hamiltonian (H0) energy", "h0", always=True),
        row("SCC charge-fluctuation energy", "scc", always=True),
        row("Spin polarization energy", "spin"),
        row("Long-range exchange energy", "exchange"),
    ):
        if entry:
            lines.append(entry)
    lines.append(rule)
    lines.append(row("Electronic energy", "electronic", always=True))
    lines.append(row("Repulsive energy", "repulsive", always=True))
    lines.append(rule)
    lines.append(row("TOTAL energy", "total", always=True))
    if components.get("band") is not None:
        lines.append("")
        lines.append("   %-33s = %18.10f" % (
            "Band-structure energy (occ. sum)", components["band"]))
    return "\n".join(line for line in lines if line is not None)


def point_charge_dipole(charges, coords_bohr):
    """Dipole vector (a.u.) of net atomic point charges at coords (Bohr)."""
    dipole = [0.0, 0.0, 0.0]
    for charge, xyz in zip(charges, coords_bohr):
        for axis in range(3):
            dipole[axis] += float(charge) * float(xyz[axis])
    return dipole


def format_charge_block(symbols, charges, dipole_au, title):
    """Mulliken charge table + point-charge dipole, OpenQP style."""
    lines = ["   %s" % title]
    rule = "   " + "-" * 34
    lines.append(rule)
    lines.append("      #   Element        Charge")
    lines.append(rule)
    for index, (symbol, charge) in enumerate(zip(symbols, charges), start=1):
        lines.append("   %5d   %-7s  %12.6f" % (index, symbol, charge))
    lines.append(rule)
    lines.append("   %-7s  %12.6f" % ("Sum", sum(float(q) for q in charges)))
    magnitude = sum(component ** 2 for component in dipole_au) ** 0.5
    lines.append("")
    lines.append("   Dipole moment (point-charge model)")
    lines.append("      X = %10.6f   Y = %10.6f   Z = %10.6f  a.u." %
                 tuple(dipole_au))
    lines.append("      |mu| = %10.6f a.u. = %10.6f Debye" %
                 (magnitude, magnitude * AU_TO_DEBYE))
    return "\n".join(lines)


def format_other_lines(trace):
    lines = trace.get("other", []) + trace.get("warnings", [])
    if not lines:
        return ""
    return "\n".join(["   Additional native output:"]
                     + ["   %s" % line.strip() for line in lines])


def map_spectrum_label(label, state_prefix="S"):
    """Map a printed spectrum label to the public state label.

    Closed-shell TDDFTB rows already carry ``S<n>`` labels.  Spin-flip and
    MRSF rows are labeled ``root <n>``; within each spin manifold the public
    label is zero-based from root 1 (root 1 = S0/T0).
    """
    label = label.strip()
    if state_prefix and label.startswith("root"):
        number = _int(label.split()[-1])
        if number is None:
            return label
        return "%s%d" % (state_prefix, number - 1)
    return label


def spectrum_from_ground(spectrum, state_prefix="S"):
    """Rows out of the lowest state, keyed by public final-state label."""
    if not spectrum:
        return {}
    rows = {}
    ground_labels = {"S0", "root 1"}
    for row in spectrum.get("rows", []):
        if row["initial"].strip() in ground_labels:
            final = map_spectrum_label(row["final"], state_prefix)
            rows[final] = row
    return rows


def final_energy_header(summary):
    """One header line naming the extra Final-Energy columns."""
    if not summary:
        return ""
    return "  %12s  %12s  %12s" % ("VEE (eV)", "|mu| (a.u.)", "f")


def final_energy_annotation(summary, index):
    """Extra Final-Energy columns (VEE in eV, |mu|, f) for state row ``index``.

    Numbers only -- the column names come from :func:`final_energy_header`.
    ``index`` is the row number in ``mol.energies``: for spin-flip/MRSF runs
    row 0 is the internal reference and roots start at row 1; for closed-shell
    TDDFTB row 0 is S0 itself.
    """
    if not summary:
        return ""
    labels = summary.get("state_labels") or []
    energies = summary.get("state_energies") or []
    offset = 1 if summary.get("method") in ("mrsf", "sf") else 0
    position = index - offset
    if position < 0 or position >= min(len(labels), len(energies)):
        return ""
    vee = (energies[position] - energies[0]) * HARTREE_TO_EV
    text = "  %12.6f" % vee
    if position == 0:
        return text
    entry = (summary.get("oscillator") or {}).get(labels[position])
    if entry is not None:
        text += "  %12.4f  %12.4f" % (
            entry.get("dipole") or 0.0, entry.get("oscillator") or 0.0)
    return text


def final_energy_label(summary, index, default):
    """Public state label for Final-Energy row ``index`` (TDDFTB: S0..Sn)."""
    if not summary or summary.get("method") in ("mrsf", "sf"):
        return default
    labels = summary.get("state_labels") or []
    if 0 <= index < len(labels):
        return labels[index]
    return default


def _display_irrep(irrep):
    """Human display form of an irrep label ('ag' -> 'Ag'; 'mixed' as-is)."""
    text = str(irrep or "")
    if not text or text in {"mixed", "?"}:
        return text
    return text[:1].upper() + text[1:]


def format_mo_table(mo):
    """Molecular-orbital table: index, occupation, energies, irrep, markers."""
    energies = mo.get("energies") or []
    if not len(energies):
        return ""
    occupations = mo.get("occupations") or []
    labels = mo.get("labels") or []
    noca = int(mo.get("noca") or 0)
    nocb = int(mo.get("nocb") or 0)
    has_labels = bool(labels)
    lines = ["   Molecular orbitals (%s)"
             % mo.get("reference", "SCC reference")]
    rule = "   " + "-" * (60 if has_labels else 50)
    lines.append(rule)
    header = "      MO   Occ      Energy (Ha)     Energy (eV)"
    if has_labels:
        header += "     Irrep"
    lines.append(header)
    lines.append(rule)
    open_shell = noca > nocb
    for index, energy in enumerate(energies, start=1):
        occupation = (occupations[index - 1]
                      if index - 1 < len(occupations) else 0)
        row = "   %5d   %3d   %14.8f   %12.4f" % (
            index, int(occupation), energy, energy * HARTREE_TO_EV)
        if has_labels:
            label = labels[index - 1] if index - 1 < len(labels) else ""
            row += "   %7s" % _display_irrep(label)
        if open_shell:
            if nocb < index <= noca:
                row += "   (SOMO)"
            elif index == nocb and nocb > 0:
                row += "   (HOMO)"
            elif index == noca + 1:
                row += "   (LUMO)"
        else:
            if index == noca and noca > 0:
                row += "   (HOMO)"
            elif index == noca + 1:
                row += "   (LUMO)"
        lines.append(row)
    lines.append(rule)
    return "\n".join(lines)


def spin_flip_pair(index, noca, nocb):
    """Map a 1-based spin-flip amplitude index to (occ, vir) MO numbers.

    The exported SF/MRSF response vectors are ``amplitudes(noca, n_virt_beta)``
    flattened column-major: index = (vir - nocb - 1) * noca + occ, with occ an
    alpha-occupied MO (1..noca) and vir a beta-virtual MO (nocb+1..nbf).
    """
    position = int(index) - 1
    return position % int(noca) + 1, int(nocb) + position // int(noca) + 1


def format_state_configurations(summary):
    """AE-style per-state dominant configurations (coeff, OCC -> VIR)."""
    configurations = summary.get("configurations") or {}
    if not any(configurations.values()):
        return ""
    labels = summary.get("state_labels") or []
    energies = summary.get("state_energies") or []
    spins = summary.get("spin_squares") or []
    state_irreps = summary.get("state_irreps") or {}
    base = energies[0] if energies else None
    lines = ["     Spin-adapted spin-flip excitations", ""]
    for position, label in enumerate(labels):
        entries = configurations.get(label)
        if not entries:
            continue
        header = "   State %-6s" % label
        if base is not None and position < len(energies):
            header += "  VEE = %10.6f eV" % (
                (energies[position] - base) * HARTREE_TO_EV)
        if position < len(spins) and spins[position] is not None:
            header += "     <S^2> = %7.4f" % spins[position]
        if state_irreps.get(label):
            header += "     Sym = %s" % _display_irrep(state_irreps[label])
        lines.append(header)
        lines.append("         #      Coeff        OCC       VIR")
        lines.append("       ----   ---------     -----     -----")
        for entry in entries:
            lines.append("      %5d   %9.6f     %5d  ->  %4d" % (
                entry["index"], entry["coeff"], entry["occ"], entry["vir"]))
        lines.append("")
    return "\n".join(lines).rstrip()


def format_excited_state_summary(summary):
    """The OpenQP-style excited-state summary + state-to-state tables.

    ``summary`` is the dict the adapter stores on the molecule:
    ``reference_energy``, ``reference_label``, ``state_labels`` (public),
    ``state_energies`` (total Hartree per root), ``spin_squares``,
    ``oscillator`` ({public label: spectrum row}), ``transitions``
    (spectrum rows with public labels), ``approximation``.
    """
    labels = summary.get("state_labels") or []
    energies = summary.get("state_energies") or []
    spins = summary.get("spin_squares") or []
    oscillator = summary.get("oscillator") or {}

    def dipole_columns(entry):
        """X, Y, Z, Abs, f columns for one ground->state spectrum row."""
        if entry is None:
            return "%11s%11s%11s%10s%11s" % ("", "", "", "n/a", "n/a")
        xyz = entry.get("dipole_xyz")
        text = ""
        if xyz is not None:
            text += "%11.4f%11.4f%11.4f" % tuple(xyz)
        else:
            text += "%11s%11s%11s" % ("", "", "")
        text += "%10.4f" % (entry.get("dipole") or 0.0)
        text += "%11.4f" % (entry.get("oscillator") or 0.0)
        return text

    lines = []
    configurations_text = format_state_configurations(summary)
    if configurations_text:
        lines.extend([configurations_text, ""])
    state_irreps = summary.get("state_irreps") or {}
    lines.extend(["     Summary table", ""])
    lines.append("   State   Sym       Energy        Excitation   <S^2>"
                 "          Transition dipole (a.u.)        Oscillator")
    lines.append("                     Hartree           eV               "
                 "      X          Y          Z       Abs.    strength")
    rule = "   " + "-" * 108
    lines.append(rule)
    base = energies[0] if len(energies) else None
    for position, label in enumerate(labels):
        energy = energies[position] if position < len(energies) else None
        row = "   %-6s%-8s" % (
            label, _display_irrep(state_irreps.get(label, "")))
        row += ("%18.10f" % energy) if energy is not None else " " * 18
        if energy is not None and base is not None:
            row += "%12.6f" % ((energy - base) * HARTREE_TO_EV)
        else:
            row += " " * 12
        spin = spins[position] if position < len(spins) else None
        row += ("%8.3f" % spin) if spin is not None else " " * 8
        if label == labels[0]:
            row += "%11s" % "--"
        else:
            row += dipole_columns(oscillator.get(label))
        lines.append(row.rstrip())
    lines.append(rule)
    reference_energy = summary.get("reference_energy")
    if reference_energy is not None:
        lines.append("   %-22s %18.10f" % (
            summary.get("reference_label", "Reference (internal)"),
            reference_energy))
    transitions = summary.get("transitions") or []
    if transitions:
        lines.append("")
        lines.append("     State-to-state transitions")
        if summary.get("approximation"):
            lines.append("     (%s)" % summary["approximation"])
        lines.append("")
        lines.append("    Initial    Final       dE (eV)          X"
                     "          Y          Z       Abs.          f")
        lines.append("   " + "-" * 90)
        for row in transitions:
            text = "   %7s   %6s   %12.6f" % (
                row["initial"], row["final"], row["de_ev"] or 0.0)
            xyz = row.get("dipole_xyz")
            if xyz is not None:
                text += "%11.4f%11.4f%11.4f" % tuple(xyz)
            else:
                text += "%11s%11s%11s" % ("", "", "")
            text += "%11.4f%11.6f" % (row["dipole"] or 0.0,
                                      row["oscillator"] or 0.0)
            lines.append(text)
    return "\n".join(lines)
