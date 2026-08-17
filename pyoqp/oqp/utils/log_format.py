"""Shared, compatibility-preserving formatting for the OpenQP text log.

The text log predates the Python driver and is consumed by user scripts.  The
``PyOQP`` field prefix and native solver markers are therefore stable public
markers.  This module standardizes the surrounding section grammar, alignment,
units, and precision without renaming those markers.
"""

from __future__ import annotations


LOG_PREFIX = "PyOQP"
LOG_RULE_WIDTH = 72
LOG_FIELD_WIDTH = 28
ENERGY_DECIMALS = 10

RUN = "RUN"
INPUT_REFERENCE = "INPUT AND REFERENCE"
CONVERGENCE = "CONVERGENCE AND ITERATIONS"
ENERGIES = "ENERGIES AND STATES"
GRADIENTS_PROPERTIES = "GRADIENTS AND PROPERTIES"
TERMINATION = "TIMING AND TERMINATION"
PROGRESS = "CALCULATION PROGRESS"

LOG_SECTION_ORDER = (
    RUN,
    INPUT_REFERENCE,
    PROGRESS,
    CONVERGENCE,
    ENERGIES,
    GRADIENTS_PROPERTIES,
    TERMINATION,
)


# All run types accepted by the Runner or retained as a legacy diagnostic are
# listed explicitly.  Keeping this map complete makes a newly added route fail
# a source-only test until its log structure is considered.
RUNTYPE_LOG_CATEGORIES = {
    "energy": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES, TERMINATION),
    "ekt": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES, TERMINATION),
    "soc": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
            GRADIENTS_PROPERTIES, TERMINATION),
    "grad": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
             GRADIENTS_PROPERTIES, TERMINATION),
    "hess": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
             GRADIENTS_PROPERTIES, TERMINATION),
    "thermo": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
               GRADIENTS_PROPERTIES, TERMINATION),
    "nac": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
            GRADIENTS_PROPERTIES, TERMINATION),
    "nacme": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
              GRADIENTS_PROPERTIES, TERMINATION),
    "bp": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
           GRADIENTS_PROPERTIES, TERMINATION),
    "prop": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
             GRADIENTS_PROPERTIES, TERMINATION),
    "data": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
             GRADIENTS_PROPERTIES, TERMINATION),
    "optimize": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
                 GRADIENTS_PROPERTIES, TERMINATION),
    "meci": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
             GRADIENTS_PROPERTIES, TERMINATION),
    "mecp": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
             GRADIENTS_PROPERTIES, TERMINATION),
    "tci": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
            GRADIENTS_PROPERTIES, TERMINATION),
    "mep": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
            GRADIENTS_PROPERTIES, TERMINATION),
    "ts": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
           GRADIENTS_PROPERTIES, TERMINATION),
    "irc": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
            GRADIENTS_PROPERTIES, TERMINATION),
    "neb": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
            GRADIENTS_PROPERTIES, TERMINATION),
    "namd": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
             GRADIENTS_PROPERTIES, TERMINATION),
    "md": (INPUT_REFERENCE, PROGRESS, CONVERGENCE, ENERGIES,
           GRADIENTS_PROPERTIES, TERMINATION),
}


_SECTION_CATEGORIES = {
    "start": RUN,
    "calculation": INPUT_REFERENCE,
    "input": INPUT_REFERENCE,
    "guess": INPUT_REFERENCE,
    "basis_overlap": INPUT_REFERENCE,
    "symmetry": INPUT_REFERENCE,
    "dftb": INPUT_REFERENCE,
    "method": INPUT_REFERENCE,
    "scf": CONVERGENCE,
    "tdhf": CONVERGENCE,
    "fci": ENERGIES,
    "correlation": CONVERGENCE,
    "casscf_macroiteration": CONVERGENCE,
    "opt": CONVERGENCE,
    "QM/MM": CONVERGENCE,
    "cons_sphere": CONVERGENCE,
    "baeka": CONVERGENCE,
    "penalty": CONVERGENCE,
    "auglag": CONVERGENCE,
    "hybrid": CONVERGENCE,
    "ubp": CONVERGENCE,
    "mecp": CONVERGENCE,
    "mep": CONVERGENCE,
    "num_hess": CONVERGENCE,
    "hess_worker": CONVERGENCE,
    "num_nacv": CONVERGENCE,
    "nacv_worker": CONVERGENCE,
    "text": CONVERGENCE,
    "dftb_runtime": CONVERGENCE,
    "energy": ENERGIES,
    "dftd": INPUT_REFERENCE,
    "dftb_state_summary": ENERGIES,
    "grad": GRADIENTS_PROPERTIES,
    "nacme": GRADIENTS_PROPERTIES,
    "nacm": GRADIENTS_PROPERTIES,
    "nacv": GRADIENTS_PROPERTIES,
    "dcv": GRADIENTS_PROPERTIES,
    "bp": GRADIENTS_PROPERTIES,
    "read_hess": GRADIENTS_PROPERTIES,
    "freq": GRADIENTS_PROPERTIES,
    "freq_modes": GRADIENTS_PROPERTIES,
    "thermo": GRADIENTS_PROPERTIES,
    "end": TERMINATION,
}


def section_category(section=None):
    """Return the standard category for one legacy ``dump_log`` section."""

    return _SECTION_CATEGORIES.get(section, PROGRESS)


def format_log_section(title, category=PROGRESS):
    """Return one uniform section heading while retaining the legacy title."""

    rule = "=" * LOG_RULE_WIDTH
    lines = ["", "   " + rule, f"   {LOG_PREFIX} LOG | {category}"]
    if title:
        lines.append(f"   {str(title).strip()}")
    lines.append("   " + rule)
    return "\n".join(lines) + "\n"


def format_log_fields(rows, *, prefix=LOG_PREFIX, label_width=LOG_FIELD_WIDTH,
                      indent=3):
    """Format ordered key/value rows with the established ``PyOQP`` marker."""

    lead = " " * indent
    marker = (str(prefix) + " ") if prefix else ""
    return "\n".join(
        f"{lead}{marker}{str(label) + ':':<{label_width}} {value}"
        for label, value in rows
    )


def format_module_banner(title, info=""):
    """Render the Python method banner in the native 40-column grammar."""

    bar = " " * 20 + "+" * 40
    lines = ["", bar, " " * 23 + "MODULE: " + str(title)]
    if info:
        lines.append(" " * 23 + str(info))
    lines.append(bar)
    return "\n".join(lines) + "\n"


def format_energy(value, width=16):
    """Format a total energy consistently with the native SCF precision."""

    return f"{float(value):<{width}.{ENERGY_DECIMALS}f}"


def format_unit(label, unit):
    """Return an additive, parser-friendly unit declaration."""

    return format_log_fields(((f"{label} unit", unit),))
