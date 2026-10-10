"""Stable TDHF/Davidson timing labels for GPU benchmark branches.

This module is intentionally dependency-light so GPU/METC branches can import it
from source-level tests and analysis scripts without requiring a built OpenQP.
Runtime Fortran/CUDA timers should emit the same ``OQP_TIMER`` line format so
benchmark collectors can compare CPU baselines, METC kernels, persistent-buffer
branches, and later CUDA variants with one parser.

``OQP_TIMER`` lines are written to the OpenQP run log (``<input>.log`` /
``infos%log_filename``), not necessarily to process stdout. When collecting
under a scheduler, point the parser at the ``.log`` file. If timer lines are
instead captured from stdout (e.g. a tee'd log), launchers may prepend a
per-line prefix such as a Slurm ``--label`` rank tag (``0: OQP_TIMER ...``);
the parser tolerates such leading prefixes and locates the ``OQP_TIMER`` token.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

_TIMER_PREFIX = "OQP_TIMER"


@dataclass(frozen=True)
class TDHFTimerLabel:
    label: str
    description: str
    private_gpu_group: bool = True


def davidson_timer_manifest() -> tuple[TDHFTimerLabel, ...]:
    """Return stable labels for TDHF/MRSF Davidson benchmark timings."""

    return (
        TDHFTimerLabel("tdhf.response.total", "Total TDHF/SF/MRSF response wall time."),
        TDHFTimerLabel("tdhf.davidson.total", "Total Davidson solver wall time."),
        TDHFTimerLabel("tdhf.davidson.sigma_build", "Davidson sigma-vector build wall time."),
        TDHFTimerLabel("tdhf.davidson.metc_contract", "METC/ERI contraction contribution to sigma build."),
        TDHFTimerLabel("tdhf.davidson.eri_buffer", "ERI buffer setup/reuse/transfer preparation time."),
        TDHFTimerLabel("tdhf.davidson.orthonormalize", "Davidson subspace orthonormalization time."),
        TDHFTimerLabel("tdhf.davidson.residual_check", "Residual/error check and convergence bookkeeping time."),
    )


def _validate_label(label: str) -> None:
    labels = {timer.label for timer in davidson_timer_manifest()}
    if label not in labels:
        raise ValueError(f"Unknown TDHF/Davidson timer label: {label}")


def format_timer_line(
    label: str,
    elapsed_seconds: float,
    metadata: Mapping[str, object] | None = None,
) -> str:
    """Format a machine-parseable OpenQP timer line.

    Values are intentionally whitespace-free ``key=value`` tokens because they
    can be parsed robustly from normal OpenQP logs and Slurm stdout captures.
    """

    _validate_label(label)
    fields = [
        _TIMER_PREFIX,
        f"label={label}",
        f"seconds={elapsed_seconds:.6f}",
    ]
    for key in sorted(metadata or {}):
        value = str((metadata or {})[key])
        if any(ch.isspace() for ch in key) or any(ch.isspace() for ch in value):
            raise ValueError("Timer metadata keys and values must not contain whitespace")
        fields.append(f"{key}={value}")
    return " ".join(fields)


def parse_timer_line(line: str) -> dict[str, object]:
    """Parse an ``OQP_TIMER`` line produced by :func:`format_timer_line`.

    Any tokens before the ``OQP_TIMER`` marker are treated as a harmless
    launcher prefix (e.g. a Slurm ``--label`` rank tag ``0: OQP_TIMER ...``)
    and ignored, so timer lines parse from raw ``.log`` files as well as from
    scheduler stdout captures.
    """

    tokens = line.strip().split()
    try:
        start = tokens.index(_TIMER_PREFIX)
    except ValueError:
        raise ValueError("Not an OQP_TIMER line")
    parsed: dict[str, object] = {}
    for token in tokens[start + 1:]:
        if "=" not in token:
            raise ValueError(f"Malformed timer token: {token}")
        key, value = token.split("=", 1)
        parsed[key] = value
    if "label" not in parsed or "seconds" not in parsed:
        raise ValueError("Timer line requires label and seconds")
    _validate_label(str(parsed["label"]))
    seconds = parsed["seconds"]
    parsed["seconds"] = float(str(seconds))
    return parsed


def parse_timer_lines(log_text: str) -> list[dict[str, object]]:
    """Extract all machine-parseable ``OQP_TIMER`` records from log text.

    A line qualifies when ``OQP_TIMER`` appears as a whitespace-delimited
    token, allowing a leading launcher prefix while avoiding false matches on
    prose that merely mentions the marker as a substring.
    """

    records: list[dict[str, object]] = []
    for line in log_text.splitlines():
        if _TIMER_PREFIX in line.split():
            records.append(parse_timer_line(line))
    return records


def summarize_timer_records(records: list[dict[str, object]]) -> dict[str, dict[str, float | int]]:
    """Group parsed timer records by label for compact benchmark reports.

    TODO(benchmark): for cross-backend tables, also group by a backend/kernel
    metadata key (e.g. ``kernel=cpu_baseline`` vs ``kernel=cuda``) instead of
    collapsing all records under one label, and report ``seconds_min`` /
    ``seconds_max`` alongside the mean. Out of scope for this observability
    branch; the parser already preserves metadata tokens to enable it.
    """

    summary: dict[str, dict[str, float | int]] = {}
    for record in records:
        label = str(record["label"])
        _validate_label(label)
        bucket = summary.setdefault(label, {"count": 0, "seconds_total": 0.0, "seconds_mean": 0.0})
        bucket["count"] = int(bucket["count"]) + 1
        bucket["seconds_total"] = float(bucket["seconds_total"]) + float(str(record["seconds"]))
    for bucket in summary.values():
        bucket["seconds_mean"] = float(bucket["seconds_total"]) / int(bucket["count"])
    return summary


def format_timer_summary_csv(summary: Mapping[str, Mapping[str, float | int]]) -> str:
    """Format timer-summary rows for data snapshots and manuscript tables."""

    lines = ["label,count,seconds_total,seconds_mean"]
    for timer in davidson_timer_manifest():
        if timer.label not in summary:
            continue
        row = summary[timer.label]
        lines.append(
            f"{timer.label},{int(row['count'])},{float(row['seconds_total']):.6f},{float(row['seconds_mean']):.6f}"
        )
    return "\n".join(lines) + "\n"
