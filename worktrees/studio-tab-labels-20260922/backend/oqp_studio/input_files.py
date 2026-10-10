"""Input and calculation-log names within one job directory."""

from __future__ import annotations

from pathlib import Path

SUPPORTED_INPUT_SUFFIXES = frozenset({".inp", ".oqp"})


def find_input_file(job_dir: Path) -> Path:
    """Return the input file prepared for a job, including a custom name."""
    for name in ("input.oqp", "input.inp"):
        candidate = job_dir / name
        if candidate.is_file():
            return candidate
    candidates = sorted(
        path for path in job_dir.iterdir()
        if path.is_file() and path.suffix.lower() in SUPPORTED_INPUT_SUFFIXES
    )
    return candidates[0] if candidates else job_dir / "input.oqp"


def calculation_log(input_file: Path) -> Path:
    """OpenQP writes its record next to an input with the same stem."""
    return input_file.with_suffix(".log")
