"""Execution adapters that run an OpenQP input file and stream its log.

Adapters share one contract (`Runner`): given a job directory containing
an .oqp or .inp input, launch OpenQP, append stdout/stderr to `job.log`, and
report the exit code. The GUI never cares which adapter ran the job.
"""

from .base import Runner, RunnerUnavailable
from .bundled import BundledRunner
from .local import LocalRunner
from .wsl import WslRunner

_REGISTRY = {
    "local": LocalRunner,
    "bundled": BundledRunner,
    "wsl": WslRunner,
}


def get_runner(name: str) -> Runner:
    try:
        return _REGISTRY[name]()
    except KeyError:
        raise ValueError(f"unknown runner '{name}' (available: {sorted(_REGISTRY)})") from None


def available_runners() -> dict[str, bool]:
    """Map of runner name -> whether it can execute on this machine."""
    return {name: cls().is_available() for name, cls in _REGISTRY.items()}


__all__ = ["Runner", "RunnerUnavailable", "available_runners", "get_runner"]
