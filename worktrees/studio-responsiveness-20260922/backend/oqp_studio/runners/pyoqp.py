from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

from .base import Runner


class PyOqpRunner(Runner):
    """Runs OpenQP through the pyoqp Python API in a worker subprocess.

    Uses `oqp.pyoqp.Runner` (the same layer the `openqp` CLI uses) and dumps
    `runner.results()` to `results.json` in the job directory, so the GUI
    gets structured energies/geometry without log scraping. A subprocess —
    not an in-server import — keeps a crashed native kernel or an MPI init
    conflict from taking down the Studio backend.
    """

    name = "pyoqp"

    def is_available(self) -> bool:
        return importlib.util.find_spec("oqp") is not None

    def build_command(self, input_file: Path) -> list[str]:
        return [
            sys.executable,
            "-u",
            "-m",
            "oqp_studio.pyoqp_worker",
            str(input_file.parent),
        ]
