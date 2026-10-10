from __future__ import annotations

import os
from pathlib import Path

from .. import environment
from .base import Runner


class LocalRunner(Runner):
    """Runs the OpenQP command available on this computer.

    This runner deliberately uses PATH only, so selecting it never silently
    switches to the engine shipped with or downloaded by Studio.
    """

    name = "local"

    @property
    def command(self) -> str:
        command = os.environ.get("OQP_COMMAND", "openqp")
        return environment.locate(command) or command

    def is_available(self) -> bool:
        return Path(self.command).is_file()

    def build_command(self, input_file: Path) -> list[str]:
        return [self.command, str(input_file)]
