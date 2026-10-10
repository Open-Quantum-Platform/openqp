from __future__ import annotations

from pathlib import Path

from .. import engine
from .base import Runner


class BundledRunner(Runner):
    """Runs the OpenQP engine packaged with or downloaded by Studio."""

    name = "bundled"

    @property
    def command(self) -> str:
        return engine.bundled_or_downloaded() or engine.EXECUTABLE

    def is_available(self) -> bool:
        return Path(self.command).is_file()

    def build_command(self, input_file: Path) -> list[str]:
        return [self.command, str(input_file)]
