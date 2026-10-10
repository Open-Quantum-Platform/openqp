from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

from .base import Runner


def _windows_path_to_wsl(path: Path) -> str:
    """C:\\Users\\x\\job -> /mnt/c/Users/x/job."""
    drive = path.drive.rstrip(":").lower()
    rest = str(path).replace("\\", "/").split(":", 1)[1]
    return f"/mnt/{drive}{rest}"


class WslRunner(Runner):
    """Runs OpenQP inside WSL from a Windows host.

    Bridges the GUI on native Windows to an OpenQP built in a WSL distro,
    until the native Windows build of the OpenQP core is ready. Job
    directories live on the Windows filesystem and are reached from WSL
    via /mnt/<drive>/..., so results are immediately visible to the GUI.
    """

    name = "wsl"

    @property
    def command(self) -> str:
        return os.environ.get("OQP_WSL_COMMAND", "openqp")

    def is_available(self) -> bool:
        if sys.platform != "win32" or shutil.which("wsl") is None:
            return False
        probe = subprocess.run(
            ["wsl", "--", "which", self.command],
            capture_output=True,
            timeout=15,
            check=False,
        )
        return probe.returncode == 0

    def build_command(self, input_file: Path) -> list[str]:
        return ["wsl", "--", self.command, _windows_path_to_wsl(input_file)]
