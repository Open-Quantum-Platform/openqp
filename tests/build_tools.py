"""Locate the build tools that configuration tests run directly."""
import shutil
from pathlib import Path


def find_cmake():
    """Return a cmake executable, or raise AssertionError naming the requirement.

    Homebrew installs under /opt/homebrew on Apple Silicon and /usr/local on
    Intel Macs; neither is on the PATH of a non-login shell (ssh, launchd).
    """
    for candidate in (shutil.which("cmake"), "/opt/homebrew/bin/cmake", "/usr/local/bin/cmake"):
        if candidate and Path(candidate).is_file():
            return candidate
    raise AssertionError("this test runs CMake scripts and needs a cmake executable; none found")
