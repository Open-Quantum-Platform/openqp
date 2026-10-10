"""Local execution hardware and conservative memory admission checks."""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

GiB = 1024 ** 3


def _sysctl(name: str) -> int | None:
    try:
        result = subprocess.run(
            ["/usr/sbin/sysctl" if Path("/usr/sbin/sysctl").is_file() else "sysctl", "-n", name],
            capture_output=True, text=True,
            timeout=1, check=False,
        )
        return int(result.stdout.strip()) if result.returncode == 0 else None
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


def _memory_bytes() -> tuple[int | None, int | None]:
    """Return physical and currently reusable memory, without extra packages."""
    if sys.platform == "darwin":
        total = _sysctl("hw.memsize")
        page_size = _sysctl("hw.pagesize") or 4096
        try:
            stat = subprocess.run(["/usr/bin/vm_stat"], capture_output=True, text=True,
                                  timeout=1, check=False).stdout
            pages = sum(
                int(match.group(1))
                for label in ("Pages free", "Pages inactive", "Pages speculative")
                for match in [re.search(rf"^{label}:\s+(\d+)", stat, re.MULTILINE)]
                if match
            )
            return total, pages * page_size
        except (OSError, subprocess.SubprocessError):
            return total, None
    if sys.platform.startswith("linux"):
        values: dict[str, int] = {}
        try:
            with open("/proc/meminfo", encoding="utf-8") as memory_info:
                for line in memory_info:
                    key, value = line.split(":", 1)
                    values[key] = int(value.strip().split()[0]) * 1024
            return values.get("MemTotal"), values.get("MemAvailable")
        except (OSError, ValueError):
            return None, None
    return None, None


def snapshot() -> dict:
    """Compact host facts for the Execution panel."""
    physical = _sysctl("hw.physicalcpu") if sys.platform == "darwin" else None
    logical = _sysctl("hw.logicalcpu") if sys.platform == "darwin" else os.cpu_count()
    total, available = _memory_bytes()
    return {
        "platform": sys.platform,
        "physical_cores": physical or logical or 1,
        "logical_cores": logical or physical or 1,
        "memory_total_bytes": total,
        "memory_available_bytes": available,
    }


def _atom_count(input_text: str) -> int:
    geometry = re.search(r'geom\s*=\s*"""(.*?)"""', input_text,
                         re.IGNORECASE | re.DOTALL)
    block = geometry.group(1) if geometry else input_text
    return sum(
        bool(re.match(r"\s*(?:[A-Z][a-z]?|\d+)\s+[-+0-9.]", line))
        for line in block.splitlines()
    )


def estimate_memory(input_text: str, threads: int) -> int:
    """Return a deliberately cautious RAM estimate in bytes.

    OpenQP's exact allocation depends on the basis, active space and integral
    screening, which are only known after setup.  This preflight estimate is
    an admission guard, not a replacement for the engine's own allocation.
    """
    atoms = max(1, _atom_count(input_text))
    route = input_text.lower()
    base = 512 * 1024**2 + atoms * 16 * 1024**2
    if any(method in route for method in ("ccsd", "caspt2", "nevpt2", "mrmp2", "mcqdpt2")):
        base += 2 * GiB + atoms**3 * 2 * 1024**2
    elif any(method in route for method in ("fci", "casci", "casscf")):
        base += 2 * GiB + atoms**3 * 1024**2
    elif any(method in route for method in ("tddft", "tda", "mrsf", "umrsf", "sf(")):
        base += 768 * 1024**2 + atoms**2 * 4 * 1024**2
    elif "mp2" in route:
        base += GiB + atoms**3 * 1024**2
    return base + max(1, threads) * 128 * 1024**2


def admission(input_text: str, threads: int) -> dict:
    facts = snapshot()
    estimate = estimate_memory(input_text, threads)
    available = facts["memory_available_bytes"]
    # Keep a 20% reserve for macOS, the GUI, and the operating system.
    permitted = available is None or estimate <= available * 0.8
    return {
        **facts,
        "threads": threads,
        "estimated_memory_bytes": estimate,
        "permitted": permitted,
        "reason": None if permitted else (
            "Estimated calculation memory exceeds the currently available RAM "
            "after reserving memory for macOS and OQP Studio."
        ),
    }
