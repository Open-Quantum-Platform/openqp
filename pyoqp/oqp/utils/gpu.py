"""Runtime helpers for the experimental OpenQP GPU backend.

The first GPU project targets MRSF/UMRSF TDHF Davidson METC contractions.
This module intentionally stays pure Python so input parsing and tests do not
require CUDA libraries or a built OpenQP shared library.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class GpuConfig:
    """Normalized GPU runtime configuration."""

    enabled: bool = False
    backend: str = "cuda"
    target: str = "metc"
    device: int = 0
    precision: str = "float64"
    fallback: str = "cpu"

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> "GpuConfig":
        """Create a normalized GPU config from an OpenQP config dictionary."""

        gpu = config.get("gpu", {}) or {}
        return cls(
            enabled=bool(gpu.get("enabled", False)),
            backend=str(gpu.get("backend", "cuda")).lower(),
            target=str(gpu.get("target", "metc")).lower(),
            device=int(gpu.get("device", 0)),
            precision=str(gpu.get("precision", "float64")).lower(),
            fallback=str(gpu.get("fallback", "cpu")).lower(),
        )

    def supports_metc(self, config: dict[str, Any]) -> bool:
        """Return True when this config targets the current GPU METC scope."""

        input_section = config.get("input", {}) or {}
        tdhf_section = config.get("tdhf", {}) or {}
        method = str(input_section.get("method", "hf")).lower()
        td_type = str(tdhf_section.get("type", "rpa")).lower()
        return (
            self.enabled
            and self.backend == "cuda"
            and self.target == "metc"
            and self.precision == "float64"
            and method == "tdhf"
            and td_type in {"mrsf", "umrsf"}
        )
