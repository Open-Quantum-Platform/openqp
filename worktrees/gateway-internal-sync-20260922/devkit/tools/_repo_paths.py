"""Canonical paths for tools integrated under ``<openqp>/devkit``."""

from pathlib import Path


DEVKIT_ROOT = Path(__file__).resolve().parents[1]
ENGINE_ROOT = DEVKIT_ROOT.parent
EXAMPLES_ROOT = ENGINE_ROOT / "examples"
SCHEMA_FILE = ENGINE_ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
