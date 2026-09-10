"""Retired standalone PyOQP packaging entry point.

OpenQP 1.3.1 is built only from the repository root, where ``pyproject.toml``
owns the license expression, bundled license files, native library, and Python
wrapper as one distribution.  Keeping a second setuptools build below
``pyoqp/`` would permit an incomplete artifact without those legal files.
"""

OPENQP_VERSION = "1.3.1"

raise SystemExit(
    "pyoqp/setup.py is retired; build OpenQP from the repository root with "
    "`python -m build` or `python -m pip install .`"
)
