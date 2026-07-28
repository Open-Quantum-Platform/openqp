"""Fallback entry point for OpenQP installations built without AFQMC."""

import sys


def installed_main():
    """Explain how to rebuild when the optional native AFQMC target is absent."""
    print(
        "OpenQP was installed without the optional AFQMC source checkout.\n"
        "Place AFQMC at openqp/external/afqmc or beside openqp as afqmc, then "
        "run `python -m pip install .` again.\n"
        "For another location, set "
        "-C cmake.define.OPENQP_AFQMC_SOURCE_DIR=/absolute/path/to/afqmc.",
        file=sys.stderr,
    )
    return 2
