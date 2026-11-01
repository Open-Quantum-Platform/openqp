#!/usr/bin/env python3
"""CCSD and CR-CC(2,3) with the external CCpy library on OpenQP RHF orbitals.

CCpy (https://github.com/piecuch-group/ccpy, GPL-3.0) is not part of OpenQP.
On first use ``oqp.interop.ccpy`` downloads and builds a pinned CCpy commit
into the running Python environment (needs Python >= 3.10, git, a Fortran
compiler and OpenBLAS). Set ``OPENQP_CCPY_AUTO_INSTALL=0`` to forbid that.

This example is run-only: it needs the optional CCpy backend. The same
workflow is checked against PySCF -> CCpy reference energies in
``tests/test_ccpy_interop.py``.

Usage
-----
    python ccpy_h2o.py h2o_rhf_sto-3g_ccpy.inp
"""

import sys

from oqp.interop.ccpy import driver_from_openqp
from oqp.pyoqp import Runner


def main():
    if len(sys.argv) < 2:
        sys.exit(f"usage: {sys.argv[0]} <input.inp>")
    input_file = sys.argv[1]
    project = input_file.rsplit(".", 1)[0]
    runner = Runner(project=project, input_file=input_file, log=f"{project}.log")
    runner.run()                                   # method=hf: converged RHF

    driver = driver_from_openqp(runner.mol, nfrozen=1)
    driver.run_cc(method="ccsd")
    e_ccsd = driver.system.reference_energy + driver.correlation_energy

    driver.run_hbar(method="ccsd")
    driver.run_leftcc(method="left_ccsd")
    driver.run_ccp3(method="crcc23")
    e_crcc23 = e_ccsd + driver.deltap3[0]["D"]

    print(f"E(RHF)          = {driver.system.reference_energy:.10f}")
    print(f"E(CCSD)         = {e_ccsd:.10f}")
    print(f"E(CR-CC(2,3)_D) = {e_crcc23:.10f}")


if __name__ == "__main__":
    main()
