#!/usr/bin/env python3
"""Export an OpenQP mean-field calculation as a FCIDUMP for quantum computing.

FCIDUMP is the standard hand-off point from a classical quantum-chemistry
mean field to a quantum-computing electronic-structure workflow. The file this
script writes can be loaded by, e.g.:

* Qiskit Nature  -- ``qiskit_nature.second_q.formats.fcidump.FCIDump.from_file``
* OpenFermion    -- via PySCF's ``MolecularData`` / FCIDUMP loaders
* Block2 / DMRG, Dice/SHCI, and most active-space front ends

The two-electron integrals are computed natively by OpenQP (``oqp.int2e`` ->
``OQP::ERI_AO``) and transformed to the MO basis, so no external integral
source is needed.

Usage
-----
    openqp h2.inp                       # run an HF/DFT single point first
    python export_fcidump.py h2.inp h2.FCIDUMP
"""

import os
import sys

from oqp.pyoqp import Runner
from oqp.quantum import from_openqp


def main():
    if len(sys.argv) < 3:
        sys.exit(f"usage: {sys.argv[0]} <input.inp> <out.FCIDUMP>")
    input_file, out = sys.argv[1], sys.argv[2]

    project = input_file.rsplit(".", 1)[0]
    log = f"{project}.log"
    runner = Runner(project=project, input_file=input_file, log=log)
    runner.run()
    mol = runner.mol

    # Builds h_pq from OQP::Hcore + MOs and (pq|rs) from OQP::ERI_AO.
    ham = from_openqp(mol)
    print(f"norb = {ham.n_orbitals}, nelec = {ham.n_electrons}, "
          f"ms2 = {ham.ms2}, E_core = {ham.core_energy:.10f}")

    ham.to_fcidump(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
