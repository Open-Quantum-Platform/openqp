#!/usr/bin/env python3
"""Generate the radial SAP (superposition of atomic potentials) data file for OpenQP.

This is a one-time, offline data-generation tool. It extracts Susi Lehtola's
tabulated effective atomic charges Z_eff(r) -- computed at the non-relativistic
level (see S. Lehtola, JCTC 15, 1593 (2019)) -- from PySCF's bundled
``pyscf.dft.sap_data`` table and writes them to a plain-text data file that
OpenQP loads at runtime. PySCF is therefore only a build/data dependency here;
the resulting OpenQP SAP guess is fully native Fortran at runtime.

The SAP potential of a neutral atom of nuclear charge Z is
    V_A(r) = -Z_eff(r) / r,
where Z_eff(0) = Z and Z_eff -> 0 as r -> infinity.

Usage:
    python3 tools/sap/generate_sap_data.py [zmax] [output_path]
Defaults: zmax=36 (H-Kr), output basis_sets/sap_grasp.dat
"""

import sys
import os
import numpy as np


def main():
    zmax = int(sys.argv[1]) if len(sys.argv) > 1 else 36
    here = os.path.dirname(os.path.abspath(__file__))
    default_out = os.path.normpath(
        os.path.join(here, os.pardir, os.pardir, "basis_sets", "sap_grasp.dat"))
    out_path = sys.argv[2] if len(sys.argv) > 2 else default_out

    from pyscf.dft.sap_data import sap_Zeff

    nr = sap_Zeff.shape[1]
    if zmax >= sap_Zeff.shape[0]:
        raise ValueError(f"zmax={zmax} exceeds available data (max {sap_Zeff.shape[0]-1})")

    r = sap_Zeff[0, :]
    with open(out_path, "w") as f:
        f.write("# OpenQP radial SAP effective-charge table Z_eff(r)\n")
        f.write("# Source: Lehtola, JCTC 15, 1593 (2019), via pyscf.dft.sap_data\n")
        f.write("# Potential of atom Z at distance r:  V(r) = -Z_eff(r)/r\n")
        f.write("# Format: line 1 = 'zmax nr'; then nr r-values; then zmax rows of nr Z_eff values\n")
        f.write(f"{zmax} {nr}\n")
        _write_row(f, r)
        for z in range(1, zmax + 1):
            _write_row(f, sap_Zeff[z, :])

    print(f"Wrote {out_path}: zmax={zmax}, nr={nr}, r_max={r.max():.4f} bohr")


def _write_row(f, vals):
    # 6 values per line, full double precision
    for i in range(0, len(vals), 6):
        chunk = vals[i:i + 6]
        f.write(" ".join(f"{v:.16e}" for v in chunk) + "\n")


if __name__ == "__main__":
    main()
