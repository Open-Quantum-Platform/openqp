#!/usr/bin/env python3
"""Generate the OpenQP MINAO atomic-density data file.

One-time, offline data-generation tool. Builds spherically-averaged neutral
atomic density matrices in a fixed minimal reference basis (STO-3G, Cartesian)
using PySCF's spherically-averaged atomic-HF density (scf.hf.init_guess_by_atom,
the same densities PySCF uses for its SAD guess). The resulting densities are written to a plain-text data
file that OpenQP loads at runtime and projects onto the target basis. PySCF is
therefore only a build/data dependency; the OpenQP MINAO guess is native at
runtime.

Cartesian AO ordering matches PySCF's, which is interoperable with OpenQP (the
existing PySCF guess path relies on the same compatibility).

Usage: python3 tools/minao/generate_minao_data.py [zmax] [output_path]
Defaults: zmax=36 (H-Kr), output basis_sets/minao_sto3g.dat
"""

import sys
import os
import numpy as np


SYMBOLS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
    'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
    'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
]


def atomic_density(sym, Z, basis='sto-3g'):
    from pyscf import gto, scf
    mol = gto.M(atom=f'{sym} 0 0 0', basis=basis, spin=Z % 2, verbose=0)
    mol.cart = True
    mol.build()
    # PySCF's spherically-averaged atomic-HF density (used by its SAD guess),
    # returned in this atom's AO basis. May be spin-resolved (2, nao, nao).
    dm = np.asarray(scf.hf.init_guess_by_atom(mol))
    if dm.ndim == 3:
        dm = dm.sum(axis=0)
    S = mol.intor('int1e_ovlp')
    nelec = float(np.einsum('ij,ji->', dm, S))
    return dm, mol.nao_nr(), nelec


def main():
    zmax = int(sys.argv[1]) if len(sys.argv) > 1 else 36
    here = os.path.dirname(os.path.abspath(__file__))
    default_out = os.path.normpath(
        os.path.join(here, os.pardir, os.pardir, "basis_sets", "minao_sto3g.dat"))
    out_path = sys.argv[2] if len(sys.argv) > 2 else default_out

    with open(out_path, "w") as f:
        f.write("# OpenQP MINAO atomic densities (spherically-averaged neutral atoms)\n")
        f.write("# Reference basis: STO-3G (Cartesian). Source: pyscf scf.hf.init_guess_by_atom\n")
        f.write("# Format: line 1 = 'zmax'; then per element: 'Z nao nelec' then nao*nao\n")
        f.write("#         density-matrix values (row-major, AO basis).\n")
        f.write(f"{zmax}\n")
        for z in range(1, zmax + 1):
            sym = SYMBOLS[z - 1]
            dm, nao, nelec = atomic_density(sym, z)
            f.write(f"{z} {nao} {nelec:.10f}\n")
            flat = dm.reshape(-1)
            for i in range(0, len(flat), 6):
                f.write(" ".join(f"{v:.16e}" for v in flat[i:i+6]) + "\n")
            print(f"{sym:3s} Z={z:2d} nao={nao:2d} nelec={nelec:.4f}")

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
