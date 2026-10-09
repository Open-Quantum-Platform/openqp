#!/usr/bin/env python3
"""Generate the OpenQP MINAO atomic-density data file.

One-time, offline data-generation tool. Builds spherically-averaged neutral
atomic density matrices in a fixed minimal reference basis (STO-3G, Cartesian)
using PySCF's spherically-averaged atomic-HF density (scf.hf.init_guess_by_atom,
the same densities PySCF uses for its SAD guess). The STO-3G exponents and
coefficients are read from OpenQP's own basis_sets/sto-3g.basis, so the table is
consistent with the minimal basis the runtime projects from (PySCF's bundled
STO-3G stops at I; OpenQP's file covers H-Xe). The resulting densities are written to a plain-text data
file that OpenQP loads at runtime and projects onto the target basis. PySCF is
therefore only a build/data dependency; the OpenQP MINAO guess is native at
runtime.

The densities are written in OpenQP's Cartesian AO convention: shells in file
order (s, then p, then d, as PySCF also orders them), d components in the order
xx, yy, zz, xy, xz, yz, and every component normalized to one. PySCF orders d
as xx, xy, xz, yy, yz, zz and normalizes xx to 4*pi/5 and xy to 4*pi/15
(2.513 and 0.838 for Br); tables written in that convention put the d
densities on the wrong functions (Br held 34.55 electrons instead of 35).

Usage: python3 tools/minao/generate_minao_data.py [zmax] [output_path]
Defaults: zmax=54 (H-Xe, the range of OpenQP's STO-3G), output
basis_sets/minao_sto3g.dat
"""

import sys
import os
import numpy as np


SYMBOLS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
    'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
    'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I', 'Xe',
]

ANGULAR = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}


def read_openqp_basis(path):
    """Read OpenQP's GAMESS-format basis file into PySCF shell lists.

    Returns {Z: [[l, [exp, coef], ...], ...]} in file order, which is the AO
    order of OpenQP's minimal basis (STO-3G is listed s shells first, then p,
    then d, the same order PySCF uses).
    """
    out, z, shell, nleft = {}, None, None, 0
    body = open(path).read().split('$DATA', 1)[1]
    for raw in body.splitlines():
        line = raw.strip()
        if line.startswith('!') or line.startswith('$END'):
            continue
        if not line:
            z = None
            continue
        tok = line.split()
        if z is None:
            z = LONG_NAMES.index(line.upper()) + 1
            out[z] = []
        elif nleft == 0:
            shell = [ANGULAR[tok[0].upper()]]
            nleft = int(tok[1])
            out[z].append(shell)
        else:
            shell.append([float(tok[1].replace('D', 'E')),
                          float(tok[2].replace('D', 'E'))])
            nleft -= 1
    return out


LONG_NAMES = [
    'HYDROGEN', 'HELIUM', 'LITHIUM', 'BERYLLIUM', 'BORON', 'CARBON', 'NITROGEN',
    'OXYGEN', 'FLUORINE', 'NEON', 'SODIUM', 'MAGNESIUM', 'ALUMINIUM', 'SILICON',
    'PHOSPHORUS', 'SULFUR', 'CHLORINE', 'ARGON', 'POTASSIUM', 'CALCIUM',
    'SCANDIUM', 'TITANIUM', 'VANADIUM', 'CHROMIUM', 'MANGANESE', 'IRON', 'COBALT',
    'NICKEL', 'COPPER', 'ZINC', 'GALLIUM', 'GERMANIUM', 'ARSENIC', 'SELENIUM',
    'BROMINE', 'KRYPTON', 'RUBIDIUM', 'STRONTIUM', 'YTTRIUM', 'ZIRCONIUM',
    'NIOBIUM', 'MOLYBDENUM', 'TECHNETIUM', 'RUTHENIUM', 'RHODIUM', 'PALLADIUM',
    'SILVER', 'CADMIUM', 'INDIUM', 'TIN', 'ANTIMONY', 'TELLURIUM', 'IODINE',
    'XENON',
]


def atomic_density(sym, Z, basis):
    from pyscf import gto, scf
    mol = gto.M(atom=f'{sym} 0 0 0', basis={sym: basis}, spin=Z % 2, verbose=0)
    mol.cart = True
    mol.build()
    # PySCF's spherically-averaged atomic-HF density (used by its SAD guess),
    # returned in this atom's AO basis. May be spin-resolved (2, nao, nao).
    dm = np.asarray(scf.hf.init_guess_by_atom(mol))
    if dm.ndim == 3:
        dm = dm.sum(axis=0)
    S = mol.intor('int1e_ovlp')
    # PySCF -> OpenQP Cartesian convention: unit-normalized components and the
    # d order xx, yy, zz, xy, xz, yz (OpenQP slot k holds PySCF component perm[k])
    norm = np.sqrt(np.diag(S))
    perm = np.arange(mol.nao_nr())
    loc = mol.ao_loc_nr()
    for ib in range(mol.nbas):
        l = mol.bas_angular(ib)
        if l == 2:
            perm[loc[ib]:loc[ib] + 6] = loc[ib] + np.array([0, 3, 5, 1, 2, 4])
        elif l > 2:
            raise NotImplementedError("Cartesian order for l > 2 is not mapped")
    dm = (dm * np.outer(norm, norm))[np.ix_(perm, perm)]
    S = (S / np.outer(norm, norm))[np.ix_(perm, perm)]
    nelec = float(np.einsum('ij,ji->', dm, S))
    return dm, mol.nao_nr(), nelec


def main():
    zmax = int(sys.argv[1]) if len(sys.argv) > 1 else len(SYMBOLS)
    here = os.path.dirname(os.path.abspath(__file__))
    sto3g = read_openqp_basis(os.path.normpath(
        os.path.join(here, os.pardir, os.pardir, "basis_sets", "sto-3g.basis")))
    default_out = os.path.normpath(
        os.path.join(here, os.pardir, os.pardir, "basis_sets", "minao_sto3g.dat"))
    out_path = sys.argv[2] if len(sys.argv) > 2 else default_out

    with open(out_path, "w") as f:
        f.write("# OpenQP MINAO atomic densities (spherically-averaged neutral atoms)\n")
        f.write("# Reference basis: STO-3G (Cartesian) from basis_sets/sto-3g.basis. Source: pyscf scf.hf.init_guess_by_atom\n")
        f.write("# AO convention: OpenQP Cartesian (unit-normalized; d order xx,yy,zz,xy,xz,yz)\n")
        f.write("# Format: line 1 = 'zmax'; then per element: 'Z nao nelec' then nao*nao\n")
        f.write("#         density-matrix values (row-major, AO basis).\n")
        f.write(f"{zmax}\n")
        for z in range(1, zmax + 1):
            sym = SYMBOLS[z - 1]
            dm, nao, nelec = atomic_density(sym, z, sto3g[z])
            f.write(f"{z} {nao} {nelec:.10f}\n")
            flat = dm.reshape(-1)
            for i in range(0, len(flat), 6):
                f.write(" ".join(f"{v:.16e}" for v in flat[i:i+6]) + "\n")
            print(f"{sym:3s} Z={z:2d} nao={nao:2d} nelec={nelec:.4f}")

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
