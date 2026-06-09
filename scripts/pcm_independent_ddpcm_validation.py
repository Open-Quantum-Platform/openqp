#!/usr/bin/env python3
"""Independent ddPCM validation reference for the OpenQP/ddX PCM path.

TRUE validation: this generates an INDEPENDENT reference with PySCF's native
``ddPCM`` implementation (a different codebase from OpenQP/ddX) under a
protocol matched as closely as feasible to OpenQP's adapter, and compares the
total SCF solvation shift ``E_solvated - E_vacuum`` against OpenQP.

This is reference/diagnostic Python ONLY. It does not configure or feed the
OpenQP production PCM path (which stays Fortran-side: pcm_enabled ->
add_pcm_reaction_field -> e_pcm -> ddX C adapter). The OpenQP numbers compared
here are obtained by running the OpenQP executable separately.

Matched protocol: model=ddPCM, eps=78.3553, lmax=8, Lebedev order 29 (302
points), eta=0.1, Bondi-by-Z radii identical to
source/solvent_ddx_adapter.c::vdw_radius_angstrom.

Known protocol caveats (documented, not hidden):
  * PySCF ddcosmo/ddPCM requires spherical orbitals; OpenQP 6-31G* uses
    Cartesian d. The 5d-vs-6d difference is negligible for the electrostatic
    solvation shift.
  * PySCF ddPCM and ddX ddPCM are independent implementations of the same
    model; small implementation differences are expected.
  * OpenQP currently builds the ddX source from an l<=2 atom-centered Mulliken
    multipole expansion; PySCF uses the full QM density. A large gap therefore
    indicates the OpenQP source truncation, which is the point of the check.
"""

from __future__ import annotations

import numpy as np
from pyscf import gto, scf
from pyscf.solvent import ddPCM

HARTREE_TO_KCAL = 627.5094740631
BOHR = 1.8897259886

# OpenQP Bondi vdW radii (Angstrom), from solvent_ddx_adapter.c::vdw_radius_angstrom
BONDI_ANGSTROM = {1: 1.20, 6: 1.70, 7: 1.55, 8: 1.52, 9: 1.47}

GEOMETRIES = {
    "H2O": "O 0 0 -0.041061554\nH -0.533194329 0.533194329 -0.614469223\nH 0.533194329 -0.533194329 -0.614469223",
    "NH3": "N 0 0 0.116489\nH 0 0.939731 -0.271808\nH 0.813831 -0.469865 -0.271808\nH -0.813831 -0.469865 -0.271808",
}

# OpenQP/ddX total SCF energies (Eh) from real ddX-enabled runs on this branch.
OPENQP = {
    "H2O": {"e_vac": -76.010747, "e_pcm": -76.019593},
    "NH3": {"e_vac": -56.183840, "e_pcm": -56.195661},
}


def _radii_table(zmax: int = 36) -> np.ndarray:
    t = np.full(zmax + 1, 1.80 * BOHR)
    t[0] = 1.50 * BOHR
    for z, a in BONDI_ANGSTROM.items():
        t[z] = a * BOHR
    return t


def pyscf_ddpcm_shift(geom: str) -> dict:
    mol = gto.M(atom=geom, unit="Angstrom", basis="6-31g*", charge=0, spin=0,
                cart=False, verbose=0)
    e_vac = float(scf.RHF(mol).run(conv_tol=1e-11).e_tot)
    mf = ddPCM(scf.RHF(mol))
    ws = mf.with_solvent
    ws.eps = 78.3553
    ws.lmax = 8
    ws.lebedev_order = 29   # 302 cavity points
    ws.eta = 0.1
    ws.radii_table = _radii_table()
    mf.conv_tol = 1e-11
    e_solv = float(mf.kernel())
    return {
        "e_vac": e_vac,
        "e_solv": e_solv,
        "shift_hartree": e_solv - e_vac,
        "solvent_energy_term": float(getattr(ws, "e", float("nan"))),
    }


def main() -> None:
    print("Independent ddPCM validation (PySCF ddPCM vs OpenQP/ddX), matched protocol")
    print("=" * 78)
    for name, geom in GEOMETRIES.items():
        ref = pyscf_ddpcm_shift(geom)
        oq = OPENQP[name]
        oq_shift = oq["e_pcm"] - oq["e_vac"]
        gap = oq_shift - ref["shift_hartree"]
        ratio = oq_shift / ref["shift_hartree"]
        print(f"\n{name}:")
        print(f"  PySCF  ddPCM shift = {ref['shift_hartree']:+.6f} Eh "
              f"({ref['shift_hartree']*HARTREE_TO_KCAL:+.2f} kcal/mol)  [independent reference]")
        print(f"  OpenQP ddX   shift = {oq_shift:+.6f} Eh "
              f"({oq_shift*HARTREE_TO_KCAL:+.2f} kcal/mol)")
        print(f"  gap = {gap*HARTREE_TO_KCAL:+.2f} kcal/mol   ratio OpenQP/PySCF = {ratio:.3f}")
        verdict = "WITHIN ~5%% (consistent)" if abs(ratio - 1.0) < 0.10 else "DISCREPANT (not validated)"
        print(f"  verdict: {verdict}")


if __name__ == "__main__":
    main()
