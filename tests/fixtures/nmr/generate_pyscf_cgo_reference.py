#!/usr/bin/env python3
"""Generate the PySCF common-gauge-origin (CGO) NMR shielding reference fixture.

This is the gate-4 oracle for the OpenQP coupled magnetic response (Phase 0).
It produces, for H2O / STO-3G with the gauge origin at the (OpenQP) center of
mass, the isotropic diamagnetic and paramagnetic (uncoupled AND coupled) NMR
shieldings for HF, BHHLYP, PBE0, and PBE.

The coupled paramagnetic response is solved by a damped fixed-point iteration of
the magnetic CPHF/CPKS equations using PySCF's own `gen_response(hermi=2)`
(anti-Hermitian density response = exact-exchange only; Coulomb and the
semi-local XC kernel vanish for the imaginary antisymmetric magnetic density).
This avoids a version-specific reshape bug in `pyscf.scf.cphf.solve` when a
common gauge origin is used.

Requires: pyscf, pyscf-properties. Run:
    python generate_pyscf_cgo_reference.py
It writes pyscf_cgo_reference.json next to this script. The JSON is the
committed oracle; do not recompute these numbers ad hoc elsewhere.
"""
import json
import os
import numpy as np
from pyscf import gto, scf, dft
from pyscf.prop.nmr import rhf as nmr_rhf, rks as nmr_rks
from pyscf.data import nist

# Geometry (Angstrom) and gauge origin (Bohr) match the OpenQP example/test.
ATOM = """O  0.000000000  0.000000000 -0.041061554
          H -0.533194329  0.533194329 -0.614469223
          H  0.533194329 -0.533194329 -0.614469223"""
BASIS = "sto-3g"
GAUGE_ORIGIN_BOHR = [0.0, 0.0, -0.19886422]   # OpenQP center of mass
UNIT = nist.ALPHA ** 2 * 1e6                    # a.u. -> ppm


def _coupled_mo1(mf, h1, s1):
    """Damped fixed-point solve of the magnetic CPHF (hermi=2 exchange response)."""
    C, e, occ = mf.mo_coeff, mf.mo_energy, mf.mo_occ
    nmo = len(e); vir = occ == 0; oc = occ > 0
    e_ai = 1.0 / (e[vir].reshape(-1, 1) - e[oc])
    vind = nmr_rhf.gen_vind(mf, C, occ)
    mo1 = np.zeros((3, nmo, oc.sum()))
    mo1[:, vir, :] = -h1[:, vir, :] * e_ai
    for _ in range(100):
        v = vind(mo1.ravel()).reshape(3, nmo, oc.sum())
        new = mo1.copy()
        new[:, vir, :] = -(h1[:, vir, :] + v[:, vir, :]) * e_ai
        if np.max(np.abs(new - mo1)) < 1e-10:
            mo1 = new
            break
        mo1 = new
    return mo1


def _shielding(mf, NMRcls):
    C, occ = mf.mo_coeff, mf.mo_occ
    oidx = occ > 0
    O = np.array(GAUGE_ORIGIN_BOHR)
    nmrobj = NMRcls(mf); nmrobj.gauge_orig = O
    h1ao = nmr_rhf.make_h10(nmrobj.mol, mf.make_rdm1(), gauge_orig=O)
    s1ao = nmr_rhf.make_s10(nmrobj.mol, gauge_orig=O)
    h1 = np.einsum('pi,xpq,qj->xij', C, h1ao, C[:, oidx])
    s1 = np.einsum('pi,xpq,qj->xij', C, s1ao, C[:, oidx])
    mo_unc, _ = nmr_rhf._solve_mo1_uncoupled(mf.mo_energy, occ, h1, s1)
    mo_cpl = _coupled_mo1(mf, h1, s1)
    dia = nmrobj.dia(gauge_orig=O) * UNIT
    pu = nmrobj.para(mo10=mo_unc)[0] * UNIT
    pc = nmrobj.para(mo10=mo_cpl)[0] * UNIT
    iso = lambda M: [float(np.trace(M[i]) / 3) for i in range(mf.mol.natm)]
    return iso(dia), iso(pu), iso(pc)


def main():
    mol = gto.M(atom=ATOM, basis=BASIS, unit="Angstrom", charge=0, spin=0)
    methods = [("HF", 1.0, None), ("bhandhlyp", 0.5, "bhandhlyp"),
               ("pbe0", 0.25, "pbe0"), ("pbe", 0.0, "pbe")]
    out = {
        "description": "PySCF common-gauge-origin NMR shielding reference "
                       "(gate-4 oracle for OpenQP coupled magnetic response).",
        "molecule": "H2O", "basis": BASIS, "unit_coords": "Angstrom",
        "geometry": ATOM, "gauge_origin_bohr": GAUGE_ORIGIN_BOHR,
        "shielding_unit": "ppm", "atom_order": ["O", "H", "H"],
        "pyscf_version": __import__("pyscf").__version__,
        "results": {},
    }
    for label, cx, xc in methods:
        if xc is None:
            mf = scf.RHF(mol).run(); cls = nmr_rhf.NMR
        else:
            mf = dft.RKS(mol); mf.xc = xc; mf.run(); cls = nmr_rks.NMR
        dia, pu, pc = _shielding(mf, cls)
        out["results"][label] = {
            "exact_exchange_fraction": cx,
            "sigma_dia": [round(v, 4) for v in dia],
            "sigma_para_uncoupled": [round(v, 4) for v in pu],
            "sigma_para_coupled": [round(v, 4) for v in pc],
            "sigma_total_uncoupled": [round(d + u, 4) for d, u in zip(dia, pu)],
            "sigma_total_coupled": [round(d + c, 4) for d, c in zip(dia, pc)],
        }
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "pyscf_cgo_reference.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    print("wrote", path)
    print(json.dumps(out["results"], indent=2))


if __name__ == "__main__":
    main()
