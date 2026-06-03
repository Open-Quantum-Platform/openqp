#!/usr/bin/env python3
"""Generate PySCF GIAO NMR shielding reference data for OpenQP validation.

This generator is intentionally separate from generate_pyscf_cgo_reference.py.
The CGO fixture sets ``gauge_orig`` to a fixed common origin; this fixture leaves
``gauge_orig=None`` so PySCF uses its London-orbital/GIAO integral path.

The output is an oracle for validating OpenQP's native GIAO implementation. It is
not an implementation path for OpenQP. The production OpenQP feature remains
unvalidated/gated until native GIAO integrals, response equations, origin-
translation tests, and benchmark comparisons pass.

Run:
    python tests/fixtures/nmr/generate_pyscf_giao_reference.py

Writes:
    tests/fixtures/nmr/pyscf_giao_reference.json
"""
from __future__ import annotations

import importlib.metadata
import json
import os
import time
import warnings
from dataclasses import dataclass
from typing import Iterable

import numpy as np
from pyscf import dft, gto, scf
from pyscf.data import nist
from pyscf.prop.nmr import rhf as nmr_rhf, rks as nmr_rks

warnings.filterwarnings("ignore", category=UserWarning, module=r"pyscf\.prop\..*")

ATOM = """O  0.000000000  0.000000000 -0.041061554
          H -0.533194329  0.533194329 -0.614469223
          H  0.533194329 -0.533194329 -0.614469223"""
BASIS = "sto-3g"
UNIT = nist.ALPHA ** 2 * 1e6


@dataclass(frozen=True)
class MethodCase:
    label: str
    exact_exchange_fraction: float
    xc: str | None


CASES = (
    MethodCase("HF", 1.0, None),
    MethodCase("bhandhlyp", 0.5, "bhandhlyp"),
    MethodCase("pbe0", 0.25, "pbe0"),
    MethodCase("pbe", 0.0, "pbe"),
)


def _coupled_mo1(mf, h1, s1, max_cycle: int = 100, conv_tol: float = 1e-10):
    """Fixed-point magnetic CPHF/CPKS solve for PySCF NMR h1.

    PySCF/property 0.1.0 with PySCF 2.13 can hit a reshape issue in the stock
    CPHF solve. This mirrors the existing CGO fixture generator and uses
    PySCF's own ``gen_vind`` response operator. The h1/s1 AO integrals are still
    PySCF's GIAO integrals because ``gauge_orig=None`` was used upstream.
    """
    coeff, energy, occ = mf.mo_coeff, mf.mo_energy, mf.mo_occ
    nmo = len(energy)
    occ_mask = occ > 0
    vir_mask = occ == 0
    e_ai = 1.0 / (energy[vir_mask].reshape(-1, 1) - energy[occ_mask])
    vind = nmr_rhf.gen_vind(mf, coeff, occ)
    hs = h1 - s1 * energy[occ_mask]
    mo1 = hs.copy()
    mo1[:, vir_mask, :] = -hs[:, vir_mask, :] * e_ai
    mo1[:, occ_mask, :] = -0.5 * s1[:, occ_mask, :]
    for _ in range(max_cycle):
        v1 = vind(mo1.ravel()).reshape(3, nmo, occ_mask.sum())
        new = mo1.copy()
        new[:, vir_mask, :] = -(hs[:, vir_mask, :] + v1[:, vir_mask, :]) * e_ai
        new[:, occ_mask, :] = -0.5 * s1[:, occ_mask, :]
        if np.max(np.abs(new - mo1)) < conv_tol:
            return new
        mo1 = new
    raise RuntimeError("PySCF GIAO CPHF reference did not converge")


def _shielding(mf, nmr_cls):
    coeff = mf.mo_coeff
    occ_mask = mf.mo_occ > 0
    nmrobj = nmr_cls(mf)
    nmrobj.gauge_orig = None
    h1ao = nmr_rhf.make_h10(nmrobj.mol, mf.make_rdm1(), gauge_orig=None)
    s1ao = nmr_rhf.make_s10(nmrobj.mol, gauge_orig=None)
    h1 = np.einsum("pi,xpq,qj->xij", coeff, h1ao, coeff[:, occ_mask])
    s1 = np.einsum("pi,xpq,qj->xij", coeff, s1ao, coeff[:, occ_mask])
    mo_uncoupled, _ = nmr_rhf._solve_mo1_uncoupled(mf.mo_energy, mf.mo_occ, h1, s1)
    mo_coupled = _coupled_mo1(mf, h1, s1)
    dia = nmrobj.dia(gauge_orig=None) * UNIT
    para_uncoupled = nmrobj.para(mo10=mo_uncoupled)[0] * UNIT
    para_coupled = nmrobj.para(mo10=mo_coupled)[0] * UNIT
    iso = lambda tensor: [float(np.trace(tensor[i]) / 3.0) for i in range(mf.mol.natm)]
    return {
        "sigma_dia_tensor": dia.tolist(),
        "sigma_para_uncoupled_tensor": para_uncoupled.tolist(),
        "sigma_para_coupled_tensor": para_coupled.tolist(),
        "sigma_total_uncoupled_tensor": (dia + para_uncoupled).tolist(),
        "sigma_total_coupled_tensor": (dia + para_coupled).tolist(),
        "sigma_dia": iso(dia),
        "sigma_para_uncoupled": iso(para_uncoupled),
        "sigma_para_coupled": iso(para_coupled),
        "sigma_total_uncoupled": iso(dia + para_uncoupled),
        "sigma_total_coupled": iso(dia + para_coupled),
    }


def _run_case(mol, case: MethodCase):
    start = time.perf_counter()
    if case.xc is None:
        mf = scf.RHF(mol).run(verbose=0)
        nmr_cls = nmr_rhf.NMR
    else:
        mf = dft.RKS(mol)
        mf.xc = case.xc
        mf.run(verbose=0)
        nmr_cls = nmr_rks.NMR
    result = _shielding(mf, nmr_cls)
    result["wall_time_seconds"] = round(time.perf_counter() - start, 6)
    result["exact_exchange_fraction"] = case.exact_exchange_fraction
    return result


def generate(cases: Iterable[MethodCase] = CASES):
    mol = gto.M(atom=ATOM, basis=BASIS, unit="Angstrom", charge=0, spin=0)
    out = {
        "description": "PySCF GIAO/London-orbital NMR shielding reference for OpenQP native GIAO validation.",
        "provenance": "Generated with pyscf.prop.nmr using gauge_orig=None; use as oracle only, not as OpenQP implementation.",
        "molecule": "H2O",
        "basis": BASIS,
        "basis_format": "sto-3g",
        "unit_coords": "Angstrom",
        "geometry": ATOM,
        "gauge_formulation": "giao",
        "gauge_origin": None,
        "shielding_unit": "ppm",
        "atom_order": ["O", "H", "H"],
        "pyscf_version": __import__("pyscf").__version__,
        "pyscf_properties_version": importlib.metadata.version("pyscf-properties"),
        "results": {},
    }
    for case in cases:
        out["results"][case.label] = _run_case(mol, case)
    return out


def main():
    out = generate()
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pyscf_giao_reference.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")
    print("wrote", path)
    compact = {
        method: {
            key: [round(v, 4) for v in vals]
            for key, vals in data.items()
            if key.startswith("sigma_") and key.endswith(("coupled", "uncoupled", "dia"))
        }
        for method, data in out["results"].items()
    }
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
