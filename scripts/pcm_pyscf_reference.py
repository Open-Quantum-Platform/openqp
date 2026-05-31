#!/usr/bin/env python3
"""PySCF reference checks for OpenQP/ddX PCM source diagnostics.

This script is intentionally diagnostic-only. It does not configure OpenQP PCM
physics. It validates the vacuum RHF density/electrostatic-potential path used
as the reference for ddX source-term work:

* run OpenQP H2O RHF/6-31g* with PCM off and collect energy/density/overlap;
* run PySCF for the same molecule, cartesian 6-31g* basis;
* build the ddX cavity through the OpenQP C adapter;
* evaluate PySCF's exact QM electrostatic potential on those cavity points;
* optionally run OpenQP PCM diagnostics and compare the exact phi_cav summary.

Run from the repository root after a ddX-enabled OpenQP build/install, e.g.:

  OPENQP_ROOT=$PWD PYTHONPATH=$PWD/pyoqp \
  DYLD_LIBRARY_PATH=/tmp/ddx-install/lib:$PWD/lib:/tmp/oqp-build-pcm/source \
  OMP_NUM_THREADS=1 python scripts/pcm_pyscf_reference.py
"""

from __future__ import annotations

import argparse
import ctypes
import json
import math
import os
import re
import tempfile
from pathlib import Path

import numpy as np
from pyscf import gto, scf, solvent

from oqp.pyoqp import Runner

ANGSTROM_TO_BOHR = 1.8897259885789233
DEFAULT_GEOM_ANG = [
    ("O", (0.000000000, 0.000000000, -0.041061554)),
    ("H", (-0.533194329, 0.533194329, -0.614469223)),
    ("H", (0.533194329, -0.533194329, -0.614469223)),
]
NUCLEAR_CHARGES = {"H": 1.0, "O": 8.0}


def write_openqp_input(path: Path, pcm_enabled: bool) -> None:
    pcm_line = "enabled=true\nbackend=ddx\nmode=reference_scf\nmodel=ddpcm\nepsilon=78.3553" if pcm_enabled else "enabled=false"
    system_lines = "\n".join(
        f"   {int(NUCLEAR_CHARGES[sym])} {xyz[0]: .9f} {xyz[1]: .9f} {xyz[2]: .9f}"
        for sym, xyz in DEFAULT_GEOM_ANG
    )
    path.write_text(
        f"""[input]
system=
{system_lines}
charge=0
runtype=energy
basis=6-31g*
method=hf

[guess]
type=huckel

[scf]
multiplicity=1
type=rhf
conv=1.0e-10

[pcm]
{pcm_line}
""",
        encoding="utf-8",
    )


def unpack_lower(packed: list[float] | np.ndarray, n: int) -> np.ndarray:
    packed = np.asarray(packed, dtype=float)
    mat = np.zeros((n, n), dtype=float)
    k = 0
    for i in range(n):
        for j in range(i + 1):
            mat[i, j] = packed[k]
            mat[j, i] = packed[k]
            k += 1
    return mat


def parse_pcm_diag(log_text: str) -> dict[str, float | int | str]:
    out: dict[str, float | int | str] = {}
    for key in (
        "e_pcm",
        "half_tr_dv",
        "q_cav_sum",
        "q_cav_absnorm",
        "source_charge_sum",
        "phi_source_vs_exact_rms",
        "phi_source_vs_exact_max",
        "phi_cav_sum",
        "phi_cav_min",
        "phi_cav_max",
    ):
        m = re.findall(rf"PCM diag {key}=\s*(-?\d+\.\d+(?:[eE][-+]?\d+)?)", log_text)
        if m:
            out[key] = float(m[-1])
    m = re.findall(r"PCM diag ncav=\s*(\d+)", log_text)
    if m:
        out["ncav"] = int(m[-1])
    m = re.findall(r"PCM diag psi_source=(\S+)", log_text)
    if m:
        out["psi_source"] = m[-1]
    return out


def run_openqp(inp: Path, project: Path) -> tuple[dict, str]:
    runner = Runner(
        project=str(project),
        input_file=str(inp),
        log=str(project.with_suffix(".log")),
        silent=1,
        usempi=False,
    )
    runner.run(test_mod=True)
    return runner.results(), project.with_suffix(".log").read_text(encoding="utf-8")


def run_pyscf() -> tuple[gto.Mole, scf.hf.RHF, np.ndarray]:
    mol = gto.M(
        atom=DEFAULT_GEOM_ANG,
        unit="Angstrom",
        basis="6-31g*",
        charge=0,
        spin=0,
        cart=True,
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.conv_tol = 1.0e-12
    mf.kernel()
    return mol, mf, mf.make_rdm1()


def run_pyscf_pcm_reference() -> dict[str, float | str]:
    """Run PySCF's PCM wrapper as a protocol-near, not identical, reference."""
    mol = gto.M(
        atom=DEFAULT_GEOM_ANG,
        unit="Angstrom",
        basis="6-31g*",
        charge=0,
        spin=0,
        cart=True,
        verbose=0,
    )
    vac = scf.RHF(mol)
    vac.conv_tol = 1.0e-12
    e_vac = float(vac.kernel())
    pcm_mf = solvent.PCM(scf.RHF(mol))
    pcm_mf.conv_tol = 1.0e-12
    pcm_mf.with_solvent.eps = 78.3553
    e_pcm_total = float(pcm_mf.kernel())
    solvent_energy = getattr(pcm_mf.with_solvent, "e", None)
    return {
        "model": "PySCF solvent.PCM; cavity/protocol not identical to OpenQP/ddX",
        "vacuum_total_energy_hartree": e_vac,
        "pcm_total_energy_hartree": e_pcm_total,
        "pcm_total_minus_vacuum_hartree": e_pcm_total - e_vac,
        "with_solvent_e_hartree": float(solvent_energy) if solvent_energy is not None else float("nan"),
    }


def ddx_cavity(lib_path: Path, xyz_bohr: np.ndarray, charges: np.ndarray, epsilon: float) -> np.ndarray:
    lib = ctypes.CDLL(str(lib_path))
    fn = lib.oqp_ddx_pcm_cavity
    fn.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_double,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_char_p,
        ctypes.c_int,
    ]
    max_cav = 20000
    ncav = ctypes.c_int(0)
    cav = np.zeros(3 * max_cav, dtype=np.double)
    msg = ctypes.create_string_buffer(256)
    rc = fn(
        len(charges),
        np.ascontiguousarray(xyz_bohr, dtype=np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        np.ascontiguousarray(charges, dtype=np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_double(epsilon),
        max_cav,
        ctypes.byref(ncav),
        cav.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        msg,
        256,
    )
    if rc != 0:
        raise RuntimeError(f"oqp_ddx_pcm_cavity failed rc={rc}: {msg.value.decode(errors='replace')}")
    return cav[: 3 * ncav.value].reshape(ncav.value, 3).copy()


def pyscf_total_esp_on_points(mol: gto.Mole, dm: np.ndarray, points_bohr: np.ndarray) -> np.ndarray:
    ints = mol.intor("int1e_grids", grids=points_bohr)
    if ints.shape[0] == points_bohr.shape[0]:
        phi_elec = np.einsum("gij,ij->g", ints, dm)
    else:
        phi_elec = np.einsum("ijg,ij->g", ints, dm)
    atom_coords = mol.atom_coords(unit="Bohr")
    charges = mol.atom_charges().astype(float)
    phi_nuc = np.zeros(points_bohr.shape[0], dtype=float)
    for charge, coord in zip(charges, atom_coords):
        phi_nuc += charge / np.linalg.norm(points_bohr - coord, axis=1)
    return phi_nuc - phi_elec


def summary(values: np.ndarray) -> dict[str, float]:
    return {
        "sum": float(np.sum(values)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "rms": float(math.sqrt(float(np.mean(values * values)))),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--liboqp", default="lib/liboqp.dylib")
    parser.add_argument("--json-out", default="")
    parser.add_argument("--skip-openqp-pcm", action="store_true")
    args = parser.parse_args()

    repo = Path.cwd()
    liboqp = (repo / args.liboqp).resolve()
    symbols = [sym for sym, _ in DEFAULT_GEOM_ANG]
    xyz_ang = np.array([coord for _, coord in DEFAULT_GEOM_ANG], dtype=float)
    xyz_bohr = (xyz_ang * ANGSTROM_TO_BOHR).reshape(-1)
    charges = np.array([NUCLEAR_CHARGES[sym] for sym in symbols], dtype=float)

    with tempfile.TemporaryDirectory(prefix="oqp_pcm_pyscf_") as tmp_s:
        tmp = Path(tmp_s)
        vac_inp = tmp / "h2o_vac.inp"
        pcm_inp = tmp / "h2o_pcm.inp"
        write_openqp_input(vac_inp, pcm_enabled=False)
        write_openqp_input(pcm_inp, pcm_enabled=True)

        oqp_vac, _ = run_openqp(vac_inp, tmp / "h2o_vac")
        nbf = int(round((math.sqrt(8 * len(oqp_vac["data"]["OQP::SM"]) + 1) - 1) / 2))
        s_oqp = unpack_lower(oqp_vac["data"]["OQP::SM"], nbf)
        d_oqp = unpack_lower(oqp_vac["data"]["OQP::DM_A"], nbf)

        mol, mf, dm_py = run_pyscf()
        s_py = mol.intor("int1e_ovlp")
        cavity = ddx_cavity(liboqp, xyz_bohr, charges, epsilon=78.3553)
        phi_py = pyscf_total_esp_on_points(mol, dm_py, cavity)

        phi_py_summary = summary(phi_py)
        report: dict[str, object] = {
            "system": "H2O RHF/6-31g* cartesian basis, epsilon=78.3553",
            "openqp_energy_hartree": float(oqp_vac["energy"][0]),
            "pyscf_energy_hartree": float(mf.e_tot),
            "abs_energy_delta_hartree": float(abs(float(oqp_vac["energy"][0]) - float(mf.e_tot))),
            "nao_openqp": nbf,
            "nao_pyscf": int(mol.nao_nr()),
            "electron_count_openqp_trace_DS": float(np.einsum("ij,ij->", d_oqp, s_oqp)),
            "electron_count_pyscf_trace_DS": float(np.einsum("ij,ij->", dm_py, s_py)),
            "note_basis_order": "OpenQP and PySCF AO ordering differ; direct matrix max-diff is not a validation metric here.",
            "ddx_ncav": int(cavity.shape[0]),
            "pyscf_exact_phi_cav": phi_py_summary,
            "pyscf_pcm_reference": run_pyscf_pcm_reference(),
        }

        if not args.skip_openqp_pcm:
            _, pcm_log = run_openqp(pcm_inp, tmp / "h2o_pcm")
            diag = parse_pcm_diag(pcm_log)
            report["openqp_pcm_diag"] = diag
            if {"phi_cav_sum", "phi_cav_min", "phi_cav_max", "ncav"}.issubset(diag):
                report["openqp_vs_pyscf_phi_cav_summary_delta"] = {
                    "sum": float(float(diag["phi_cav_sum"]) - phi_py_summary["sum"]),
                    "min": float(float(diag["phi_cav_min"]) - phi_py_summary["min"]),
                    "max": float(float(diag["phi_cav_max"]) - phi_py_summary["max"]),
                    "ncav_equal": bool(int(diag["ncav"]) == int(cavity.shape[0])),
                }

    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.json_out:
        Path(args.json_out).write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
