#!/usr/bin/env python3
"""Run the NMR CGO/GIAO benchmark matrix scaffold.

This is an integrated OpenQP benchmark driver, not a production GIAO
implementation. It records the current validation state explicitly:

* trusted-reference columns are generated with PySCF CGO/GIAO NMR;
* OpenQP CGO may be run through the normal driver when OPENQP_ROOT is available;
* OpenQP GIAO is expected to remain gated until native London integral/response
  terms are connected. A gated/NotImplemented result is recorded as such, never
  converted to CGO numbers.

The decisive metric is gauge-origin translation sensitivity. For PySCF, CGO is
sampled with several fixed origins while GIAO uses ``gauge_orig=None`` and should
be invariant to translation of the molecular coordinates. For OpenQP, the same
metric is populated only after native GIAO is available.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from pyscf import dft, gto, scf
from pyscf.data import nist
from pyscf.prop.nmr import rhf as nmr_rhf, rks as nmr_rks

ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "tests" / "fixtures" / "nmr" / "giao_benchmark_matrix.json"
OUTDIR = ROOT / "tests" / "fixtures" / "nmr" / "benchmark_results"
UNIT = nist.ALPHA ** 2 * 1e6

GEOMS_ANG = {
    "h2o": """O  0.000000000  0.000000000 -0.041061554
              H -0.533194329  0.533194329 -0.614469223
              H  0.533194329 -0.533194329 -0.614469223""",
    "hf": """H 0.000000000 0.000000000 0.000000000
             F 0.000000000 0.000000000 0.916800000""",
    "ch4": """C  0.000000000  0.000000000  0.000000000
              H  0.629118000  0.629118000  0.629118000
              H -0.629118000 -0.629118000  0.629118000
              H -0.629118000  0.629118000 -0.629118000
              H  0.629118000 -0.629118000 -0.629118000""",
    "formaldehyde": """C  0.000000000  0.000000000  0.000000000
                         O  0.000000000  0.000000000  1.205000000
                         H  0.942900000  0.000000000 -0.587600000
                         H -0.942900000  0.000000000 -0.587600000""",
    "benzene": """C  1.397000000  0.000000000 0.000000000
                  C  0.698500000  1.209837000 0.000000000
                  C -0.698500000  1.209837000 0.000000000
                  C -1.397000000  0.000000000 0.000000000
                  C -0.698500000 -1.209837000 0.000000000
                  C  0.698500000 -1.209837000 0.000000000
                  H  2.481000000  0.000000000 0.000000000
                  H  1.240500000  2.148584000 0.000000000
                  H -1.240500000  2.148584000 0.000000000
                  H -2.481000000  0.000000000 0.000000000
                  H -1.240500000 -2.148584000 0.000000000
                  H  1.240500000 -2.148584000 0.000000000""",
    "glycine": """N -1.453000000  0.000000000  0.000000000
                 C -0.248000000  0.000000000  0.000000000
                 C  0.520000000  1.238000000  0.000000000
                 O  1.733000000  1.231000000  0.000000000
                 O -0.170000000  2.334000000  0.000000000
                 H -1.990000000  0.833000000  0.000000000
                 H -1.990000000 -0.833000000  0.000000000
                 H -0.214000000 -0.548000000  0.935000000
                 H -0.214000000 -0.548000000 -0.935000000
                 H  0.337000000  3.142000000  0.000000000""",
}

METHODS = {
    "hf": (None, 1.0),
    "pbe": ("pbe", 0.0),
    "pbe0": ("pbe0", 0.25),
    "bhhlyp": ("bhandhlyp", 0.5),
}


def _mol(name: str, basis: str, translation_bohr=(0.0, 0.0, 0.0)):
    mol = gto.M(atom=GEOMS_ANG[name], basis=basis, unit="Angstrom", charge=0, spin=0, verbose=0)
    if any(abs(x) > 0 for x in translation_bohr):
        coords = mol.atom_coords(unit="Bohr") + np.asarray(translation_bohr, dtype=float)
        atoms = [(mol.atom_symbol(i), coords[i]) for i in range(mol.natm)]
        mol = gto.M(atom=atoms, basis=basis, unit="Bohr", charge=0, spin=0, verbose=0)
    return mol


def _coupled_mo1(mf, h1, s1, max_cycle=100, conv_tol=1e-10):
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
    raise RuntimeError("PySCF NMR CPHF reference did not converge")


def _pyscf_shielding(mol, method: str, gauge: str, origin_bohr=(0.0, 0.0, 0.0)):
    xc, cx = METHODS[method]
    t0 = time.perf_counter()
    if xc is None:
        mf = scf.RHF(mol).run(verbose=0)
        cls = nmr_rhf.NMR
    else:
        mf = dft.RKS(mol)
        mf.xc = xc
        mf.run(verbose=0)
        cls = nmr_rks.NMR
    nmrobj = cls(mf)
    gauge_orig = None if gauge == "giao" else np.asarray(origin_bohr, dtype=float)
    coeff = mf.mo_coeff
    occ_mask = mf.mo_occ > 0
    h1ao = nmr_rhf.make_h10(nmrobj.mol, mf.make_rdm1(), gauge_orig=gauge_orig)
    s1ao = nmr_rhf.make_s10(nmrobj.mol, gauge_orig=gauge_orig)
    h1 = np.einsum("pi,xpq,qj->xij", coeff, h1ao, coeff[:, occ_mask])
    s1 = np.einsum("pi,xpq,qj->xij", coeff, s1ao, coeff[:, occ_mask])
    mo_coupled = _coupled_mo1(mf, h1, s1)
    tensor = (nmrobj.dia(gauge_orig=gauge_orig) + nmrobj.para(mo10=mo_coupled)[0]) * UNIT
    wall = time.perf_counter() - t0
    iso = np.asarray([float(np.trace(tensor[i]) / 3.0) for i in range(mol.natm)])
    return {
        "tensor_components_ppm": tensor.tolist(),
        "isotropic_shielding_ppm": iso.tolist(),
        "wall_time_seconds": wall,
        "exact_exchange_fraction": cx,
    }


def _run_openqp(name: str, basis: str, method: str, gauge: str):
    root = os.environ.get("OPENQP_ROOT")
    if not root or not ((Path(root) / "lib" / "liboqp.dylib").exists() or (Path(root) / "lib" / "liboqp.so").exists()):
        return {"status": "skipped", "reason": "OPENQP_ROOT does not point to a built OpenQP package"}
    xc, _ = METHODS[method]
    functional = f"functional={method}\n" if xc else ""
    dftgrid = "[dftgrid]\nrad_type=becke\n\n" if xc else ""
    system = "\n".join(
        f"   {gto.charge(sym):d} {float(x): .10f} {float(y): .10f} {float(z): .10f}"
        for sym, (x, y, z) in gto.format_atom(GEOMS_ANG[name], unit=1)
    )
    text = (
        f"[input]\nsystem=\n{system}\ncharge=0\nruntype=energy\n{functional}basis={basis}\nmethod=hf\n\n"
        f"[guess]\ntype=huckel\n\n[scf]\nmultiplicity=1\ntype=rhf\n\n{dftgrid}"
        f"[properties]\nscf_prop=nmr\nnmr_gauge={gauge}\n"
    )
    with tempfile.TemporaryDirectory() as wd:
        inp = Path(wd) / f"{name}_{method}_{gauge}.inp"
        inp.write_text(text)
        env = dict(os.environ)
        env["OPENQP_ROOT"] = root
        env["PYTHONPATH"] = str(ROOT / "pyoqp") + os.pathsep + env.get("PYTHONPATH", "")
        t0 = time.perf_counter()
        proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(inp)], cwd=wd, env=env, capture_output=True, text=True)
        wall = time.perf_counter() - t0
        log = inp.with_suffix(".log")
        stdout_stderr = proc.stdout + proc.stderr + (log.read_text() if log.exists() else "")
        if "NotImplementedError" in stdout_stderr or "not yet implemented" in stdout_stderr:
            return {"status": "gated", "reason": "OpenQP native GIAO path is still NotImplemented", "wall_time_seconds": wall}
        if proc.returncode != 0:
            return {"status": "failed", "returncode": proc.returncode, "tail": stdout_stderr[-2000:], "wall_time_seconds": wall}
        rows = []
        for line in stdout_stderr.splitlines():
            m = re.match(r"\s*(\d+)\s+[\d.]+" + r"\s+(-?\d+\.\d+)" * 5 + r"\s*$", line)
            if m:
                vals = [float(m.group(j)) for j in range(2, 7)]
                rows.append(vals[4])
        return {"status": "ok", "isotropic_shielding_ppm": rows, "wall_time_seconds": wall}


def _write_outputs(records, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    json_path = outdir / "nmr_giao_benchmark_results.json"
    csv_path = outdir / "nmr_giao_benchmark_results.csv"
    md_path = outdir / "nmr_giao_benchmark_results.md"
    json_path.write_text(json.dumps(records, indent=2) + "\n")
    fields = ["system", "basis", "method", "backend", "gauge", "status", "max_origin_delta_ppm", "wall_time_seconds", "reference_delta_ppm", "notes"]
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for r in records["rows"]:
            writer.writerow({k: r.get(k, "") for k in fields})
    lines = ["# NMR CGO/GIAO benchmark matrix", "", f"Status: {records['status']}", "", "|System|Basis|Method|Backend|Gauge|Status|Max origin delta/ppm|Reference delta/ppm|Wall/s|Notes|", "|---|---:|---|---|---|---|---:|---:|---:|---|"]
    for r in records["rows"]:
        lines.append("|{system}|{basis}|{method}|{backend}|{gauge}|{status}|{max_origin_delta_ppm}|{reference_delta_ppm}|{wall_time_seconds}|{notes}|".format(**{k: r.get(k, "") for k in ["system","basis","method","backend","gauge","status","max_origin_delta_ppm","reference_delta_ppm","wall_time_seconds","notes"]}))
    md_path.write_text("\n".join(lines) + "\n")
    return json_path, csv_path, md_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true", help="run every matrix case; default runs a quick subset")
    parser.add_argument("--run-openqp", action="store_true", help="also run OpenQP CGO/GIAO through the normal driver")
    parser.add_argument("--outdir", default=str(OUTDIR))
    args = parser.parse_args()

    matrix = json.loads(MATRIX.read_text())
    translations = matrix["gauge_origin_translations_bohr"]
    cases = []
    sections = matrix["small_regression_systems"] + (matrix["origin_dependence_systems"] if args.all else [])
    if not args.all:
        sections = [matrix["small_regression_systems"][0], matrix["origin_dependence_systems"][0]]
    for item in sections:
        for method in item["methods"]:
            cases.append((item["name"], item["basis"], method.lower()))

    rows = []
    for name, basis, method in cases:
        if method not in METHODS:
            continue
        ref_giao0 = None
        for gauge in ("cgo", "giao"):
            vals = []
            tensors = None
            wall = 0.0
            status = "ok"
            notes = "PySCF trusted reference"
            for trans in translations:
                try:
                    mol = _mol(name, basis, translation_bohr=trans if gauge == "giao" else (0.0, 0.0, 0.0))
                    result = _pyscf_shielding(mol, method, gauge, origin_bohr=trans)
                    vals.append(np.asarray(result["isotropic_shielding_ppm"]))
                    tensors = result["tensor_components_ppm"]
                    wall += result["wall_time_seconds"]
                except Exception as exc:
                    status = "failed"
                    notes = type(exc).__name__ + ": " + str(exc)[:160]
                    break
            if vals:
                arr = np.vstack(vals)
                max_delta = float(np.max(np.ptp(arr, axis=0)))
                if gauge == "giao":
                    ref_giao0 = arr[0]
                ref_delta = "0.000000" if gauge == "giao" or ref_giao0 is None else f"{float(np.max(np.abs(arr[0] - ref_giao0))):.6f}"
            else:
                max_delta = ""
                ref_delta = ""
            rows.append({
                "system": name, "basis": basis, "method": method, "backend": "pyscf", "gauge": gauge,
                "status": status, "max_origin_delta_ppm": f"{max_delta:.6f}" if isinstance(max_delta, float) else max_delta,
                "reference_delta_ppm": ref_delta, "wall_time_seconds": f"{wall:.3f}", "notes": notes,
                "tensor_components_ppm": tensors,
            })
        if args.run_openqp:
            for gauge in ("cgo", "giao"):
                r = _run_openqp(name, basis, method, gauge)
                rows.append({
                    "system": name, "basis": basis, "method": method, "backend": "openqp", "gauge": gauge,
                    "status": r["status"], "max_origin_delta_ppm": "pending" if gauge == "giao" else "not_sampled",
                    "reference_delta_ppm": "pending", "wall_time_seconds": f"{r.get('wall_time_seconds', 0.0):.3f}",
                    "notes": r.get("reason", r.get("tail", ""))[:160],
                })
    records = {
        "status": "reference_ready_openqp_giao_gated",
        "matrix_file": str(MATRIX),
        "metrics": matrix["metrics"],
        "origin_translations_bohr": translations,
        "rows": rows,
    }
    paths = _write_outputs(records, Path(args.outdir))
    for p in paths:
        print(p)


if __name__ == "__main__":
    main()
