#!/usr/bin/env python3
"""Validate OpenQP/ddX PCM e_pcm against independent PySCF ddPCM for all rows.

Runs each closed-shell RHF/HF benchmark through (a) the real OpenQP/ddX path and
(b) PySCF native ddPCM at matched protocol (ddPCM, eps=78.3553, lmax=8, Lebedev
302, eta=0.1, Bondi-by-Z radii), and compares the reaction-field energy:
OpenQP e_pcm  vs  PySCF with_solvent.e.

Reference/diagnostic only; production PCM stays Fortran-side. Requires a
ddX-enabled OpenQP build (lib/liboqp + LD_LIBRARY_PATH to libddx) and pyscf.
"""
from __future__ import annotations
import json, os, re, subprocess, sys, tempfile
from pathlib import Path
import numpy as np
from pyscf import gto, scf
from pyscf.solvent import ddPCM

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data" / "pcm_literature_benchmarks.json"
BOHR = 1.8897259886
BONDI_A = {1: 1.20, 6: 1.70, 7: 1.55, 8: 1.52, 9: 1.47}
Z2SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F"}
K = 627.5094740631


def radii_table(zmax=36):
    t = np.full(zmax + 1, 1.80 * BOHR); t[0] = 1.50 * BOHR
    for z, a in BONDI_A.items():
        t[z] = a * BOHR
    return t


def to_pyscf_geom(system):
    lines = []
    for ln in system:
        p = ln.split()
        lines.append(f"{Z2SYM[int(p[0])]} {p[1]} {p[2]} {p[3]}")
    return "\n".join(lines)


def pyscf_with_solvent_e(system):
    mol = gto.M(atom=to_pyscf_geom(system), unit="Angstrom", basis="6-31g*",
                charge=0, spin=0, cart=False, verbose=0)
    mf = ddPCM(scf.RHF(mol)); ws = mf.with_solvent
    ws.eps = 78.3553; ws.lmax = 8; ws.lebedev_order = 29; ws.eta = 0.1
    ws.radii_table = radii_table()
    mf.conv_tol = 1e-11
    mf.kernel()
    return float(ws.e)


def openqp_e_pcm(system):
    inp = ("[input]\nsystem=\n" + "\n".join("   " + l for l in system) +
           "\ncharge=0\nruntype=energy\nbasis=6-31g*\nmethod=hf\n\n[guess]\ntype=huckel\n\n"
           "[scf]\nmultiplicity=1\ntype=rhf\n\n[pcm]\nenabled=true\nbackend=ddx\n"
           "mode=reference_scf\nmodel=ddpcm\nepsilon=78.3553\n")
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "m.inp"; p.write_text(inp)
        env = dict(os.environ)
        env["OPENQP_ROOT"] = str(ROOT)
        env["PYTHONPATH"] = str(ROOT / "pyoqp") + os.pathsep + env.get("PYTHONPATH", "")
        env["OMP_NUM_THREADS"] = "1"
        subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(p)], env=env,
                       capture_output=True, text=True)
        log = (p.with_suffix(".log")).read_text()
    m = re.findall(r"PCM solvent energy.*?=\s*(-?\d+\.\d+)", log)
    return float(m[-1]) if m else None


def main():
    data = json.loads(DATA.read_text())
    seen = set()
    print(f"{'molecule':10} {'OpenQP e_pcm':>13} {'PySCF wse':>11} {'diff(kcal)':>10} {'ratio':>6}")
    for b in data["benchmarks"]:
        if b["id"].endswith("_hf_631gs_water") and "_rhf_hf_" in b["id"]:
            mol = b["molecule"]
            if mol in seen:
                continue
            seen.add(mol)
            oq = openqp_e_pcm(b["system"])
            ps = pyscf_with_solvent_e(b["system"])
            diff = (oq - ps) * K if oq is not None else float("nan")
            ratio = oq / ps if oq is not None else float("nan")
            print(f"{mol:10} {oq:>13.6f} {ps:>11.6f} {diff:>10.2f} {ratio:>6.3f}")


if __name__ == "__main__":
    main()
