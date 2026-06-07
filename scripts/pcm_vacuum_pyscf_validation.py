#!/usr/bin/env python3
"""Validate OpenQP vacuum energies against PySCF for the PCM benchmark matrix."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

from pyscf import dft, gto, scf

ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / "tests"
DATA = ROOT / "tests" / "data" / "pcm_literature_benchmarks.json"

ATOM_MAP = {
    "1": "H",
    "6": "C",
    "7": "N",
    "8": "O",
    "9": "F",
}
XC_MAP = {
    "pbe": "pbe,pbe",
    # OpenQP's BHHLYP label is Becke half-and-half LYP; PySCF accepts bhandhlyp.
    "bhhlyp": "bhandhlyp",
}


def pyscf_geom(system: list[str]) -> str:
    lines = []
    for line in system:
        parts = line.split()
        sym = ATOM_MAP.get(parts[0], parts[0])
        lines.append(" ".join([sym, *parts[1:4]]))
    return "; ".join(lines)


def run_pyscf(row: dict) -> float:
    mol = gto.M(
        atom=pyscf_geom(row["system"]),
        basis=row["basis"],
        charge=row["charge"],
        spin=row["multiplicity"] - 1,
        unit="Angstrom",
        cart=True,
        verbose=0,
    )
    functional = row.get("functional")
    if functional:
        mf = dft.RKS(mol)
        mf.xc = XC_MAP.get(functional, functional)
        mf.grids.atom_grid = (99, 590)
    else:
        mf = scf.RHF(mol)
    mf.conv_tol = 1.0e-11
    mf.max_cycle = 100
    return float(mf.kernel())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", type=Path, default=ROOT / "tests" / "data" / "pcm_vacuum_pyscf_validation.json")
    parser.add_argument("--hf-tol", type=float, default=2.0e-7)
    parser.add_argument("--dft-tol", type=float, default=1.0e-4)
    args = parser.parse_args()

    sys.path.insert(0, str(TESTS))
    import test_pcm_literature_benchmarks as benchmod  # noqa: PLC0415

    data = json.loads(DATA.read_text(encoding="utf-8"))
    records = []
    failed = []
    for row in data["benchmarks"]:
        proc, log = benchmod._run(benchmod._input_text(row, pcm_on=False))
        if proc.returncode != 0:
            raise RuntimeError(f"OpenQP vacuum failed for {row['id']}\n{log[-2000:]}")
        e_oqp = benchmod._total_energy(log)
        if e_oqp is None or not math.isfinite(e_oqp):
            raise RuntimeError(f"OpenQP vacuum missing for {row['id']}")
        e_pyscf = run_pyscf(row)
        delta = abs(e_oqp - e_pyscf)
        tol = args.dft_tol if row.get("functional") else args.hf_tol
        status = "pass" if delta < tol else "fail"
        rec = {
            "id": row["id"],
            "functional": row.get("functional"),
            "openqp_vacuum_total_energy_hartree": e_oqp,
            "pyscf_vacuum_total_energy_hartree": e_pyscf,
            "abs_delta_hartree": delta,
            "tolerance_hartree": tol,
            "status": status,
        }
        records.append(rec)
        if status != "pass":
            failed.append(rec)
        print(f"{status:4s} {row['id']:32s} delta={delta:.12g}")

    out = {
        "schema": "openqp-pyscf-vacuum-validation-v1",
        "note": "OpenQP vs PySCF vacuum validation for PCM benchmark geometries. HF rows are the strict cross-code gate; DFT rows are protocol-near because OpenQP LibXC+SG1/MHL grid settings are not bit-identical to the PySCF grid used here.",
        "records": records,
        "summary": {"total": len(records), "failed": len(failed), "max_delta_hartree": max(r["abs_delta_hartree"] for r in records)},
    }
    args.json_out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if failed:
        print(json.dumps({"failed": failed}, indent=2))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
