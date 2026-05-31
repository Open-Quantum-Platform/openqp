#!/usr/bin/env python3
"""Build and populate the closed-shell OpenQP/ddX PCM benchmark matrix.

The references produced here are *trusted regression baselines* for the stabilized
OpenQP/ddX closed-shell diagnostic convention, not literature values. They are
intended to turn pending-reference skips into hard pass/fail regression tests for
this branch.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data" / "pcm_literature_benchmarks.json"
TESTS = ROOT / "tests"

EPS_WATER = 78.3553
BASIS = "6-31g*"
TOL = 1.0e-7

MOLECULES = {
    "h2o": {
        "name": "H2O",
        "source": "OpenQP examples/HF/H2O geometry",
        "system": [
            "8   0.000000000   0.000000000  -0.041061554",
            "1  -0.533194329   0.533194329  -0.614469223",
            "1   0.533194329  -0.533194329  -0.614469223",
        ],
    },
    "nh3": {
        "name": "NH3",
        "source": "standard C3v equilibrium geometry",
        "system": [
            "7   0.000000000   0.000000000   0.116489000",
            "1   0.000000000   0.939731000  -0.271808000",
            "1   0.813831000  -0.469865000  -0.271808000",
            "1  -0.813831000  -0.469865000  -0.271808000",
        ],
    },
    "hf": {
        "name": "HF",
        "source": "linear HF, r(HF)=0.917 Angstrom",
        "system": [
            "9   0.000000000   0.000000000   0.000000000",
            "1   0.000000000   0.000000000   0.917000000",
        ],
    },
    "ch4": {
        "name": "CH4",
        "source": "tetrahedral methane, r(CH)=1.09 Angstrom",
        "system": [
            "6   0.000000000   0.000000000   0.000000000",
            "1   0.629311793   0.629311793   0.629311793",
            "1  -0.629311793  -0.629311793   0.629311793",
            "1  -0.629311793   0.629311793  -0.629311793",
            "1   0.629311793  -0.629311793  -0.629311793",
        ],
    },
    "co": {
        "name": "CO",
        "source": "linear CO, r(CO)=1.128 Angstrom",
        "system": [
            "6   0.000000000   0.000000000   0.000000000",
            "8   0.000000000   0.000000000   1.128000000",
        ],
    },
    "co2": {
        "name": "CO2",
        "source": "linear OCO, r(CO)=1.160 Angstrom",
        "system": [
            "8   0.000000000   0.000000000  -1.160000000",
            "6   0.000000000   0.000000000   0.000000000",
            "8   0.000000000   0.000000000   1.160000000",
        ],
    },
    "hcn": {
        "name": "HCN",
        "source": "linear HCN, r(CH)=1.064 Angstrom, r(CN)=1.156 Angstrom",
        "system": [
            "1   0.000000000   0.000000000  -1.064000000",
            "6   0.000000000   0.000000000   0.000000000",
            "7   0.000000000   0.000000000   1.156000000",
        ],
    },
    "h2co": {
        "name": "H2CO",
        "source": "planar formaldehyde, r(CO)=1.210 Angstrom, r(CH)=1.100 Angstrom, angle H-C-H about 116.5 deg",
        "system": [
            "6   0.000000000   0.000000000   0.000000000",
            "8   0.000000000   0.000000000   1.210000000",
            "1   0.941985000   0.000000000  -0.568834000",
            "1  -0.941985000   0.000000000  -0.568834000",
        ],
    },
    "ch3oh": {
        "name": "CH3OH",
        "source": "staggered methanol approximate equilibrium geometry",
        "system": [
            "6   0.000000000   0.000000000   0.000000000",
            "8   1.430000000   0.000000000   0.000000000",
            "1   1.750000000   0.900000000   0.000000000",
            "1  -0.360000000   1.020000000   0.000000000",
            "1  -0.360000000  -0.510000000   0.883346000",
            "1  -0.360000000  -0.510000000  -0.883346000",
        ],
    },
    "ch3cn": {
        "name": "CH3CN",
        "source": "linear C-C-N skeleton acetonitrile approximate geometry",
        "system": [
            "6   0.000000000   0.000000000   0.000000000",
            "6   0.000000000   0.000000000   1.460000000",
            "7   0.000000000   0.000000000   2.620000000",
            "1   1.029000000   0.000000000  -0.363000000",
            "1  -0.514500000   0.891109000  -0.363000000",
            "1  -0.514500000  -0.891109000  -0.363000000",
        ],
    },
}

METHODS = [
    ("hf", None, "RHF/HF"),
    ("hf", "bhhlyp", "RHF/BHHLYP"),
    ("hf", "pbe", "RHF/PBE"),
]


def git_text(args: list[str]) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()
    except Exception:
        return "unknown"


def benchmark_rows() -> list[dict]:
    rows = []
    for mol_id, mol in MOLECULES.items():
        for method, functional, label in METHODS:
            suffix = functional if functional else "hf"
            row = {
                "id": f"{mol_id}_rhf_{suffix}_631gs_water",
                "category": "trusted_reference_regression",
                "tier": "tier2_closed_shell_qm_scf_pcm_10mol",
                "molecule": mol["name"],
                "geometry_source": mol["source"],
                "system": mol["system"],
                "charge": 0,
                "multiplicity": 1,
                "scf_type": "rhf",
                "method": method,
                "basis": BASIS,
                "epsilon": EPS_WATER,
                "quantity": "e_pcm",
                "reference_value": None,
                "tolerance": TOL,
                "reference_kind": "openqp_ddx_closed_shell_diagnostic_baseline",
                "reference_source": "pending population from stabilized OpenQP/ddX closed-shell diagnostic baseline",
                "doi_or_url": "n/a; OpenQP/ddX same-branch trusted regression baseline, not a literature value",
                "status": "pending_reference",
                "benchmark_label": f"{mol['name']} {label}/{BASIS} ddPCM water",
            }
            if functional:
                row["functional"] = functional
            rows.append(row)
    return rows


def reset_table() -> None:
    data = json.loads(DATA.read_text(encoding="utf-8"))
    data["_about"] = (
        "Reference targets for ddX-enabled PCM validation. PRODUCTION PCM PHYSICS IS FORTRAN-SIDE; "
        "these values are comparison targets only and never feed the SCF. The tier2 closed-shell rows "
        "are OpenQP/ddX trusted regression baselines from the stabilized diagnostic convention, not "
        "literature benchmarks. A row is a hard pass/fail gate when reference_value is non-null and "
        "status is 'verified'."
    )
    data["benchmarks"] = benchmark_rows()
    DATA.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def populate_table() -> None:
    sys.path.insert(0, str(TESTS))
    import test_pcm_literature_benchmarks as benchmod  # noqa: PLC0415

    data = json.loads(DATA.read_text(encoding="utf-8"))
    by_id = {row["id"]: row for row in data["benchmarks"]}
    branch = git_text(["branch", "--show-current"])
    commit = git_text(["rev-parse", "--short=12", "HEAD"])
    for row in data["benchmarks"]:
        proc_vac, log_vac = benchmod._run(benchmod._input_text(row, pcm_on=False))
        if proc_vac.returncode != 0:
            raise RuntimeError(f"vacuum run failed for {row['id']}\n{log_vac[-2000:]}")
        e_vac = benchmod._total_energy(log_vac)
        if e_vac is None or not math.isfinite(e_vac):
            raise RuntimeError(f"missing vacuum energy for {row['id']}")

        proc_pcm, log_pcm = benchmod._run(benchmod._input_text(row, pcm_on=True))
        if proc_pcm.returncode != 0:
            raise RuntimeError(f"PCM run failed for {row['id']}\n{log_pcm[-3000:]}")
        e_pcm = benchmod._pcm_energy(log_pcm)
        diag = benchmod._pcm_diag(log_pcm)
        if e_pcm is None or not math.isfinite(e_pcm):
            raise RuntimeError(f"missing PCM energy for {row['id']}")
        for key in ("fock_q_scale", "fd_fock_scale_mean", "phi_source_vs_exact_rms"):
            if key not in diag:
                raise RuntimeError(f"missing diagnostic {key} for {row['id']}")

        row["vacuum_total_energy"] = e_vac
        row["reference_value"] = e_pcm
        row["phi_source_vs_exact_rms_threshold"] = max(0.05, round(abs(diag.get("phi_source_vs_exact_rms", 0.0)) * 1.25 + 1.0e-12, 12))
        row["phi_source_vs_exact_max_threshold"] = max(0.10, round(abs(diag.get("phi_source_vs_exact_max", 0.0)) * 1.25 + 1.0e-12, 12))
        row["reference_source"] = (
            f"OpenQP/ddX closed-shell diagnostic baseline generated on {branch} at {commit}; "
            "same Fortran/C ddX protocol as the regression test. This is a trusted regression "
            "baseline, not an independent literature value."
        )
        row["status"] = "verified"
        row["empirical_verification"] = {
            "generator": "scripts/pcm_benchmark_matrix.py --populate",
            "openqp_branch": branch,
            "openqp_commit": commit,
            "pcm_total_energy_hartree": benchmod._total_energy(log_pcm),
            "pcm_solvent_energy_hartree": e_pcm,
            "diagnostics": diag,
        }
    DATA.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reset", action="store_true", help="replace benchmark rows with the 10 molecule x 3 method matrix")
    parser.add_argument("--populate", action="store_true", help="run OpenQP vacuum/PCM jobs and populate verified regression references")
    args = parser.parse_args()
    if not args.reset and not args.populate:
        parser.error("choose --reset and/or --populate")
    if args.reset:
        reset_table()
    if args.populate:
        populate_table()
    print(json.dumps({"benchmarks": len(benchmark_rows()), "data": str(DATA)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
