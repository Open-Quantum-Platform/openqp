#!/usr/bin/env python3
"""Profile CPU/GPU METC runtime modes for a fixed polyene size.

This harness is intentionally profiling-only: it reuses the existing polyene
input generator and parser, runs CPU and GPU modes with explicit environment
controls, and records wall time plus OQP_GPU_PROFILE counters emitted by the C/CUDA
runtime.  It does not implement new scientific methods.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import statistics
import sys
from pathlib import Path

from polyene_ladder_metc import parse_nbf, parse_scalars, run_once, write_input

PROFILE_RE = re.compile(r"OQP_GPU_PROFILE\s+(.*)")
KERNEL_PROFILE_RE = re.compile(r"OQP_GPU_KERNEL_PROFILE\s+(.*)")


def parse_keyvals(text: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for part in text.split():
        if "=" in part:
            key, value = part.split("=", 1)
            fields[key] = value
    return fields


def parse_profile(combined: str) -> dict[str, str]:
    rows = []
    for line in combined.splitlines():
        m = PROFILE_RE.search(line)
        if not m:
            continue
        fields = parse_keyvals(m.group(1))
        rows.append(fields)
    if not rows:
        return {}
    out: dict[str, str] = {"profile_records": str(len(rows))}
    numeric_sums: dict[str, float] = {}
    numeric_keys = {
        "h2d_bytes", "d2h_bytes", "h2d_calls", "d2h_calls", "cuda_malloc_calls",
        "cuda_free_calls", "alloc_s", "free_s", "upload_d3_s", "upload_batch_s",
        "zero_f3_s", "download_f3_s", "session_wall_s", "contract_wall_s",
        "contract_calls", "contract_host_overhead_s", "upload_d3_n", "upload_batch_n", "zero_f3_n",
        "download_f3_n", "kernel_launches", "kernel_s", "sync_s",
    }
    for row in rows:
        for key, value in row.items():
            if key in numeric_keys:
                try:
                    numeric_sums[key] = numeric_sums.get(key, 0.0) + float(value)
                except ValueError:
                    pass
            else:
                out.setdefault(key, value)
        if "device_bytes_peak" in row:
            try:
                out["device_bytes_peak"] = str(max(int(out.get("device_bytes_peak", "0")), int(row["device_bytes_peak"])))
            except ValueError:
                out["device_bytes_peak"] = row["device_bytes_peak"]
    for key, value in numeric_sums.items():
        if key.endswith("_bytes") or key.endswith("_calls") or key.endswith("_n") or key == "kernel_launches":
            out[key] = str(int(value))
        else:
            out[key] = f"{value:.9f}"
    return out


def parse_kernel_profiles(combined: str, mode: str, repeat_index: int) -> list[dict[str, str]]:
    rows = []
    session_index = 0
    for line in combined.splitlines():
        if PROFILE_RE.search(line):
            session_index += 1
        m = KERNEL_PROFILE_RE.search(line)
        if not m:
            continue
        row = parse_keyvals(m.group(1))
        row["mode"] = mode
        row["repeat_index"] = str(repeat_index)
        row["session_index"] = str(session_index)
        rows.append(row)
    return rows


def max_scalar_delta(cpu_vals: dict, vals: dict) -> float:
    keys = set(cpu_vals) & set(vals)
    if not keys:
        return float("nan")
    return max(abs(cpu_vals[k] - vals[k]) for k in keys)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-cmd", default="python -m oqp.pyoqp")
    ap.add_argument("--m", type=int, default=6)
    ap.add_argument("--basis", default="6-31g")
    ap.add_argument("--nstate", type=int, default=3)
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--kernel-csv", default=None)
    ap.add_argument("--log-dir", required=True)
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()

    tokens = shlex.split(args.run_cmd)
    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    inp = log_dir / f"polyene_m{args.m}_profile.inp"
    write_input(inp, args.m, args.basis, args.nstate)

    commit = os.popen("git rev-parse --short HEAD").read().strip()
    host = os.uname().nodename
    job = os.environ.get("SLURM_JOB_ID", "")
    gpu_name = os.popen("nvidia-smi --query-gpu=name --format=csv,noheader | head -1").read().strip()

    modes = [
        ("cpu_1t", {"OQP_GPU_METC": "0", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}),
        ("gpu_clean_reference", {
            "OQP_GPU_METC": "1", "OQP_GPU_STRICT": "0", "OQP_GPU_PROFILE": "1",
            "OQP_GPU_METC_VARIANT": "reference",
            "OQP_GPU_F3_CHECK": "0", "OQP_GPU_METC_F3_CHECK": "0",
            "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        }),
        ("gpu_clean_combined_coulomb", {
            "OQP_GPU_METC": "1", "OQP_GPU_STRICT": "0", "OQP_GPU_PROFILE": "1",
            "OQP_GPU_METC_VARIANT": "combined_coulomb",
            "OQP_GPU_F3_CHECK": "0", "OQP_GPU_METC_F3_CHECK": "0",
            "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        }),
        ("gpu_clean_combined_coulomb_warp", {
            "OQP_GPU_METC": "1", "OQP_GPU_STRICT": "0", "OQP_GPU_PROFILE": "1",
            "OQP_GPU_METC_VARIANT": "combined_coulomb_warp",
            "OQP_GPU_F3_CHECK": "0", "OQP_GPU_METC_F3_CHECK": "0",
            "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        }),
        ("gpu_clean_two_phase_accum", {
            "OQP_GPU_METC": "1", "OQP_GPU_STRICT": "0", "OQP_GPU_PROFILE": "1",
            "OQP_GPU_METC_VARIANT": "two_phase_accum",
            "OQP_GPU_F3_CHECK": "0", "OQP_GPU_METC_F3_CHECK": "0",
            "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        }),
        ("gpu_f3check_two_phase_accum", {
            "OQP_GPU_METC": "1", "OQP_GPU_STRICT": "1", "OQP_GPU_PROFILE": "1",
            "OQP_GPU_METC_VARIANT": "two_phase_accum",
            "OQP_GPU_F3_CHECK": "1", "OQP_GPU_METC_F3_CHECK_ATOL": "1e-7",
            "OQP_GPU_METC_F3_CHECK_RTOL": "1e-7", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        }),
        ("gpu_f3check_combined_coulomb_warp", {
            "OQP_GPU_METC": "1", "OQP_GPU_STRICT": "1", "OQP_GPU_PROFILE": "1",
            "OQP_GPU_METC_VARIANT": "combined_coulomb_warp",
            "OQP_GPU_F3_CHECK": "1", "OQP_GPU_METC_F3_CHECK_ATOL": "1e-7",
            "OQP_GPU_METC_F3_CHECK_RTOL": "1e-7", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        }),
        ("gpu_f3check_combined_coulomb", {
            "OQP_GPU_METC": "1", "OQP_GPU_STRICT": "1", "OQP_GPU_PROFILE": "1",
            "OQP_GPU_METC_VARIANT": "combined_coulomb",
            "OQP_GPU_F3_CHECK": "1", "OQP_GPU_METC_F3_CHECK_ATOL": "1e-7",
            "OQP_GPU_METC_F3_CHECK_RTOL": "1e-7", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        }),
    ]

    rows = []
    kernel_rows = []
    cpu_vals = None
    nbf = ""
    for mode, env in modes:
        walls = []
        vals = None
        status = "pass"
        notes = "ok"
        prof: dict[str, str] = {}
        for r in range(args.repeat):
            rc, wall, stderr, log_text, combined = run_once(
                tokens, inp, env, log_dir / f"{mode}_r{r}.log"
            )
            if rc != 0:
                status = "fail"
                notes = f"rc={rc}"
                break
            walls.append(wall)
            vals = parse_scalars(log_text)
            nbf = nbf or parse_nbf(log_text)
            parsed = parse_profile(combined)
            kernel_rows.extend(parse_kernel_profiles(combined, mode, r))
            for key, value in parsed.items():
                prof[key] = value
        if mode == "cpu_1t" and vals:
            cpu_vals = vals
        delta = max_scalar_delta(cpu_vals or vals or {}, vals or {}) if vals else float("nan")
        med = statistics.median(walls) if walls else float("nan")
        row = {
            "mode": mode, "variant": env.get("OQP_GPU_METC_VARIANT", ""),
            "m": args.m, "formula": f"C{2*args.m}H{2*args.m+2}",
            "basis": args.basis, "nstate": args.nstate, "nbf": nbf,
            "median_wall_s": f"{med:.6f}" if walls else "", "repeat": len(walls),
            "scalar_max_abs_delta_vs_cpu_1t": f"{delta:.12g}" if vals else "",
            "status": status, "notes": notes, "gpu_name": gpu_name, "host": host,
            "slurm_job_id": job, "commit": commit,
        }
        row.update(prof)
        rows.append(row)

    fields = [
        "mode", "variant", "m", "formula", "basis", "nstate", "nbf", "median_wall_s", "repeat",
        "scalar_max_abs_delta_vs_cpu_1t", "profile_records", "target", "nf", "nmatrix",
        "nthreads", "max_ncur", "total_bytes", "device_bytes_peak", "h2d_bytes",
        "d2h_bytes", "h2d_calls", "d2h_calls", "cuda_malloc_calls", "cuda_free_calls",
        "alloc_s", "free_s", "upload_d3_s", "upload_batch_s", "zero_f3_s",
        "download_f3_s", "session_wall_s", "contract_wall_s", "contract_calls",
        "contract_host_overhead_s", "upload_d3_n", "upload_batch_n", "zero_f3_n",
        "download_f3_n", "kernel_launches", "kernel_s", "sync_s", "f3_check",
        "timing_semantics", "status", "notes", "gpu_name", "host", "slurm_job_id", "commit",
    ]
    with open(args.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    kernel_csv = args.kernel_csv or str(Path(args.out_csv).with_suffix(".kernels.csv"))
    kernel_fields = [
        "mode", "repeat_index", "session_index", "schema", "kernel", "call_count",
        "total_s", "mean_s", "median_s", "p95_s", "min_s", "max_s", "block_dim",
        "grid_blocks_mean", "grid_blocks_max", "ncur_min", "ncur_max", "elements_total",
        "elements_mean", "elements_min", "elements_max", "est_flops", "est_bytes",
        "est_atomics", "est_arith_intensity", "sample_count", "timing_note",
    ]
    with open(kernel_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=kernel_fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(kernel_rows)
    Path(args.summary).write_text("\n".join([
        f"GPU METC profile breakdown job {job}",
        f"commit: {commit}",
        f"polyene m={args.m} basis={args.basis} nstate={args.nstate}",
        "modes: CPU 1t, GPU clean reference, GPU clean combined-coulomb, GPU clean combined-coulomb-warp, GPU clean two-phase-accum, f3-check variants",
        f"csv: {args.out_csv}",
        f"kernel_csv: {kernel_csv}",
        f"log_dir: {args.log_dir}",
    ]) + "\n")
    print(Path(args.summary).read_text())
    return 0 if all(r["status"] == "pass" for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
