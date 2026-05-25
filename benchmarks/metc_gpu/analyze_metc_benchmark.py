#!/usr/bin/env python3
"""Generate METC GPU manuscript tables and figures from benchmark CSV.

Input CSV columns are produced by benchmarks/metc_gpu/metc_cuda_benchmark.cu.
This script writes:
  - processed summary CSV
  - LaTeX table fragments
  - PDF figures for Overleaf/manuscript use
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.csv)
    df["case"] = df.apply(lambda r: f"{r.mode.upper()} nbf={int(r.nbf)}, nf={int(r.nf)}, pass={int(r['pass'])}", axis=1)
    df["passed"] = (df["ierr"] == 0) & (df["max_abs_diff"] < 1e-9)
    df.to_csv(args.outdir / "metc_gpu_benchmark_summary.csv", index=False)

    validation_cols = ["mode", "nbf", "nf", "nmatrix", "ncur", "pass", "max_abs_diff", "rms_diff", "rel_l2_diff", "passed"]
    validation = df[validation_cols].copy()
    validation["mode"] = validation["mode"].str.upper()
    latex = validation.to_latex(index=False, escape=False, float_format=lambda x: f"{x:.3e}")
    (args.outdir / "table_metc_validation.tex").write_text(latex)

    perf_cols = ["mode", "nbf", "nf", "ncur", "pass", "cpu_ms", "gpu_total_ms", "speedup_total"]
    perf = df[perf_cols].copy()
    perf["mode"] = perf["mode"].str.upper()
    latex = perf.to_latex(index=False, escape=False, float_format=lambda x: f"{x:.3f}")
    (args.outdir / "table_metc_performance.tex").write_text(latex)

    plt.rcParams.update({"font.size": 9, "figure.dpi": 160})
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    labels = [f"{m.upper()}\n{nbf}/{nf}" for m, nbf, nf in zip(df["mode"], df["nbf"], df["nf"])]
    x = range(len(df))
    ax.bar([i - 0.18 for i in x], df["cpu_ms"], width=0.36, label="CPU reference")
    ax.bar([i + 0.18 for i in x], df["gpu_total_ms"], width=0.36, label="CUDA total")
    ax.set_yscale("log")
    ax.set_ylabel("Time / ms (log scale)")
    ax.set_xticks(list(x), labels, rotation=45, ha="right")
    ax.set_title("Isolated METC contraction timing")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.outdir / "fig_metc_timing.pdf")
    fig.savefig(args.outdir / "fig_metc_timing.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 3.4))
    colors = ["#1f77b4" if p else "#d62728" for p in df["passed"]]
    ax.bar(labels, df["max_abs_diff"], color=colors)
    ax.axhline(1e-9, color="black", linestyle="--", linewidth=1, label="1e-9 target")
    ax.set_yscale("log")
    ax.set_ylabel("Max |CPU - CUDA|")
    ax.set_title("METC numerical reproducibility")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.outdir / "fig_metc_error.pdf")
    fig.savefig(args.outdir / "fig_metc_error.png", dpi=300)
    plt.close(fig)

    with (args.outdir / "metc_gpu_benchmark_report.md").open("w") as f:
        f.write("# METC GPU benchmark report\n\n")
        f.write(f"Input CSV: `{args.csv}`\n\n")
        f.write(f"Cases: {len(df)}\n\n")
        f.write(f"Passed 1e-9 max-abs criterion: {int(df['passed'].sum())}/{len(df)}\n\n")
        f.write(f"Worst max abs diff: {df['max_abs_diff'].max():.6e}\n\n")
        f.write(f"Best total speedup: {df['speedup_total'].max():.3f}\n\n")
        f.write("Note: CUDA total includes the current prototype's per-call cudaMalloc and full tensor transfers. ")
        f.write("Kernel-resident/persistent-buffer timings should be added before final performance claims.\n")


if __name__ == "__main__":
    main()
