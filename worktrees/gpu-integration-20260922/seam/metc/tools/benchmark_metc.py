#!/usr/bin/env python3
"""C1c benchmark harness for the MRSF/UMRSF METC path (preparation only).

This script is the *scaffolding* for the eventual GPU performance benchmark.
It measures end-to-end wall time for a CPU run and a GPU-METC run, and parses
the runtime's separated-region timing line (emitted to stderr when
``OQP_GPU_METC_TIMING=1``) so the device-transfer regions (upload d3 / zero f3
/ download f3) can be reported separately from kernel/host time.

HONESTY CONTRACT
----------------
This harness refuses to print a speedup unless it actually measured a GPU run
on a machine with a CUDA toolchain.  Running it here (no GPU) yields CPU-only
numbers labelled as such.  A speedup printed by this tool always corresponds to
a real measured GPU run -- it never extrapolates or fabricates.

C1c performance numbers are NOT valid until this is run on real Fortran/CUDA
hardware against a deterministic, parity-validated case (see
``validate_metc_cuda.py`` and ``docs/gpu_metc_validation.md``).
"""

import argparse
import os
import re
import shutil
import statistics
import subprocess
import sys
import time


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_EXAMPLE = os.path.join("examples", "MRSF-TDDFT", "H2O_BHHLYP-MRSFTDDFT_ENERGY.json")

TIMING_RE = re.compile(
    r"OQP_GPU_METC_TIMING\s+"
    r"upload_d3_s=(?P<upload>[-+0-9.eE]+)\s+"
    r"zero_f3_s=(?P<zero>[-+0-9.eE]+)\s+"
    r"download_f3_s=(?P<download>[-+0-9.eE]+)")


def log(msg):
    print("[benchmark-metc] " + msg, flush=True)


def have_cuda_toolchain():
    return shutil.which("nvcc") is not None


def run_once(run_tokens, example, use_gpu, timing):
    env = dict(os.environ)
    env["OQP_GPU_METC"] = "1" if use_gpu else "0"
    if use_gpu and timing:
        env["OQP_GPU_METC_TIMING"] = "1"
    t0 = time.perf_counter()
    proc = subprocess.run(list(run_tokens) + [example], env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True)
    wall = time.perf_counter() - t0
    return proc.returncode, wall, proc.stderr


def parse_region_timing(stderr_text):
    regions = {"upload": 0.0, "zero": 0.0, "download": 0.0}
    found = False
    for m in TIMING_RE.finditer(stderr_text):
        found = True
        regions["upload"] += float(m.group("upload"))
        regions["zero"] += float(m.group("zero"))
        regions["download"] += float(m.group("download"))
    return (regions if found else None)


def bench(run_tokens, example, use_gpu, repeat, warmup, timing):
    for _ in range(warmup):
        rc, _, _ = run_once(run_tokens, example, use_gpu, timing)
        if rc != 0:
            return None, None, rc
    walls = []
    last_regions = None
    for _ in range(repeat):
        rc, wall, stderr = run_once(run_tokens, example, use_gpu, timing)
        if rc != 0:
            return None, None, rc
        walls.append(wall)
        if use_gpu and timing:
            r = parse_region_timing(stderr)
            if r is not None:
                last_regions = r
    return walls, last_regions, 0


def summarize(walls):
    return {
        "min": min(walls),
        "median": statistics.median(walls),
        "mean": statistics.mean(walls),
        "n": len(walls),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--example", default=os.path.join(REPO_ROOT, DEFAULT_EXAMPLE))
    p.add_argument("--run-cmd", default="openqp")
    p.add_argument("--repeat", type=int, default=3)
    p.add_argument("--warmup", type=int, default=1)
    p.add_argument("--require-gpu", action="store_true",
                   help="exit non-zero if no CUDA toolchain is available")
    args = p.parse_args(argv)

    run_tokens = args.run_cmd.split()
    gpu_available = have_cuda_toolchain()

    if not os.path.exists(args.example):
        log("ERROR: example not found: %s" % args.example)
        return 1

    log("CPU baseline (OQP_GPU_METC=0): repeat=%d warmup=%d" % (args.repeat, args.warmup))
    cpu_walls, _, rc = bench(run_tokens, args.example, False, args.repeat, args.warmup, False)
    if rc != 0:
        log("ERROR: CPU run exited %d" % rc)
        return 1
    cpu = summarize(cpu_walls)
    log("CPU wall (s): min=%.4f median=%.4f mean=%.4f (n=%d)"
        % (cpu["min"], cpu["median"], cpu["mean"], cpu["n"]))

    if not gpu_available:
        log("GPU UNAVAILABLE: no CUDA toolchain (nvcc) found.")
        log("Reporting CPU-only numbers. NO speedup is claimed.")
        return 1 if args.require_gpu else 0

    log("GPU METC (OQP_GPU_METC=1, OQP_GPU_METC_TIMING=1) ...")
    gpu_walls, regions, rc = bench(run_tokens, args.example, True, args.repeat, args.warmup, True)
    if rc != 0:
        log("ERROR: GPU run exited %d (an all-or-nothing abort on device "
            "failure is a deliberate, correct outcome)." % rc)
        return 1
    gpu = summarize(gpu_walls)
    log("GPU wall (s): min=%.4f median=%.4f mean=%.4f (n=%d)"
        % (gpu["min"], gpu["median"], gpu["mean"], gpu["n"]))

    if regions is not None:
        log("Separated device-transfer regions (s, summed over the run):")
        log("  upload_d3=%.6f  zero_f3=%.6f  download_f3=%.6f"
            % (regions["upload"], regions["zero"], regions["download"]))
    else:
        log("NOTE: no OQP_GPU_METC_TIMING lines were emitted; region breakdown "
            "unavailable (was the CUDA build run?).")

    # Speedup is only printed because we actually measured a real GPU run.
    speedup = cpu["median"] / gpu["median"] if gpu["median"] > 0 else float("nan")
    log("Measured median speedup (CPU/GPU): %.2fx" % speedup)
    return 0


if __name__ == "__main__":
    sys.exit(main())
