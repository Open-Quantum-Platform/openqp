#!/usr/bin/env python3
"""C1c CPU/GPU parity validation harness for the MRSF/UMRSF METC path.

This is a *correctness* harness, not a performance harness.  On a real
Fortran/CUDA machine it:

  1. (optionally) builds OpenQP with ``-DENABLE_CUDA=ON``;
  2. runs a small deterministic MRSF/UMRSF case twice -- once on the CPU
     (``OQP_GPU_METC=0``) and once on the GPU METC path (``OQP_GPU_METC=1``);
  3. compares the scalar results within a numeric tolerance.

It deliberately makes NO speedup claim and NEVER fabricates results.  If a CUDA
toolchain is unavailable it exits with status 2 (SKIPPED) rather than pretending
the GPU path was validated.

Exit codes:
    0  parity validated within tolerance
    1  parity FAILED (mismatch) or a run/build error
    2  SKIPPED (no CUDA toolchain / GPU path not available)

The output-parsing step is intentionally simple and may need to be adapted to
the exact OpenQP result format on the target machine; see ``extract_scalars``.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_EXAMPLE = os.path.join("examples", "MRSF-TDDFT", "H2O_BHHLYP-MRSFTDDFT_ENERGY.json")

EXIT_OK = 0
EXIT_FAIL = 1
EXIT_SKIP = 2

# Lines whose trailing float we treat as a comparable scalar result.  Tuned to
# be conservative; extend on the target machine to cover the quantities you care
# about (e.g. excitation energies, oscillator strengths).
SCALAR_KEYWORDS = (
    "total energy",
    "excitation energy",
    "final energy",
    "scf energy",
)

FLOAT_RE = re.compile(r"[-+]?\d+\.\d+(?:[eEdD][-+]?\d+)?")


def log(msg):
    print("[validate-metc] " + msg, flush=True)


def have_cuda_toolchain():
    return shutil.which("nvcc") is not None


def run_cmd(cmd, env=None, cwd=None):
    log("$ " + " ".join(cmd))
    proc = subprocess.run(cmd, env=env, cwd=cwd, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, universal_newlines=True)
    return proc.returncode, proc.stdout


def build(source_dir, build_dir):
    rc, out = run_cmd(["cmake", "-S", source_dir, "-B", build_dir,
                       "-DENABLE_CUDA=ON"])
    if rc != 0:
        log("cmake configure failed:\n" + out)
        return False
    rc, out = run_cmd(["cmake", "--build", build_dir, "-j"])
    if rc != 0:
        log("cmake build failed:\n" + out)
        return False
    return True


def extract_scalars(text):
    """Extract a deterministic list of (label, value) scalar results.

    Returns floats from lines that mention a known result keyword.  ``D``-style
    Fortran exponents are normalized so ``1.23D-04`` parses as a Python float.
    """
    scalars = []
    for line in text.splitlines():
        low = line.lower()
        if not any(k in low for k in SCALAR_KEYWORDS):
            continue
        matches = FLOAT_RE.findall(line)
        if not matches:
            continue
        value = float(matches[-1].replace("D", "E").replace("d", "e"))
        scalars.append((low.strip(), value))
    return scalars


def run_case(run_cmd_tokens, example, use_gpu, strict):
    env = dict(os.environ)
    env["OQP_GPU_METC"] = "1" if use_gpu else "0"
    if use_gpu and strict:
        env["OQP_GPU_METC_STRICT"] = "1"
    rc, out = run_cmd(list(run_cmd_tokens) + [example], env=env)

    # PyOQP writes the chemically relevant result table to <input>.log rather
    # than stdout. Append that log immediately after each run so CPU/GPU
    # comparisons are made from the corresponding run, even though the GPU run
    # later overwrites the same log path.
    root, _ = os.path.splitext(example)
    log_path = root + ".log"
    if os.path.exists(log_path):
        try:
            with open(log_path, "r", encoding="utf-8", errors="replace") as fh:
                out += "\n[validate-metc] Captured PyOQP log: %s\n" % log_path
                out += fh.read()
        except OSError as exc:
            out += "\n[validate-metc] WARNING: could not read PyOQP log %s: %s\n" % (log_path, exc)
    return rc, out


def compare(cpu, gpu, tol):
    if not cpu or not gpu:
        log("ERROR: no comparable scalars extracted from one or both runs.")
        log("       Adapt extract_scalars() to the OpenQP output format.")
        return False
    if len(cpu) != len(gpu):
        log("ERROR: CPU produced %d scalars but GPU produced %d." % (len(cpu), len(gpu)))
        return False
    ok = True
    for (lc, vc), (lg, vg) in zip(cpu, gpu):
        diff = abs(vc - vg)
        scale = max(1.0, abs(vc))
        rel = diff / scale
        status = "ok" if rel <= tol else "MISMATCH"
        if rel > tol:
            ok = False
        log("  %-40s cpu=%.10g gpu=%.10g |d|=%.3g rel=%.3g  %s"
            % (lc[:40], vc, vg, diff, rel, status))
    return ok


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source-dir", default=REPO_ROOT)
    p.add_argument("--build-dir", default=os.path.join(REPO_ROOT, "build-cuda"))
    p.add_argument("--example", default=os.path.join(REPO_ROOT, DEFAULT_EXAMPLE))
    p.add_argument("--run-cmd", default="openqp",
                   help="command used to run an input (default: openqp)")
    p.add_argument("--tol", type=float, default=1e-7,
                   help="relative tolerance for CPU/GPU parity (default 1e-7)")
    p.add_argument("--strict", action="store_true",
                   help="run the GPU case with OQP_GPU_METC_STRICT=1")
    p.add_argument("--skip-build", action="store_true")
    args = p.parse_args(argv)

    if not have_cuda_toolchain():
        log("SKIPPED: no CUDA toolchain (nvcc) found. The GPU METC path cannot "
            "be built or validated here. This is NOT a pass.")
        return EXIT_SKIP

    if not args.skip_build:
        if not build(args.source_dir, args.build_dir):
            log("FAILED: build error.")
            return EXIT_FAIL

    if not os.path.exists(args.example):
        log("FAILED: example not found: %s" % args.example)
        return EXIT_FAIL

    run_tokens = args.run_cmd.split()

    log("Running CPU reference (OQP_GPU_METC=0) ...")
    rc_cpu, out_cpu = run_case(run_tokens, args.example, use_gpu=False, strict=False)
    if rc_cpu != 0:
        log("FAILED: CPU run exited %d:\n%s" % (rc_cpu, out_cpu))
        return EXIT_FAIL

    log("Running GPU METC (OQP_GPU_METC=1) ...")
    rc_gpu, out_gpu = run_case(run_tokens, args.example, use_gpu=True, strict=args.strict)
    if rc_gpu != 0:
        log("FAILED: GPU run exited %d (note: an all-or-nothing abort is a "
            "deliberate, correct outcome on device failure):\n%s" % (rc_gpu, out_gpu))
        return EXIT_FAIL

    cpu = extract_scalars(out_cpu)
    gpu = extract_scalars(out_gpu)
    log("Comparing %d CPU vs %d GPU scalar(s) at rel tol %.1e ..."
        % (len(cpu), len(gpu), args.tol))
    if compare(cpu, gpu, args.tol):
        log("PASS: CPU/GPU MRSF METC parity within tolerance. "
            "(Correctness only -- no performance is claimed.)")
        return EXIT_OK
    log("FAIL: CPU/GPU parity mismatch.")
    return EXIT_FAIL


if __name__ == "__main__":
    sys.exit(main())
