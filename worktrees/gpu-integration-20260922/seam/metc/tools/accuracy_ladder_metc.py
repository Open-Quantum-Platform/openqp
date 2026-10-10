#!/usr/bin/env python3
"""Run CPU/GPU METC accuracy ladder and emit raw CSV evidence."""
import argparse, csv, os, re, shutil, subprocess, sys
from pathlib import Path

FLOAT = r"[-+]?\d+\.\d+(?:[eEdD][-+]?\d+)?"
TOTAL_RE = re.compile(r"TOTAL energy\s*=\s*(%s)" % FLOAT, re.I)
NBF_RE = re.compile(r"Number of Basis Set functions\s*=\s*(\d+)", re.I)
SUMMARY_RE = re.compile(
    r"^\s*(\d+)\s+(%s)\s+(%s)\s+(%s)\s+" % (FLOAT, FLOAT, FLOAT)
)
PY_STATE_RE = re.compile(r"PyOQP state\s+(\d+)\s+(%s)" % FLOAT, re.I)

def run(tokens, inp, use_gpu):
    env = dict(os.environ)
    env["OQP_GPU_METC"] = "1" if use_gpu else "0"
    proc = subprocess.run(tokens + [inp], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True, env=env)
    log_path = Path(inp).with_suffix(".log")
    log_text = log_path.read_text(errors="replace") if log_path.exists() else ""
    return proc.returncode, proc.stdout, proc.stderr, log_text

def parse(log_text):
    nbf = None
    m = NBF_RE.search(log_text)
    if m:
        nbf = int(m.group(1))
    values = {}
    m = TOTAL_RE.search(log_text)
    if m:
        values["scf_total_energy_hartree"] = float(m.group(1).replace("D","E").replace("d","e"))
    in_summary = False
    for line in log_text.splitlines():
        if line.strip().startswith("State") and "Oscillator" in line:
            in_summary = True
            continue
        if in_summary:
            if not line.strip():
                continue
            m = SUMMARY_RE.match(line)
            if not m:
                if line.lstrip().startswith("Transition"):
                    in_summary = False
                continue
            state = int(m.group(1))
            energy = float(m.group(2).replace("D","E").replace("d","e"))
            exc_ev = float(m.group(3).replace("D","E").replace("d","e"))
            exc_rel = float(m.group(4).replace("D","E").replace("d","e"))
            toks = line.split()
            osc = None
            try:
                osc = float(toks[-1].replace("D","E").replace("d","e"))
            except Exception:
                pass
            values[f"state_{state}_energy_hartree"] = energy
            values[f"state_{state}_excitation_ev"] = exc_ev
            values[f"state_{state}_excitation_rel_gs_ev"] = exc_rel
            if osc is not None:
                values[f"state_{state}_oscillator_strength"] = osc
    # final energy table is a robust fallback for per-state Hartree energies.
    for m in PY_STATE_RE.finditer(log_text):
        state = int(m.group(1))
        values.setdefault(f"state_{state}_energy_hartree", float(m.group(2).replace("D","E").replace("d","e")))
    return nbf, values

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-cmd", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--log-dir", required=True)
    ap.add_argument("case", nargs="+", help="label:path")
    args = ap.parse_args()
    tokens = args.run_cmd.split()
    log_dir = Path(args.log_dir); log_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for spec in args.case:
        label, inp = spec.split(":", 1)
        case_results = {}
        nbf_seen = None
        for mode, use_gpu in [("cpu", False), ("gpu", True)]:
            rc, stdout, stderr, log_text = run(tokens, inp, use_gpu)
            (log_dir / f"{label}_{mode}.stdout").write_text(stdout)
            (log_dir / f"{label}_{mode}.stderr").write_text(stderr)
            (log_dir / f"{label}_{mode}.log").write_text(log_text)
            if rc != 0:
                raise SystemExit(f"{label} {mode} failed rc={rc}; see {log_dir}")
            nbf, values = parse(log_text)
            nbf_seen = nbf_seen or nbf
            case_results[mode] = values
        keys = sorted(set(case_results["cpu"]) | set(case_results["gpu"]))
        for key in keys:
            cpu = case_results["cpu"].get(key)
            gpu = case_results["gpu"].get(key)
            if cpu is None or gpu is None:
                status = "missing"
                diff = rel = ""
            else:
                diff_f = abs(cpu - gpu)
                rel_f = diff_f / max(1.0, abs(cpu))
                diff = f"{diff_f:.12g}"; rel = f"{rel_f:.12g}"
                status = "ok" if rel_f <= 1e-7 else "mismatch"
            rows.append({"case": label, "input": inp, "nbf": nbf_seen, "quantity": key,
                         "cpu": cpu, "gpu": gpu, "abs_diff": diff, "rel_diff": rel, "status": status})
    with open(args.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["case","input","nbf","quantity","cpu","gpu","abs_diff","rel_diff","status"])
        w.writeheader(); w.writerows(rows)
    bad = [r for r in rows if r["status"] != "ok"]
    print(f"wrote {args.out_csv} with {len(rows)} rows; bad={len(bad)}")
    return 1 if bad else 0
if __name__ == "__main__":
    sys.exit(main())
