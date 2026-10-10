#!/usr/bin/env python3
"""Recompute micro/fock_builds for TRAH rows of an existing db.csv from the
per-cell logs, using the corrected trah_micro_sum() parser. Non-destructive:
writes <db>.fixed.csv. Use when the harness parser was updated after a run.

  python reparse_fock.py --db results/database/db.csv --logs results/database
"""
import argparse, csv, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from build_database import trah_micro_sum  # corrected parser


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--logs", required=True)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.db)))
    fields = rows[0].keys()
    changed = 0
    for r in rows:
        if r["converger"] != "trah" or r["converged"] != "True":
            continue
        try:
            niter = int(float(r["niter"]))
        except (TypeError, ValueError):
            continue
        if niter <= 0:
            continue
        func = r["func"] or "hf"
        log = os.path.join(a.logs, f"{r['system']}_{r['ref']}_{func}_trah.log")
        if not os.path.exists(log):
            print("MISSING", log); continue
        micro = trah_micro_sum(open(log).read())
        fock = niter + micro
        if str(r.get("micro")) != str(micro) or str(r.get("fock_builds")) != str(fock):
            changed += 1
        r["micro"] = micro
        r["fock_builds"] = fock
    out = a.db.replace(".csv", ".fixed.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(fields))
        w.writeheader()
        w.writerows(rows)
    print(f"reparsed {changed} TRAH rows -> {out}")


if __name__ == "__main__":
    main()
