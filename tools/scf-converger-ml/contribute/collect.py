#!/usr/bin/env python3
"""Auto-accumulate SCF-convergence data from your OpenQP runs.

Point it at a directory of OpenQP output `.log` files (your normal calculations) and it
parses each into the standard contribution row, appending to one CSV you can PR back:

    system,ref,func,converger,natom,charge,converged,niter,energy,basis

It is incremental/idempotent: re-running only adds new (file) records. Two modes:

  # 1) harvest existing logs (passive: just run your jobs normally, then collect)
  python contribute/collect.py --logs ~/my_openqp_runs --out data/contrib/yourname.csv

  # 2) continuous watch (auto-accumulate as new logs appear)
  python contribute/collect.py --logs ~/my_openqp_runs --out data/contrib/yourname.csv --watch 60

Then open a Pull Request adding your CSV (and any new geometries to geometries/). Hard cases
where a converger FAILS are especially valuable. See README for the data policy.
"""
import argparse, csv, glob, os, re, time

_TM = {'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
       'La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg'}
FIELDS = ["system","ref","func","converger","natom","charge","transition_metal","converged","niter","energy","basis"]


def parse_log(path):
    t = open(path, errors="ignore").read()
    if "MODULE" not in t and "energy is" not in t:
        return None
    def g(pat, d=""):
        m = re.search(pat, t, re.I); return m.group(1).strip() if m else d
    basis = g(r"basis\s*=?\s*([\w\-\*\(\)]+)")
    func = g(r"functional\s*=?\s*([\w\-]+)") or "hf"
    ref = (g(r"\bscf\.?type\s*=?\s*(rhf|rohf|uhf)") or g(r"\b(RHF|ROHF|UHF) energy")).lower() or "rhf"
    charge = g(r"charge\s*=?\s*(-?\d+)", "0")
    conv = "convergence achieved" in t
    es = re.findall(r"energy is\s+(-?\d+\.\d+)\s+after\s+(\d+)\s+iterations", t)
    # TRAH macro 'micro' column -> add response builds to the cost
    micro = sum(int(m.group(1)) for m in re.finditer(
        r"^\s*\d+\s+-?\d+\.\d+\s+[\dEe.+-]+\s+[\dEe.+-]+\s+[\d.]+\s+(\d+)\s+\w+", t, re.M))
    niter = (sum(int(n) for _, n in es) + micro) if es else -1
    energy = es[-1][0] if es else "nan"
    conver = g(r"Converger\s*=\s*([\w\-]+)") or g(r"diis type\s*:\s*(\w+)") or "cdiis"
    # natom + TM from the echoed geometry (lines of: <El> x y z)
    atoms = re.findall(r"^\s*([A-Z][a-z]?)\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+", t, re.M)
    natom = len(atoms); tm = 1 if any(a in _TM for a in atoms) else 0
    sysname = os.path.splitext(os.path.basename(path))[0]
    return dict(system=sysname, ref=ref, func=func, converger=conver.lower(), natom=natom,
                charge=charge, transition_metal=tm, converged=conv, niter=niter,
                energy=energy, basis=basis)


def collect(logs, out):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    seen = set()
    if os.path.exists(out):
        for r in csv.DictReader(open(out)):
            seen.add((r["system"], r["ref"], r["func"], r["converger"]))
    new = 0
    write_header = not os.path.exists(out)
    with open(out, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if write_header:
            w.writeheader()
        for lg in sorted(glob.glob(os.path.join(logs, "**", "*.log"), recursive=True)):
            r = parse_log(lg)
            if not r:
                continue
            k = (r["system"], r["ref"], r["func"], r["converger"])
            if k in seen:
                continue
            w.writerow(r); seen.add(k); new += 1
    return new


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--logs", required=True, help="directory of OpenQP .log files")
    ap.add_argument("--out", default="data/contrib/contribution.csv")
    ap.add_argument("--watch", type=int, default=0, help="seconds; keep harvesting new logs")
    a = ap.parse_args()
    n = collect(a.logs, a.out); print("added %d records -> %s" % (n, a.out))
    while a.watch:
        time.sleep(a.watch)
        n = collect(a.logs, a.out)
        if n:
            print("added %d new -> %s" % (n, a.out))


if __name__ == "__main__":
    main()
