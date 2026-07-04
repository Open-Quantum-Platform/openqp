#!/usr/bin/env python3
"""Self-contained, resumable, parallel SCF-converger DATABASE driver.

For each (system, reference, functional, converger) runs OpenQP and records the
convergence cost (total Fock builds across escalation stages, wall time, converged?,
final energy) plus simple per-system features. Output: results/database/db.csv.
The labeled target for an ML method-selector is, per (system,ref,func), which
converger reached convergence in the fewest Fock builds.

Resumable (skips rows already in db.csv). Parallel (POOL concurrent openqp jobs).
Telegram ping on completion. Designed to run nohup-detached on zeus:

  PATH=$HOME/venv-oqp-bench/bin:/usr/local/bin:$PATH OPENQP_ROOT=$HOME/clone/openqp-bench \
    nohup python3 harness/build_database.py --root . --pool 4 --omp 2 >db.log 2>&1 &
"""
import argparse, csv, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from systems import SYSTEMS, PROFILES  # noqa: E402

BASIS = "def2-svp"
# Functionals and convergers are per-profile (systems.PROFILES) and driven per system in
# main(), so editing a system's profile (or adding a new one) actually changes its DB rows.
# Previously two hard-coded globals here overrode the profiles, silently ignoring hard-case
# coverage such as pbe / bhhlyp / vdiis.
# Optional progress pings: set TG_TOKEN and TG_CHAT in the environment (never hardcode).
TG_TOKEN = os.environ.get("TG_TOKEN", "")
TG_CHAT = os.environ.get("TG_CHAT", "")


def telegram(msg):
    if not (TG_TOKEN and TG_CHAT):
        return
    try:
        subprocess.run(["curl", "-s", f"https://api.telegram.org/bot{TG_TOKEN}/sendMessage",
                        "-d", f"chat_id={TG_CHAT}", "--data-urlencode", f"text={msg}"],
                       capture_output=True, timeout=20)
    except Exception:
        pass


def read_xyz(root, geom):
    lines = open(os.path.join(root, geom)).read().splitlines()
    n = int(lines[0].split()[0])
    return [ln for ln in lines[2:2 + n] if ln.strip()]


def conv_block(c):
    if c in ("cdiis", "ediis", "adiis", "vdiis"):
        # DIIS family: native top-level converger is 'diis' (default), diis_type picks the
        # subtype. vDIIS is a subtype too (native auto-applies its level shift).
        return ["maxdiis=10", f"diis_type={c}", "maxit=200"]
    if c == "soscf":
        return ["converger_type=soscf", "maxit=200"]
    if c == "trah":
        return ["converger_type=trah", "maxit=100"]
    raise ValueError(c)


def make_inp(root, s, ref, mult, func, c, wd):
    atoms = read_xyz(root, s["geom"])
    L = ["[input]", "system="] + [" " + a for a in atoms] + ["",
         f"charge={s['charge']}", "method=hf", f"basis={BASIS}", "runtype=energy"]
    if func:
        L += [f"functional={func}", "d4=false"]
    L += ["", "[guess]", "type=huckel", "save_mol=true", "",
          "[scf]", f"type={ref}", f"multiplicity={mult}", "conv=1.0e-8"] + conv_block(c) + [
          "", "[dftgrid]", "rad_npts=96", "ang_npts=302", "pruned="]
    p = os.path.join(wd, f"{s['id']}_{ref}_{func or 'hf'}_{c}.inp")
    open(p, "w").write("\n".join(L) + "\n")
    return p


def trah_micro_sum(txt):
    """Sum the 'micro' column of every TRAH macro table in the log.

    The native TRAH macro table is
        Macro  Energy  |grad|  rho  trust  micro  step
    and each micro-iteration is one Hessian-vector (response) Fock build. niter
    only counts macro iterations, so it UNDERCOUNTS TRAH's true cost by exactly
    this sum -- which is why labeling by niter wrongly favors TRAH.
    """
    total, in_tbl = 0, False
    for ln in txt.splitlines():
        if "Macro" in ln and "micro" in ln:
            in_tbl = True
            continue
        if not in_tbl:
            continue
        s = ln.strip()
        # The table ends only at its trailing dashed rule / SCF-convergence banner.
        # Crucially do NOT end on interleaved 'DFT XC integration time' or '===='
        # lines (DFT runs print XC timing between macro rows) -- ending there would
        # silently zero the micro count for every DFT TRAH cell.
        if s.startswith("---") or "SCF convergence" in ln or ("Final" in ln and "energy" in ln):
            in_tbl = False
            continue
        tok = ln.split()
        if len(tok) >= 7:                       # data row: idx E |grad| rho trust micro step
            try:
                int(tok[0]); float(tok[1]); total += int(tok[5])
            except (ValueError, IndexError):
                pass
        # '===' separators, 'DFT XC integration time', 'CONVERGED' rows, blanks:
        # ignore but stay in the table.
    return total


def run_cell(root, s, ref, mult, func, c, wd, omp):
    inp = make_inp(root, s, ref, mult, func, c, wd)
    log = inp[:-4] + ".log"
    env = dict(os.environ, OMP_NUM_THREADS=str(omp))
    t = time.time()
    try:
        subprocess.run(["openqp", inp], cwd=wd, env=env, capture_output=True, timeout=2400)
    except subprocess.TimeoutExpired:
        return dict(converged=False, niter=-1, micro=-1, fock_builds=-1,
                    time=round(time.time() - t, 2), E="nan")
    dt = time.time() - t
    txt = open(log).read() if os.path.exists(log) else ""
    es = re.findall(r"energy is\s+(-?\d+\.\d+)\s+after\s+(\d+)\s+iterations", txt)
    conv = "convergence achieved" in txt
    niter = sum(int(n) for _, n in es) if es else -1
    E = es[-1][0] if es else "nan"
    # True Fock-equivalent cost: each SCF iteration is ~1 Fock build; TRAH adds one
    # response build per micro-iteration. Sum micro from ANY TRAH macro table in the log,
    # not only when TRAH is the primary converger: the robustness ladder lets a DIIS/SOSCF
    # primary fall through to a TRAH escalation stage whose micro builds are real cost.
    # trah_micro_sum() returns 0 when the log has no TRAH table (pure DIIS/SOSCF runs).
    micro = trah_micro_sum(txt)
    fock_builds = (niter + micro) if niter > 0 else -1
    return dict(converged=conv, niter=niter, micro=micro, fock_builds=fock_builds,
                time=round(dt, 2), E=E)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="."); ap.add_argument("--pool", type=int, default=4)
    ap.add_argument("--omp", type=int, default=2); ap.add_argument("--maxtier", type=int, default=3)
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    wd = os.path.join(root, "results", "database"); os.makedirs(wd, exist_ok=True)
    csv_path = os.path.join(wd, "db.csv")
    fields = ["system", "tier", "ref", "mult", "func", "converger", "natom", "charge",
              "converged", "niter", "micro", "fock_builds", "time", "E"]
    done = set()
    if os.path.exists(csv_path):
        for r in csv.DictReader(open(csv_path)):
            done.add((r["system"], r["ref"], r["func"], r["converger"]))
    else:
        with open(csv_path, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=fields).writeheader()

    jobs = []
    for s in SYSTEMS:
        if s["tier"] > a.maxtier:
            continue
        prof = PROFILES[s["profile"]]
        natom = len(read_xyz(root, s["geom"]))
        for ref, mult in s["refs"].items():
            for func in prof["functionals"]:
                for c in prof["convergers"]:
                    if (s["id"], ref, func or "hf", c) in done:
                        continue
                    jobs.append((s, ref, mult, func, c, natom))
    telegram(f"DB build start: {len(jobs)} cells (pool={a.pool}) on $(hostname)")
    t0 = time.time(); n_done = 0
    with ThreadPoolExecutor(max_workers=a.pool) as ex:
        futs = {ex.submit(run_cell, root, s, ref, mult, func, c, wd, a.omp):
                (s, ref, mult, func, c, natom) for (s, ref, mult, func, c, natom) in jobs}
        for fut in as_completed(futs):
            s, ref, mult, func, c, natom = futs[fut]
            r = fut.result()
            row = dict(system=s["id"], tier=s["tier"], ref=ref, mult=mult, func=func or "hf",
                       converger=c, natom=natom, charge=s["charge"], **r)
            with open(csv_path, "a", newline="") as f:
                csv.DictWriter(f, fieldnames=fields).writerow(row)
            n_done += 1
            if n_done % 20 == 0:
                telegram(f"DB build: {n_done}/{len(jobs)} cells ({int(time.time()-t0)}s)")
    telegram(f"DB build DONE: {len(jobs)} cells in {int(time.time()-t0)}s -> {csv_path}")
    print("done", csv_path)


if __name__ == "__main__":
    main()
