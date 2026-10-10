#!/usr/bin/env python3
"""ML method-selector + descriptive analysis over the converger database (db.csv).

Per (system, ref, func) the *label* is the converger with the fewest Fock builds
(among those that converged). We (1) report descriptive win-statistics, and (2)
train a classifier on cheap per-system features to predict the best converger,
evaluated by leave-one-system-out CV. The headline metric is *regret*: extra Fock
builds vs the oracle, compared with the always-C-DIIS baseline.

  python3 harness/build_selector.py --db results/database/db.csv
"""
import argparse, csv, math
from collections import defaultdict


def load(db):
    rows = list(csv.DictReader(open(db)))
    for r in rows:
        r["niter"] = int(r["niter"]); r["converged"] = (r["converged"] == "True")
        r["tier"] = int(r["tier"]); r["natom"] = int(r["natom"]); r["mult"] = int(r["mult"])
        r["charge"] = int(r["charge"])
        # True Fock-equivalent cost: for TRAH this folds in the expensive response/micro
        # builds that niter (macro count) omits, so ranking by niter falsely favors TRAH.
        # Older DBs without the column fall back to niter.
        fb = r.get("fock_builds", "")
        r["fock_builds"] = int(float(fb)) if str(fb) not in ("", "-1", "None") else r["niter"]
    return rows


def features(cell_rows):
    r = cell_rows[0]
    return {
        "natom": r["natom"], "charge": r["charge"], "mult": r["mult"], "tier": r["tier"],
        "is_dft": 1 if r["func"] != "hf" else 0,
        "open_shell": 1 if (r["ref"] != "rhf" or r["mult"] > 1) else 0,
        "log_natom": math.log(max(r["natom"], 1)),
    }


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--db", default="results/database/db.csv")
    a = ap.parse_args()
    rows = load(a.db)
    cells = defaultdict(list)
    for r in rows:
        cells[(r["system"], r["ref"], r["func"])].append(r)

    # best (oracle) converger per cell, plus cdiis baseline cost
    X, y, oracle_cost, cdiis_cost, keys = [], [], [], [], []
    wins = defaultdict(int)
    for key, crows in cells.items():
        ok = [r for r in crows if r["converged"] and r["niter"] > 0]
        if not ok:
            continue
        best = min(ok, key=lambda r: r["fock_builds"])
        cd = next((r for r in crows if r["converger"] == "cdiis" and r["converged"]), None)
        wins[best["converger"]] += 1
        X.append(features(crows)); y.append(best["converger"]); keys.append(key)
        oracle_cost.append(best["fock_builds"])
        cdiis_cost.append(cd["fock_builds"] if cd else best["fock_builds"] * 3)

    n = len(y)
    print(f"cells with >=1 converged converger: {n}")
    print("oracle best-converger wins:", dict(wins))
    print(f"total Fock builds  oracle={sum(oracle_cost)}  always-cdiis={sum(cdiis_cost)}"
          f"  (cdiis overhead {100*(sum(cdiis_cost)-sum(oracle_cost))/max(sum(oracle_cost),1):.0f}%)")

    try:
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.feature_extraction import DictVectorizer
        import numpy as np
    except Exception as e:
        print("sklearn unavailable (%s); descriptive only." % e); return

    vec = DictVectorizer(sparse=False)
    Xm = vec.fit_transform(X); ya = np.array(y)
    # leave-one-out CV (small dataset)
    pred = []
    for i in range(n):
        tr = [j for j in range(n) if j != i]
        clf = RandomForestClassifier(n_estimators=200, random_state=0)
        clf.fit(Xm[tr], ya[tr]); pred.append(clf.predict(Xm[i:i+1])[0])
    acc = sum(p == t for p, t in zip(pred, y)) / n
    # selector regret: cost of predicted converger per cell
    sel_cost = []
    for i, key in enumerate(keys):
        crows = cells[key]
        pr = next((r for r in crows if r["converger"] == pred[i] and r["converged"] and r["niter"] > 0), None)
        sel_cost.append(pr["fock_builds"] if pr else cdiis_cost[i])
    print(f"\nLOO selector accuracy: {acc:.2f}")
    print(f"total Fock builds  selector={sum(sel_cost)}  oracle={sum(oracle_cost)}"
          f"  always-cdiis={sum(cdiis_cost)}")
    print(f"selector vs cdiis: {100*(sum(cdiis_cost)-sum(sel_cost))/max(sum(cdiis_cost),1):.0f}% fewer Fock builds")


if __name__ == "__main__":
    main()
