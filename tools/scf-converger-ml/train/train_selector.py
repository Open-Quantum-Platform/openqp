#!/usr/bin/env python3
"""Train the SCF method-selector and distill it to a tiny plain-Python model.

Input : an OpenQP per-converger database CSV with columns
        system,tier,ref,mult,func,converger,natom,charge,converged,niter,time,E
Label : per (system,ref,func) the converger with the fewest Fock builds (converged).
Output: model/scf_selector_model.py  -- a dependency-free predict(features)->converger
        that OpenQP's converger_type=ml loads at runtime.

Usage : python train/train_selector.py --db data/db_openqp_postfix.csv
"""
import argparse, csv, collections, os, glob

TM = set(range(21, 31)) | set(range(39, 49)) | set(range(57, 81))
# symbol->Z covering the SAME ranges the runtime _scf_features treats as transition metals
# (21-30, 39-48, 57-80). Stopping at Cd (48) made the geometry-based fallback label common
# 5d/lanthanide systems (Hf, Pt, Au, Hg, ...) as transition_metal=0 while production predicts
# with 1, corrupting any split on that feature.
_SYM2Z = {
    'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,
    'Y':39,'Zr':40,'Nb':41,'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,
    'La':57,'Ce':58,'Pr':59,'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,
    'Ho':67,'Er':68,'Tm':69,'Yb':70,'Lu':71,'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,
    'Ir':77,'Pt':78,'Au':79,'Hg':80}


def has_tm(system, geomdir):
    for fn in glob.glob(os.path.join(geomdir, system + ".xyz")):
        for ln in open(fn).read().splitlines()[2:]:
            s = ln.split()
            if s and s[0] in _SYM2Z:
                return 1
    return 0


def _int(r, k, d=0):
    try:
        return int(float(r.get(k, d)))
    except (TypeError, ValueError):
        return d


def features(r, geomdir):
    mult = _int(r, "mult", 1)
    tm = r.get("transition_metal")
    tm = int(tm) if (tm not in (None, "")) else has_tm(r["system"], geomdir)
    ref = (r.get("ref", "rhf") or "rhf").lower()
    return {
        "natom": _int(r, "natom"), "charge": _int(r, "charge"), "mult": mult,
        "is_dft": 0 if r.get("func", "hf") in ("hf", "") else 1,
        "open_shell": 1 if (ref != "rhf" or mult > 1) else 0,
        # Distinguish the SCF reference: ROHF and UHF converge differently, and
        # ROHF (ROKS) is the MRSF-TDDFT reference, so it gets its own feature.
        "is_rohf": 1 if ref == "rohf" else 0,
        "is_uhf": 1 if ref == "uhf" else 0,
        "transition_metal": tm,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", nargs="*",
                    default=["data/db_openqp_postfix.csv", "data/contrib/*.csv"],
                    help="CSV files or globs, unioned: the canonical post-fix OpenQP DB "
                         "(true Fock-cost labels) plus community contributions. Obsolete/"
                         "different-schema studies live in data/studies/ and are excluded.")
    ap.add_argument("--geom", default="geometries")
    ap.add_argument("--out", default="model/scf_selector_model.py")
    a = ap.parse_args()
    files = [fn for pat in a.db for fn in (glob.glob(pat) if any(c in pat for c in "*?[") else [pat])]
    rows = []
    for fn in files:
        if os.path.exists(fn):
            rows += list(csv.DictReader(open(fn)))
    print("data files:", len(files), "| rows:", len(rows))
    cells = collections.defaultdict(dict)
    feat = {}
    for r in rows:
        key = (r["system"], r["ref"], r["func"])
        if r["converged"] == "True" and int(float(r["niter"])) > 0:
            # Cost = true Fock-equivalent builds when present (TRAH's response/micro
            # builds counted), else fall back to niter for older contributed DBs.
            fb = r.get("fock_builds", "")
            cost = int(float(fb)) if fb not in (None, "", "-1") else int(float(r["niter"]))
            cells[key][r["converger"]] = cost
        feat.setdefault(key, features(r, a.geom))
    X, y, keys = [], [], []
    for k, cc in cells.items():
        if not cc:
            continue
        X.append(feat[k]); y.append(min(cc, key=cc.get)); keys.append(k)
    print("labeled cells:", len(y), "| label distribution:", dict(collections.Counter(y)))

    from sklearn.tree import DecisionTreeClassifier, export_text
    from sklearn.feature_extraction import DictVectorizer
    import numpy as np
    vec = DictVectorizer(sparse=False); Xm = vec.fit_transform(X); ya = np.array(y)
    # leave-one-system-out CV
    sysarr = np.array([k[0] for k in keys]); pred = [None]*len(y)
    for s in set(sysarr):
        te = sysarr == s; tr = ~te
        clf = DecisionTreeClassifier(max_depth=4, random_state=0).fit(Xm[tr], ya[tr])
        for i in np.where(te)[0]:
            pred[i] = clf.predict(Xm[i:i+1])[0]
    acc = sum(p == t for p, t in zip(pred, y)) / len(y)
    print("LOO-by-system accuracy: %.2f" % acc)

    # final tree on all data + distill to plain python
    clf = DecisionTreeClassifier(max_depth=4, random_state=0).fit(Xm, ya)
    print("\n" + export_text(clf, feature_names=list(vec.feature_names_)))
    _distill(clf, vec, a.out)
    print("distilled ->", a.out)


def _distill(clf, vec, out):
    t = clf.tree_; names = list(vec.feature_names_); cls = clf.classes_
    lines = ['"""Distilled SCF method-selector (auto-generated by train_selector.py).',
             'predict(features:dict)->str ; features: natom,charge,mult,is_dft,open_shell,is_rohf,is_uhf,transition_metal."""',
             '', 'def predict(f):']
    def rec(node, d):
        ind = "    " * (d + 1)
        if t.children_left[node] == t.children_right[node]:
            import numpy as np
            lines.append(ind + 'return "%s"' % cls[int(np.argmax(t.value[node][0]))])
            return
        feat = names[t.feature[node]]; thr = t.threshold[node]
        lines.append(ind + 'if f.get("%s",0) <= %.6g:' % (feat, thr))
        rec(t.children_left[node], d + 1)
        lines.append(ind + 'else:')
        rec(t.children_right[node], d + 1)
    rec(0, 0)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
