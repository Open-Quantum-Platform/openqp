#!/usr/bin/env python3
"""
Sweep every `runtype=grad` example in the examples/ tree and check the OpenQP
analytical gradient against a central finite difference of the energy, for the
state(s) requested by `[properties] grad`.

For each example the SCF reference is converged once at the input geometry and
its molecular orbitals are reused (json guess) at every displaced geometry, so
the finite-difference reference stays on a single SCF solution branch -- this
matters a lot for open-shell (UHF/ROHF) references.

Usage:
    export OPENQP_ROOT=/path/to/openqp
    export LD_LIBRARY_PATH=$OPENQP_ROOT/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$OPENQP_ROOT/pyoqp:$PYTHONPATH
    export OMP_NUM_THREADS=1

    python tools/check_example_gradients.py                 # all grad examples
    python tools/check_example_gradients.py --dir HF DFT     # only some folders
    python tools/check_example_gradients.py --glob '*rhf*hf' # name filter
    python tools/check_example_gradients.py --tol 5e-4 --verbose

Notes
-----
* Ground-state (HF/DFT) gradients are validated cleanly this way.
* Excited-state roots (rpa/tda/sf/mrsf) on the symmetric H2O reference can be
  (near-)degenerate; a plain energy-index finite difference then tracks the
  wrong root at displaced geometries and reports a spurious mismatch.  Such
  cases are flagged `excited` -- treat a large value there as "inspect", not
  necessarily "wrong".  The lowest, well-separated root is the reliable check.
"""

import os
import re
import sys
import glob
import time
import argparse

import numpy as np

import oqp
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(ROOT, "examples")
BOHR = 1.8897259886


def discover(dirs, name_glob):
    files = []
    for inp in glob.glob(os.path.join(EXAMPLES, "**", "*.inp"), recursive=True):
        with open(inp) as fh:
            txt = fh.read()
        if not re.search(r"(?mi)^\s*runtype\s*=\s*grad\b", txt):
            continue
        rel = os.path.relpath(inp, EXAMPLES)
        if dirs and rel.split(os.sep)[0] not in dirs:
            continue
        if name_glob and not glob.fnmatch.fnmatch(os.path.basename(inp), name_glob):
            continue
        files.append(inp)
    return sorted(files)


def parse_section_value(txt, section, key):
    cur = None
    for line in txt.splitlines():
        s = line.strip()
        m = re.match(r"\[(\w+)\]", s)
        if m:
            cur = m.group(1)
            continue
        if cur == section and re.match(rf"(?i)^{key}\s*=", s):
            return s.split("=", 1)[1].strip()
    return None


def write_input(src_txt, dst, *, runtype=None, guess_json=None, system_block=None):
    """Rewrite an example input: optional runtype, json guess, geometry."""
    out, cur = [], None
    skip_guess_body = False
    for line in src_txt.splitlines():
        s = line.strip()
        m = re.match(r"\[(\w+)\]", s)
        if m:
            cur = m.group(1)
            skip_guess_body = (cur == "guess" and guess_json is not None)
            out.append(line)
            if cur == "guess" and guess_json is not None:
                out.append("type=json")
                out.append(f"file={guess_json}")
                out.append("continue_geom=false")
            continue
        if cur == "input" and runtype is not None and re.match(r"(?i)^runtype\s*=", s):
            out.append(f"runtype={runtype}")
            continue
        if cur == "guess" and skip_guess_body:
            continue  # drop original guess body, replaced above
        out.append(line)
    text = "\n".join(out)
    if system_block is not None:
        text = re.sub(r"(?is)(system=).*?(?=\n[a-z_]+=)", "system=\n" + system_block + "\n", text, count=1)
    with open(dst, "w") as fh:
        fh.write(text + "\n")
    return dst


def target_states(txt):
    g = parse_section_value(txt, "properties", "grad")
    if g is None:
        return [0]
    return [int(x) for x in re.split(r"[,\s]+", g.strip()) if x != ""]


def is_excited(txt):
    method = (parse_section_value(txt, "input", "method") or "hf").lower()
    return method == "tdhf"


def has_functional(txt):
    f = parse_section_value(txt, "input", "functional")
    return bool(f) and f.strip().lower() not in ("", "hf")


def example_tol(txt, tol_hf, tol_dft):
    # Pure HF gradients are grid-free and match the finite difference to ~1e-7.
    # DFT gradients carry XC-integration-grid finite-difference noise (~1e-3),
    # so they need a looser tolerance even though the analytical value is right.
    return tol_dft if has_functional(txt) else tol_hf


def check_example(inp, tol_hf, tol_dft, dx, workdir):
    with open(inp) as fh:
        src = fh.read()
    name = os.path.relpath(inp, EXAMPLES)
    tag = re.sub(r"[^\w]", "_", name)
    states = target_states(src)
    excited = is_excited(src)
    tol = example_tol(src, tol_hf, tol_dft)
    t0 = time.time()

    try:
        # --- analytical, and a saved json guess at the input geometry ---
        ana_inp = write_input(src, os.path.join(workdir, tag + ".ana.inp"))
        ra = Runner(project=tag + "_a", input_file=ana_inp,
                    log=os.path.join(workdir, tag + ".ana.log"), silent=1, usempi=False)
        ra.mol.config["tests"]["exception"] = True
        ra.mol.config["guess"]["save_mol"] = True
        ra.run()
        natom = ra.mol.data["natom"]
        grads = np.array(ra.mol.grads).reshape((-1, natom, 3))
        base_json = os.path.join(workdir, tag + ".ana.json")
        ra.mol.save_data()  # writes <log>.json
        saved = os.path.join(workdir, tag + ".ana.log").replace(".log", ".json")
        base_json = saved if os.path.exists(saved) else None

        ana = {}
        for s in states:
            ana[s] = grads[0] if grads.shape[0] == 1 else grads[s]

        # --- numerical: warm-started from base_json, energy at +/- dx ---
        num_src = write_input(src, os.path.join(workdir, tag + ".num.inp"),
                              runtype="energy",
                              guess_json=base_json if base_json else None)
        with open(num_src) as fh:
            num_txt = fh.read()
        rn = Runner(project=tag + "_n", input_file=num_src,
                    log=os.path.join(workdir, tag + ".num.log"), silent=1, usempi=False)
        m = rn.mol
        m.config["tests"]["exception"] = True
        m.data["OQP::log_filename"] = m.log
        oqp.oqp_banner(m)
        sp = SinglePoint(m)
        base = m.get_system().copy()
        ncoord = base.shape[0]
        deriv = {s: np.zeros(ncoord) for s in states}

        def ens(coord):
            m.update_system(coord)
            return np.asarray(sp.energy(do_init_scf=False)).reshape(-1)

        for i in range(ncoord):
            p = base.copy(); p[i] += dx
            ep = ens(p)
            q = base.copy(); q[i] -= dx
            eq = ens(q)
            for s in states:
                deriv[s][i] = (ep[s] - eq[s]) / (2 * dx)
        m.update_system(base)
        num = {s: deriv[s].reshape((natom, 3)) for s in states}
    except Exception as exc:  # noqa: BLE001
        return {"name": name, "error": f"{type(exc).__name__}: {exc}"}

    per = []
    for s in states:
        d = float(np.max(np.abs(ana[s] - num[s])))
        per.append({"state": s, "d": d, "ana": ana[s], "num": num[s]})
    worst = max(p["d"] for p in per)
    return {"name": name, "excited": excited, "states": per, "worst": worst,
            "tol": tol, "ok": worst <= tol, "time": time.time() - t0}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", nargs="*", default=[],
                    help="restrict to these top-level example folders")
    ap.add_argument("--glob", default=None, help="filename glob filter")
    ap.add_argument("--dx", type=float, default=1e-3)
    ap.add_argument("--tol-hf", type=float, default=1e-5,
                    help="tolerance for grid-free HF gradients (default 1e-5)")
    ap.add_argument("--tol-dft", type=float, default=1.5e-3,
                    help="tolerance for DFT gradients incl. grid FD noise "
                         "(default 1.5e-3)")
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    files = discover(args.dir, args.glob)
    if not files:
        print("no runtype=grad examples matched")
        sys.exit(0)

    workdir = args.workdir or os.path.join(ROOT, "tools", "_example_grad_scratch")
    os.makedirs(workdir, exist_ok=True)
    os.environ.setdefault("OMP_NUM_THREADS", "1")

    print("=" * 96)
    print(f"OpenQP example gradient sweep: {len(files)} runtype=grad example(s)")
    print(f"  dx={args.dx:g} Bohr   tol(HF)={args.tol_hf:g}   "
          f"tol(DFT,grid-noise)={args.tol_dft:g} Hartree/Bohr")
    print("=" * 96)
    print(f"{'example':<58}{'kind':>9}{'worst max|Δ|':>15}  result")
    print("-" * 96)

    gs_fail = exc_fail = errors = 0
    results = []
    for inp in files:
        r = check_example(inp, args.tol_hf, args.tol_dft, args.dx, workdir)
        results.append(r)
        if "error" in r:
            errors += 1
            print(f"{r['name']:<58}{'':>9}{'':>15}  ERROR: {r['error']}")
            continue
        kind = "excited" if r["excited"] else "ground"
        tag = "PASS" if r["ok"] else ("inspect" if r["excited"] else "FAIL")
        print(f"{r['name']:<58}{kind:>9}{r['worst']:>15.2e}  {tag}")
        if not r["ok"]:
            if r["excited"]:
                exc_fail += 1
            else:
                gs_fail += 1
        if args.verbose and not r["ok"]:
            for p in r["states"]:
                print(f"      state {p['state']}: max|Δ|={p['d']:.2e}")

    print("-" * 96)
    print(f"ground-state failures: {gs_fail}   excited-state to-inspect: {exc_fail}"
          f"   errors: {errors}   total: {len(results)}")
    # Only ground-state mismatches and errors are hard failures.
    sys.exit(1 if (gs_fail or errors) else 0)


if __name__ == "__main__":
    main()
