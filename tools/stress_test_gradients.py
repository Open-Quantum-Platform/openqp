#!/usr/bin/env python3
"""
Stress test for OpenQP analytical gradients on a spread of systems that are
prone to gradient errors: larger molecules, second-row atoms (d functions in
the basis), polarized / correlation-consistent bases, and open-shell radicals
and triplets (UHF and ROHF).

For each system it compares the OpenQP analytical gradient against
  (a) a central finite difference of OpenQP's own energy, and
  (b) PySCF's analytical gradient (if PySCF is installed) -- an independent
      reference.  OpenQP uses cartesian d/f functions for Pople basis sets, so
      PySCF is run with cart=True; the SCF energies then agree to ~1e-9 and the
      gradient comparison is unambiguous.

A case PASSES when the analytical gradient agrees with at least one trustworthy
reference within tolerance.  For open-shell systems the SCF can converge to
different solutions in the two codes (or be bistable under finite-difference
displacements); the script reports the SCF energy match so such cases are easy
to spot.

Usage:
    export OPENQP_ROOT=/path/to/openqp
    export LD_LIBRARY_PATH=$OPENQP_ROOT/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$OPENQP_ROOT/pyoqp:$PYTHONPATH
    export OMP_NUM_THREADS=1                # 1 = most reproducible
    python tools/stress_test_gradients.py
    python tools/stress_test_gradients.py --dx 5e-4 --tol 2e-4 --verbose
    python tools/stress_test_gradients.py --only co_rhf oh_uhf

Exit code is non-zero if any case fails against every available reference.
"""

import os
import sys
import time
import argparse

import numpy as np

import oqp
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

try:
    from pyscf import gto, scf, dft as pyscf_dft
    HAVE_PYSCF = True
except ImportError:
    HAVE_PYSCF = False


# Each case: geometry (element symbol + xyz, Angstrom), charge, multiplicity
# (2S+1), scf type, basis, optional DFT functional.  Geometries are sensible
# but deliberately not at the exact minimum so gradients are non-trivial.
CASES = {
    # ---- closed-shell, first/second row, various bases ----
    "co_rhf":   dict(scf="rhf",  mult=1, charge=0, basis="cc-pvdz", func=None,
                     atoms="C 0.0 0.0 0.0; O 0.0 0.0 1.200"),
    "hf_rhf":   dict(scf="rhf",  mult=1, charge=0, basis="cc-pvdz", func=None,
                     atoms="H 0.0 0.0 0.0; F 0.0 0.0 0.950"),
    "h2co_rhf": dict(scf="rhf",  mult=1, charge=0, basis="6-31g*",  func=None,
                     atoms="C 0.0 0.0 -0.520; O 0.0 0.0 0.680; "
                           "H 0.0 0.943 -1.080; H 0.0 -0.943 -1.080"),
    "nh3_rhf":  dict(scf="rhf",  mult=1, charge=0, basis="6-31g*",  func=None,
                     atoms="N 0.0 0.0 0.130; H 0.94 0.0 -0.30; "
                           "H -0.47 0.814 -0.30; H -0.47 -0.814 -0.30"),
    "h2s_rhf":  dict(scf="rhf",  mult=1, charge=0, basis="6-31g*",  func=None,
                     atoms="S 0.0 0.0 0.100; H 0.96 0.0 -0.86; H -0.96 0.0 -0.86"),
    "c2h2_rhf": dict(scf="rhf",  mult=1, charge=0, basis="6-31g*",  func=None,
                     atoms="H 0.0 0.0 -1.660; C 0.0 0.0 -0.600; "
                           "C 0.0 0.0 0.600; H 0.0 0.0 1.660"),
    # ---- open-shell radicals (doublet): UHF and ROHF ----
    "oh_uhf":   dict(scf="uhf",  mult=2, charge=0, basis="6-31g*",  func=None,
                     atoms="O 0.0 0.0 0.0; H 0.0 0.0 0.980"),
    "oh_rohf":  dict(scf="rohf", mult=2, charge=0, basis="6-31g*",  func=None,
                     atoms="O 0.0 0.0 0.0; H 0.0 0.0 0.980"),
    "nh2_uhf":  dict(scf="uhf",  mult=2, charge=0, basis="6-31g*",  func=None,
                     atoms="N 0.0 0.0 0.140; H 0.80 0.0 -0.50; H -0.80 0.0 -0.50"),
    "ch3_uhf":  dict(scf="uhf",  mult=2, charge=0, basis="6-31g*",  func=None,
                     atoms="C 0.0 0.0 0.0; H 1.05 0.0 0.10; "
                           "H -0.52 0.91 0.10; H -0.52 -0.91 0.10"),
    # ---- open-shell triplets ----
    "o2_uhf":   dict(scf="uhf",  mult=3, charge=0, basis="6-31g*",  func=None,
                     atoms="O 0.0 0.0 0.0; O 0.0 0.0 1.210"),
    "ch2_rohf": dict(scf="rohf", mult=3, charge=0, basis="6-31g*",  func=None,
                     atoms="C 0.0 0.0 0.150; H 0.87 0.0 -0.58; H -0.87 0.0 -0.58"),
    # ---- a DFT case, to confirm the (already-correct) DFT path is unchanged ----
    "h2co_bhhlyp": dict(scf="rhf", mult=1, charge=0, basis="6-31g*", func="bhhlyp",
                        atoms="C 0.0 0.0 -0.520; O 0.0 0.0 0.680; "
                              "H 0.0 0.943 -1.080; H 0.0 -0.943 -1.080"),
}


def parse_atoms(s):
    syms, xyz = [], []
    for tok in s.split(";"):
        p = tok.split()
        syms.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, np.array(xyz)


def oqp_input(case, runtype, workdir, tag):
    syms, xyz = parse_atoms(case["atoms"])
    lines = [f"{s}  {x:.10f}  {y:.10f}  {z:.10f}" for s, (x, y, z) in zip(syms, xyz)]
    func = f"functional={case['func']}\n" if case["func"] else ""
    method = "hf"
    text = f"""[input]
system=
{chr(10).join('   ' + l for l in lines)}
charge={case['charge']}
runtype={runtype}
basis={case['basis']}
{func}method={method}

[guess]
type=huckel

[scf]
type={case['scf']}
multiplicity={case['mult']}
conv=1e-10
maxit=200

[properties]
grad=0
"""
    path = os.path.join(workdir, f"{tag}.inp")
    with open(path, "w") as fh:
        fh.write(text)
    return path


def oqp_analytical(case, workdir, tag):
    inp = oqp_input(case, "grad", workdir, tag + ".ana")
    r = Runner(project=tag, input_file=inp,
               log=os.path.join(workdir, tag + ".ana.log"), silent=1, usempi=False)
    r.mol.config["tests"]["exception"] = True
    r.run()
    natom = r.mol.data["natom"]
    g = np.array(r.mol.get_grad()).reshape(natom, 3)
    e = float(np.asarray(r.mol.energies).reshape(-1)[0])
    return g, e, natom


def oqp_numerical(case, workdir, tag, dx):
    inp = oqp_input(case, "energy", workdir, tag + ".num")
    r = Runner(project=tag, input_file=inp,
               log=os.path.join(workdir, tag + ".num.log"), silent=1, usempi=False)
    m = r.mol
    m.config["tests"]["exception"] = True
    m.data["OQP::log_filename"] = m.log
    oqp.oqp_banner(m)
    sp = SinglePoint(m)
    base = m.get_system().copy()
    ncoord = base.shape[0]
    natom = m.data["natom"]
    g = np.zeros(ncoord)

    def energy(c):
        m.update_system(c)
        return float(np.asarray(sp.energy()).reshape(-1)[0])

    for i in range(ncoord):
        p = base.copy(); p[i] += dx
        mn = base.copy(); mn[i] -= dx
        g[i] = (energy(p) - energy(mn)) / (2 * dx)
    m.update_system(base)
    return g.reshape(natom, 3)


def pyscf_reference(case):
    if not HAVE_PYSCF:
        return None, None
    syms, xyz = parse_atoms(case["atoms"])
    atom = [[s, tuple(c)] for s, c in zip(syms, xyz)]
    spin = case["mult"] - 1  # 2S = (2S+1) - 1
    mol = gto.M(atom=atom, basis=case["basis"], unit="Angstrom",
                cart=True, spin=spin, charge=case["charge"], verbose=0)
    if case["func"]:
        xc = {"bhhlyp": "bhandhlyp"}.get(case["func"], case["func"])
        mf = (pyscf_dft.RKS if case["scf"] == "rhf" else
              pyscf_dft.ROKS if case["scf"] == "rohf" else pyscf_dft.UKS)(mol)
        mf.xc = xc
        mf.grids.level = 5
    else:
        mf = {"rhf": scf.RHF, "uhf": scf.UHF, "rohf": scf.ROHF}[case["scf"]](mol)
    mf.conv_tol = 1e-11
    e = mf.kernel()
    g = mf.nuc_grad_method().kernel()
    return g, e


def run_case(name, dx, tol, workdir):
    case = CASES[name]
    t0 = time.time()
    try:
        ana, e_oqp, natom = oqp_analytical(case, workdir, name)
        num = oqp_numerical(case, workdir, name, dx)
    except Exception as exc:  # noqa: BLE001
        return {"name": name, "error": f"{type(exc).__name__}: {exc}"}

    d_num = float(np.max(np.abs(ana - num)))

    py = pyscf_reference(case) if HAVE_PYSCF else (None, None)
    g_py, e_py = py
    d_py = e_diff = None
    if g_py is not None:
        d_py = float(np.max(np.abs(ana - g_py)))
        e_diff = abs(e_oqp - e_py)

    # Verdict: agree with numerical, or (energies match) agree with PySCF.
    ok_num = d_num <= tol
    ok_py = (d_py is not None and d_py <= tol and (e_diff is None or e_diff < 1e-5))
    ok = ok_num or ok_py
    return {"name": name, "natom": natom, "scf": case["scf"], "mult": case["mult"],
            "basis": case["basis"], "func": case["func"], "e_oqp": e_oqp,
            "e_py": e_py, "e_diff": e_diff, "d_num": d_num, "d_py": d_py,
            "ok": ok, "ok_num": ok_num, "ok_py": ok_py, "ana": ana, "num": num,
            "g_py": g_py, "time": time.time() - t0}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", default=[], help="subset of case names")
    ap.add_argument("--dx", type=float, default=1e-3, help="FD step, Bohr")
    ap.add_argument("--tol", type=float, default=2e-4,
                    help="max abs tolerance, Hartree/Bohr")
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    names = args.only or list(CASES)
    bad = [n for n in names if n not in CASES]
    if bad:
        ap.error(f"unknown case(s): {', '.join(bad)}; choose from {', '.join(CASES)}")

    workdir = args.workdir or os.path.join(ROOT, "tools", "_stress_scratch")
    os.makedirs(workdir, exist_ok=True)
    os.environ.setdefault("OMP_NUM_THREADS", "1")

    print("=" * 92)
    print("OpenQP gradient stress test"
          + ("" if HAVE_PYSCF else "   (PySCF not found -> numerical reference only)"))
    print(f"  dx={args.dx:g} Bohr   tol={args.tol:g} Hartree/Bohr")
    print("=" * 92)
    hdr = (f"{'case':<14}{'scf':>5}{'mult':>5}{'basis':>9}{'func':>8}"
           f"{'max|ana-num|':>14}{'max|ana-pyscf|':>16}{'dE(OQP-py)':>12}  result")
    print(hdr)
    print("-" * len(hdr))

    results = []
    for n in names:
        r = run_case(n, args.dx, args.tol, workdir)
        results.append(r)
        if "error" in r:
            print(f"{n:<14}{'':>5}{'':>5}{'':>9}{'':>8}{'':>14}{'':>16}{'':>12}  ERROR: {r['error']}")
            continue
        dpy = "n/a" if r["d_py"] is None else f"{r['d_py']:.2e}"
        ede = "n/a" if r["e_diff"] is None else f"{r['e_diff']:.1e}"
        print(f"{r['name']:<14}{r['scf']:>5}{r['mult']:>5}{r['basis']:>9}"
              f"{str(r['func']):>8}{r['d_num']:>14.2e}{dpy:>16}{ede:>12}  "
              f"{'PASS' if r['ok'] else 'FAIL'}")
        if args.verbose or not r["ok"]:
            for i in range(r["ana"].shape[0]):
                row = (f"    atom {i+1}: ana={r['ana'][i].round(6)}  "
                       f"num={r['num'][i].round(6)}")
                if r["g_py"] is not None:
                    row += f"  pyscf={r['g_py'][i].round(6)}"
                print(row)

    print("-" * len(hdr))
    failed = sum(1 for r in results if r.get("error") or not r.get("ok"))
    notes = []
    for r in results:
        if r.get("error") or r["ok"]:
            continue
        if r["e_diff"] is not None and r["e_diff"] >= 1e-5:
            notes.append(f"  {r['name']}: OQP and PySCF SCF energies differ by "
                         f"{r['e_diff']:.1e} (different SCF solution) -- compare "
                         f"OQP analytical vs OQP numerical instead: {r['d_num']:.2e}")
    if notes:
        print("Notes:")
        print("\n".join(notes))
    if failed:
        print(f"\n{failed}/{len(results)} case(s) need attention.")
        sys.exit(1)
    print(f"\nAll {len(results)} cases agree with a trusted reference within {args.tol:g}.")
    sys.exit(0)


if __name__ == "__main__":
    main()
