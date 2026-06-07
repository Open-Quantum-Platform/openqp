#!/usr/bin/env python3
"""
Extensive excited-state analytical-vs-numerical gradient validation.

Checks the analytical gradient of every excited-state method/functional
combination against a 5-point central finite difference of the (state) energy,
for states 1..N (default 5), at a low-symmetry H2O/6-31G* geometry that lifts
the C2v degeneracies so the energy-index finite difference tracks a consistent
root.

Tolerances are method-aware:
  * pure-HF references (no XC grid)  -> 1e-5 Hartree/Bohr
  * DFT references (XC-grid FD floor) -> 5e-4 Hartree/Bohr

See tools/EXCITED_GRADIENT_VALIDATION.md for the recorded results and the
discussion of near-degenerate-root reference artifacts (the upper member of a
near-degenerate pair cannot be validated by an energy-index finite
difference; that is a property of the reference, not the analytical gradient).

Usage:
    export OPENQP_ROOT=/path/to/openqp
    export LD_LIBRARY_PATH=$OPENQP_ROOT/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$OPENQP_ROOT/pyoqp:$PYTHONPATH
    python tools/excited_gradient_matrix.py                 # full matrix
    python tools/excited_gradient_matrix.py rpa tda         # families by name
    python tools/excited_gradient_matrix.py --nstate 5 --dx 1e-3
"""

import os
import sys
import argparse

import numpy as np

import oqp
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint

WD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_excited_grad_scratch")

# Low-symmetry water (Angstrom): no C2v degeneracies for the low roots.
GEOM = """   8   0.07000000   0.02000000   0.01000000
   1   0.95000000   0.06000000  -0.03000000
   1  -0.25000000   0.93000000   0.10000000"""

# label, scf type, multiplicity, functional, tdhf type
MATRIX = [
    ("RPA/b3lyp5", "rhf", 1, "b3lyp5", "rpa"),
    ("RPA/bhhlyp", "rhf", 1, "bhhlyp", "rpa"),
    ("RPA/hf", "rhf", 1, "", "rpa"),
    ("RPA/pbe", "rhf", 1, "pbe", "rpa"),
    ("TDA/b3lyp5", "rhf", 1, "b3lyp5", "tda"),
    ("TDA/bhhlyp", "rhf", 1, "bhhlyp", "tda"),
    ("TDA/hf", "rhf", 1, "", "tda"),
    ("TDA/pbe", "rhf", 1, "pbe", "tda"),
    ("SF/bhhlyp", "rohf", 3, "bhhlyp", "sf"),
    ("SF/b3lyp5", "rohf", 3, "b3lyp5", "sf"),
    ("SF/hf", "rohf", 3, "", "sf"),
    ("MRSF/bhhlyp", "rohf", 3, "bhhlyp", "mrsf"),
    ("MRSF/b3lyp5", "rohf", 3, "b3lyp5", "mrsf"),
    ("MRSF/hf", "rohf", 3, "", "mrsf"),
    ("MRSF/cam-b3lyp", "rohf", 3, "cam-b3lyp", "mrsf"),
]


def _inp(scftype, mult, func, tdtype, nstate, target, runtype):
    funcline = f"functional={func}\n" if func else ""
    return f"""[input]
system=
{GEOM}
charge=0
runtype={runtype}
basis=6-31g*
{funcline}method=tdhf
[guess]
type=huckel
[scf]
type={scftype}
multiplicity={mult}
conv=1e-11
[tdhf]
type={tdtype}
nstate={nstate}
conv=1e-10
zvconv=1e-10
[properties]
grad={target}
"""


def _analytic(scftype, mult, func, tdtype, nstate, target):
    p = os.path.join(WD, "g.inp")
    open(p, "w").write(_inp(scftype, mult, func, tdtype, nstate, target, "grad"))
    r = Runner(project="g", input_file=p, log=os.path.join(WD, "g.log"),
               silent=1, usempi=False)
    r.mol.config["tests"]["exception"] = True
    r.run()
    nat = r.mol.data["natom"]
    return np.array(r.mol.grads).reshape((-1, nat, 3))[target]


def _numeric(scftype, mult, func, tdtype, nstate, target, dx):
    p = os.path.join(WD, "n.inp")
    open(p, "w").write(_inp(scftype, mult, func, tdtype, nstate, target, "energy"))
    r = Runner(project="n", input_file=p, log=os.path.join(WD, "n.log"),
               silent=1, usempi=False)
    m = r.mol
    m.config["tests"]["exception"] = True
    m.data["OQP::log_filename"] = m.log
    oqp.oqp_banner(m)
    sp = SinglePoint(m)
    base = m.get_system().copy()
    nat = m.data["natom"]
    g = np.zeros(base.shape[0])

    def en(c):
        m.update_system(c)
        return np.asarray(sp.energy()).reshape(-1)[target]

    for i in range(base.shape[0]):
        cp = base.copy(); cp[i] += dx; ep1 = en(cp)
        cp = base.copy(); cp[i] += 2 * dx; ep2 = en(cp)
        cp = base.copy(); cp[i] -= dx; em1 = en(cp)
        cp = base.copy(); cp[i] -= 2 * dx; em2 = en(cp)
        g[i] = (-ep2 + 8 * ep1 - 8 * em1 + em2) / (12 * dx)
    m.update_system(base)
    return g.reshape((nat, 3))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("families", nargs="*", default=[],
                    help="optional filter: rpa tda sf mrsf (matches the label prefix)")
    ap.add_argument("--nstate", type=int, default=5, help="states to check (default 5)")
    ap.add_argument("--dx", type=float, default=1e-3, help="FD step, Bohr")
    ap.add_argument("--tol-hf", type=float, default=1e-5)
    ap.add_argument("--tol-dft", type=float, default=5e-4)
    args = ap.parse_args()

    os.makedirs(WD, exist_ok=True)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    compute = args.nstate + 2

    rows = MATRIX
    if args.families:
        keys = {f.lower() for f in args.families}
        rows = [r for r in MATRIX if r[0].split("/")[0].lower() in keys]

    print("=" * 78)
    print("Excited-state analytical vs numerical gradient (low-symmetry H2O/6-31G*)")
    print(f"  5-point FD dx={args.dx:g}  states 1..{args.nstate}  "
          f"tol(HF)={args.tol_hf:g} tol(DFT)={args.tol_dft:g}")
    print("=" * 78)

    failed = 0
    for label, scftype, mult, func, tdtype in rows:
        try:
            worst = 0.0
            per = []
            for s in range(1, args.nstate + 1):
                a = _analytic(scftype, mult, func, tdtype, compute, s)
                n = _numeric(scftype, mult, func, tdtype, compute, s, args.dx)
                d = float(np.max(np.abs(a - n)))
                per.append(d)
                worst = max(worst, d)
            tol = args.tol_hf if func == "" else args.tol_dft
            ok = worst <= tol
            failed += 0 if ok else 1
            print(f"[{label:<16}] worst={worst:.2e} (tol {tol:.0e}) "
                  f"{'PASS' if ok else 'CHECK'}  states: "
                  + " ".join(f"{d:.1e}" for d in per))
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"[{label:<16}] ERROR {type(exc).__name__}: {exc}")

    print("=" * 78)
    print(f"{len(rows) - failed}/{len(rows)} combinations within tolerance.")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
