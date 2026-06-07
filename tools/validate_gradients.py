#!/usr/bin/env python3
"""
Analytical vs. numerical gradient validation for OpenQP.

Runs the analytical gradient implemented in the OpenQP Fortran core and
compares it, component by component, against a numerical gradient obtained by
central finite differences of the (state) energy, for every electronic
structure method that provides gradients:

    * RHF        (closed-shell Hartree-Fock)          -- ground state
    * UHF        (unrestricted Hartree-Fock)          -- ground state
    * ROHF       (restricted open-shell Hartree-Fock) -- ground state
    * TDDFT      (linear-response TDDFT)              -- excited states 1..N
    * SF-TDDFT   (spin-flip TDDFT)                    -- states 1..N
    * MRSF-TDDFT (mixed-reference spin-flip TDDFT)    -- states 1..N

For the excited-state methods every requested state (default: the first 10)
is validated, since the gradient of each root has its own analytical
implementation (Z-vector + response density).  All state energies come from a
single energy evaluation per displaced geometry, so the cost is 6N energy
evaluations regardless of how many states are checked.

The numerical reference is the central finite difference

        dE_s/dx_i  ~=  ( E_s(x_i + dx) - E_s(x_i - dx) ) / (2 dx)

Coordinates handled through ``mol.get_system`` / ``mol.update_system`` are in
Bohr and the analytical gradient is returned in Hartree/Bohr, so the two are
directly comparable.

Usage
-----
        python tools/validate_gradients.py                 # all methods
        python tools/validate_gradients.py rhf tddft       # a subset
        python tools/validate_gradients.py --dx 5e-4       # custom FD step
        python tools/validate_gradients.py --nstate 10     # excited states
        python tools/validate_gradients.py --tol 1e-4 --verbose

Exit code is non-zero if any checked state of any method exceeds the
tolerance, so the script can be wired into CI.
"""

import os
import re
import sys
import time
import argparse

import numpy as np

import oqp
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint

# Repository root (this file lives in <root>/tools/)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(ROOT, "examples")


# key -> (label, example input file, kind)
#   kind == "hf"  : ground-state gradient only (state 0)
#   kind == "td"  : excited-state gradients (states 1..nstate)
METHODS = {
    "rhf": ("RHF", os.path.join(EXAMPLES, "HF", "H2O_RHF-HF_GRADIENT.inp"), "hf"),
    "uhf": ("UHF", os.path.join(EXAMPLES, "HF", "H2O_UHF-HF_GRADIENT.inp"), "hf"),
    "rohf": ("ROHF", os.path.join(EXAMPLES, "HF", "H2O_ROHF-HF_GRADIENT.inp"), "hf"),
    "tddft": ("TDDFT", os.path.join(EXAMPLES, "TDDFT", "H2O_B3LYP5-TDDFT_GRADIENT.inp"), "td"),
    "sf-tddft": ("SF-TDDFT", os.path.join(EXAMPLES, "SF-TDDFT", "H2O_BHHLYP-SFTDDFT_GRADIENT.inp"), "td"),
    "mrsf-tddft": ("MRSF-TDDFT", os.path.join(EXAMPLES, "MRSF-TDDFT", "H2O_BHHLYP-MRSFTDDFT_GRADIENT.inp"), "td"),
}


# Tight convergence thresholds keep the finite-difference reference clean
# (energy noise epsilon shows up in the numerical gradient as epsilon / 2dx).
SCF_CONV = "1e-10"
TD_CONV = "1e-9"
ZV_CONV = "1e-9"


def _write_input(src_inp, dst_inp, runtype, grad_list, nstate=None):
    """Copy an example input file, overriding runtype, grad list, nstate, and
    injecting tight convergence thresholds.

    Returns the path to the written file.  Editing the text (rather than the
    parsed config) guarantees the values flow through OpenQP's normal input
    parser, which is where ``nstate`` is actually consumed.
    """
    with open(src_inp) as fh:
        lines = fh.readlines()

    grad_str = ",".join(str(x) for x in grad_list)
    out = []
    section = None
    has = {"grad": False, "nstate": False, "tdhf": False, "properties": False}

    for raw in lines:
        line = raw.rstrip("\n")
        stripped = line.strip()
        m = re.match(r"\[(\w+)\]", stripped)
        if m:
            section = m.group(1)
            out.append(line)
            # Inject convergence right after the relevant section header so we
            # never touch comment lines that merely mention "[scf]"/"[tdhf]".
            if section == "scf":
                out.append(f"conv={SCF_CONV}")
            elif section == "tdhf":
                has["tdhf"] = True
                out.append(f"conv={TD_CONV}")
                out.append(f"zvconv={ZV_CONV}")
                if nstate is not None:
                    out.append(f"nstate={nstate}")
                    has["nstate"] = True
            elif section == "properties":
                has["properties"] = True
            continue
        low = stripped.lower()
        # Drop pre-existing keys we are overriding (only inside their section).
        if section == "input" and low.startswith("runtype"):
            out.append(f"runtype={runtype}")
            continue
        if section == "scf" and low.startswith("conv"):
            continue
        if section == "tdhf" and (low.startswith("conv")
                                  or low.startswith("zvconv")):
            continue
        if section == "tdhf" and low.startswith("nstate"):
            continue  # already injected after the header
        if section == "properties" and low.startswith("grad"):
            out.append(f"grad={grad_str}")
            has["grad"] = True
            continue
        out.append(line)

    if not has["properties"]:
        out.append("[properties]")
        out.append(f"grad={grad_str}")
    elif not has["grad"]:
        out.append(f"grad={grad_str}")

    with open(dst_inp, "w") as fh:
        fh.write("\n".join(out) + "\n")
    return dst_inp


def _runner(inp, workdir, tag):
    project = os.path.splitext(os.path.basename(inp))[0] + f".{tag}"
    log = os.path.join(workdir, project + ".log")
    runner = Runner(project=project, input_file=inp, log=log,
                    silent=1, usempi=False)
    runner.mol.config["tests"]["exception"] = True
    return runner


def analytical_gradients(key, workdir, nstate, compute_nstate):
    """Return {state: (natom,3)} from the analytical implementation."""
    label, src, kind = METHODS[key]
    if kind == "hf":
        states = [0]
        grad_list = [0]
        ns = None
    else:
        states = list(range(1, nstate + 1))
        grad_list = states
        ns = compute_nstate

    inp = _write_input(src, os.path.join(workdir, f"{key}.ana.inp"),
                       runtype="grad", grad_list=grad_list, nstate=ns)
    runner = _runner(inp, workdir, "ana")
    runner.run()
    natom = runner.mol.data["natom"]
    grads = np.array(runner.mol.grads).reshape((-1, natom, 3))

    out = {}
    for s in states:
        out[s] = grads[0] if kind == "hf" else grads[s]
    return out, natom, states


def numerical_gradients(key, workdir, states, dx, nstate, compute_nstate):
    """Central finite-difference gradient of every checked state's energy."""
    label, src, kind = METHODS[key]
    ns = None if kind == "hf" else compute_nstate
    inp = _write_input(src, os.path.join(workdir, f"{key}.num.inp"),
                       runtype="energy", grad_list=states, nstate=ns)
    runner = _runner(inp, workdir, "num")
    mol = runner.mol
    # Runner.run() normally performs this setup; we drive SinglePoint directly.
    mol.data["OQP::log_filename"] = mol.log
    oqp.oqp_banner(mol)
    sp = SinglePoint(mol)

    base = mol.get_system().copy()  # 3N, Bohr
    ncoord = base.shape[0]
    natom = mol.data["natom"]
    deriv = {s: np.zeros(ncoord) for s in states}

    def state_energies(coord):
        mol.update_system(coord)
        return np.asarray(sp.energy()).reshape(-1)

    for i in range(ncoord):
        plus = base.copy(); plus[i] += dx
        e_plus = state_energies(plus)
        minus = base.copy(); minus[i] -= dx
        e_minus = state_energies(minus)
        for s in states:
            deriv[s][i] = (e_plus[s] - e_minus[s]) / (2.0 * dx)

    mol.update_system(base)
    return {s: deriv[s].reshape((natom, 3)) for s in states}


def validate(key, dx, nstate, workdir):
    label, src, kind = METHODS[key]
    if not os.path.exists(src):
        return {"label": label, "error": f"input not found: {src}"}

    # Compute a couple of extra excited states so the highest *checked* root is
    # not the last Davidson/Z-vector root (improves Z-vector convergence).
    compute_nstate = nstate + 2 if kind == "td" else None

    t0 = time.time()
    try:
        ana, natom, states = analytical_gradients(key, workdir, nstate, compute_nstate)
        num = numerical_gradients(key, workdir, states, dx, nstate, compute_nstate)
    except Exception as exc:  # noqa: BLE001 - report any failure cleanly
        return {"label": label, "error": f"{type(exc).__name__}: {exc}"}

    per_state = []
    for s in states:
        diff = ana[s] - num[s]
        per_state.append({
            "state": s,
            "max_abs": float(np.max(np.abs(diff))),
            "rms": float(np.sqrt(np.mean(diff ** 2))),
            "ana": ana[s],
            "num": num[s],
        })
    return {"label": label, "kind": kind, "natom": natom,
            "states": per_state, "time": time.time() - t0}


def format_grad_block(state, ana, num):
    lines = [f"    state {state}: analytical | numerical | diff (Hartree/Bohr)"]
    fa, fn = ana.reshape(-1, 3), num.reshape(-1, 3)
    for at in range(fa.shape[0]):
        for c, axis in enumerate("xyz"):
            a, n = fa[at, c], fn[at, c]
            lines.append(f"      atom {at + 1} {axis}: "
                         f"{a:>15.8f} | {n:>15.8f} | {a - n:>+12.2e}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("methods", nargs="*", default=[],
                        help="subset of methods (default: all of "
                             f"{', '.join(METHODS)})")
    parser.add_argument("--dx", type=float, default=1e-3,
                        help="finite-difference step in Bohr (default: 1e-3)")
    parser.add_argument("--tol", type=float, default=1e-4,
                        help="max abs error tolerance, Hartree/Bohr (default: 1e-4)")
    parser.add_argument("--nstate", type=int, default=10,
                        help="number of excited states to validate for the TD "
                             "methods (default: 10)")
    parser.add_argument("--workdir", default=None,
                        help="scratch directory for input/log files")
    parser.add_argument("--verbose", action="store_true",
                        help="print full per-component gradient tables")
    args = parser.parse_args()

    keys = args.methods or list(METHODS)
    unknown = [k for k in keys if k not in METHODS]
    if unknown:
        parser.error(f"unknown method(s): {', '.join(unknown)}. "
                     f"choose from {', '.join(METHODS)}")

    workdir = args.workdir or os.path.join(ROOT, "tools",
                                           "_grad_validation_scratch")
    os.makedirs(workdir, exist_ok=True)
    os.environ.setdefault("OMP_NUM_THREADS", "1")

    print("=" * 80)
    print("OpenQP analytical vs. numerical gradient validation")
    print(f"  dx = {args.dx:g} Bohr   tol = {args.tol:g} Hartree/Bohr   "
          f"excited states checked = 1..{args.nstate}")
    print("=" * 80)

    results = []
    for key in keys:
        res = validate(key, args.dx, args.nstate, workdir)
        results.append(res)
        if "error" in res:
            print(f"\n[{res['label']:<11}] ERROR: {res['error']}")
            continue
        worst = max(s["max_abs"] for s in res["states"])
        ok = worst <= args.tol
        print(f"\n[{res['label']:<11}] {len(res['states'])} state(s)  "
              f"worst max|Δ|={worst:.2e}  ({res['time']:.1f}s)  "
              f"-> {'PASS' if ok else 'FAIL'}")
        for st in res["states"]:
            st_ok = st["max_abs"] <= args.tol
            print(f"    state {st['state']:>2}: max|Δ|={st['max_abs']:.2e}  "
                  f"rms={st['rms']:.2e}  -> {'pass' if st_ok else 'FAIL'}")
            if args.verbose or not st_ok:
                print(format_grad_block(st["state"], st["ana"], st["num"]))

    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"{'method':<12}{'states':>8}{'worst max|Δ|':>16}   result")
    failed = 0
    for res in results:
        if "error" in res:
            failed += 1
            print(f"{res['label']:<12}{'-':>8}{'-':>16}   ERROR")
            continue
        worst = max(s["max_abs"] for s in res["states"])
        ok = worst <= args.tol
        failed += 0 if ok else 1
        print(f"{res['label']:<12}{len(res['states']):>8}"
              f"{worst:>16.2e}   {'PASS' if ok else 'FAIL'}")
    print("=" * 80)

    if failed:
        print(f"{failed} method(s) failed / errored.")
        sys.exit(1)
    print("All methods agree with numerical gradients within tolerance.")
    sys.exit(0)


if __name__ == "__main__":
    main()
