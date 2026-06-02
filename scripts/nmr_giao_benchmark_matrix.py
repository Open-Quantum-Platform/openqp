#!/usr/bin/env python3
"""Run the NMR CGO/GIAO benchmark matrix scaffold.

This is an integrated OpenQP benchmark driver, not a production GIAO
implementation. It records the current validation state explicitly:

* trusted-reference columns are generated with PySCF CGO/GIAO NMR;
* OpenQP CGO may be run through the pip-installed OpenQP driver;
* OpenQP GIAO is expected to remain gated until native London integral/response
  terms are connected. A gated/NotImplemented result is recorded as such, never
  converted to CGO numbers.

The decisive metric is gauge-origin translation sensitivity. For PySCF, CGO is
sampled with several fixed origins while GIAO uses ``gauge_orig=None`` and should
be invariant to translation of the molecular coordinates. For OpenQP, the same
metric is populated only after native GIAO is available.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import resource
import re
import subprocess
import sys
import tempfile
import threading
import time
from contextlib import contextmanager
from pathlib import Path

try:
    import psutil
except ImportError:  # pragma: no cover - exercised on minimal environments
    psutil = None

np = None
dft = None
gto = None
scf = None
nmr_rhf = None
nmr_rks = None
UNIT = None

ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "tests" / "fixtures" / "nmr" / "giao_benchmark_matrix.json"
OUTDIR = ROOT / "tests" / "fixtures" / "nmr" / "benchmark_results"


def _ensure_pyscf_imports():
    """Import PySCF lazily after early OpenQP runtime checks."""
    global np, dft, gto, scf, nmr_rhf, nmr_rks, UNIT
    if np is not None:
        return
    import numpy as _np
    from pyscf import dft as _dft, gto as _gto, scf as _scf
    from pyscf.data import nist
    from pyscf.prop.nmr import rhf as _nmr_rhf, rks as _nmr_rks

    np = _np
    dft = _dft
    gto = _gto
    scf = _scf
    nmr_rhf = _nmr_rhf
    nmr_rks = _nmr_rks
    UNIT = nist.ALPHA ** 2 * 1e6


def _ru_maxrss_mb(who=resource.RUSAGE_SELF) -> float:
    """Return ru_maxrss in MiB with macOS/Linux unit normalization."""
    value = float(resource.getrusage(who).ru_maxrss)
    # macOS reports bytes; Linux reports KiB.
    return value / (1024.0 * 1024.0) if sys.platform == "darwin" else value / 1024.0


class PeakMemorySampler:
    """Poll RSS while a benchmark block runs and report peak MiB.

    If psutil is unavailable, fall back to ru_maxrss. The fallback is a process
    high-water mark, so it is conservative rather than a clean per-block peak.
    """

    def __init__(self, pid: int | None = None, interval: float = 0.01):
        self.pid = pid or os.getpid()
        self.interval = interval
        self.peak_mb = 0.0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._proc = psutil.Process(self.pid) if psutil is not None else None

    def __enter__(self):
        if self._proc is None:
            self.peak_mb = _ru_maxrss_mb()
            return self
        self.peak_mb = self._rss_mb()
        self._thread = threading.Thread(target=self._poll, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        if self._proc is None:
            self.peak_mb = max(self.peak_mb, _ru_maxrss_mb())
            return False
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1.0)
        self.peak_mb = max(self.peak_mb, self._rss_mb())
        return False

    def _rss_mb(self) -> float:
        try:
            return float(self._proc.memory_info().rss) / (1024.0 * 1024.0)
        except Exception:
            return self.peak_mb

    def _poll(self):
        while not self._stop.is_set():
            self.peak_mb = max(self.peak_mb, self._rss_mb())
            self._stop.wait(self.interval)


@contextmanager
def _child_memory_delta():
    """Fallback peak-memory context for subprocesses without psutil."""
    before = _ru_maxrss_mb(resource.RUSAGE_CHILDREN)
    result = {"peak_memory_mb": 0.0}
    yield result
    after = _ru_maxrss_mb(resource.RUSAGE_CHILDREN)
    result["peak_memory_mb"] = max(0.0, after - before)


GEOMS_ANG = {
    "h2o": """O  0.000000000  0.000000000 -0.041061554
              H -0.533194329  0.533194329 -0.614469223
              H  0.533194329 -0.533194329 -0.614469223""",
    "hf": """H 0.000000000 0.000000000 0.000000000
             F 0.000000000 0.000000000 0.916800000""",
    "ch4": """C  0.000000000  0.000000000  0.000000000
              H  0.629118000  0.629118000  0.629118000
              H -0.629118000 -0.629118000  0.629118000
              H -0.629118000  0.629118000 -0.629118000
              H  0.629118000 -0.629118000 -0.629118000""",
    "formaldehyde": """C  0.000000000  0.000000000  0.000000000
                         O  0.000000000  0.000000000  1.205000000
                         H  0.942900000  0.000000000 -0.587600000
                         H -0.942900000  0.000000000 -0.587600000""",
    "benzene": """C  1.397000000  0.000000000 0.000000000
                  C  0.698500000  1.209837000 0.000000000
                  C -0.698500000  1.209837000 0.000000000
                  C -1.397000000  0.000000000 0.000000000
                  C -0.698500000 -1.209837000 0.000000000
                  C  0.698500000 -1.209837000 0.000000000
                  H  2.481000000  0.000000000 0.000000000
                  H  1.240500000  2.148584000 0.000000000
                  H -1.240500000  2.148584000 0.000000000
                  H -2.481000000  0.000000000 0.000000000
                  H -1.240500000 -2.148584000 0.000000000
                  H  1.240500000 -2.148584000 0.000000000""",
    "glycine": """N -1.453000000  0.000000000  0.000000000
                 C -0.248000000  0.000000000  0.000000000
                 C  0.520000000  1.238000000  0.000000000
                 O  1.733000000  1.231000000  0.000000000
                 O -0.170000000  2.334000000  0.000000000
                 H -1.990000000  0.833000000  0.000000000
                 H -1.990000000 -0.833000000  0.000000000
                 H -0.214000000 -0.548000000  0.935000000
                 H -0.214000000 -0.548000000 -0.935000000
                 H  0.337000000  3.142000000  0.000000000""",
}

METHODS = {
    "hf": (None, 1.0),
    "pbe": ("pbe", 0.0),
    "pbe0": ("pbe0", 0.25),
    "bhhlyp": ("bhandhlyp", 0.5),
}


def _mol(name: str, basis: str, translation_bohr=(0.0, 0.0, 0.0)):
    mol = gto.M(atom=GEOMS_ANG[name], basis=basis, unit="Angstrom", charge=0, spin=0, verbose=0)
    if any(abs(x) > 0 for x in translation_bohr):
        coords = mol.atom_coords(unit="Bohr") + np.asarray(translation_bohr, dtype=float)
        atoms = [(mol.atom_symbol(i), coords[i]) for i in range(mol.natm)]
        mol = gto.M(atom=atoms, basis=basis, unit="Bohr", charge=0, spin=0, verbose=0)
    return mol


def _coupled_mo1(mf, h1, s1, max_cycle=100, conv_tol=1e-10):
    coeff, energy, occ = mf.mo_coeff, mf.mo_energy, mf.mo_occ
    nmo = len(energy)
    occ_mask = occ > 0
    vir_mask = occ == 0
    e_ai = 1.0 / (energy[vir_mask].reshape(-1, 1) - energy[occ_mask])
    vind = nmr_rhf.gen_vind(mf, coeff, occ)
    hs = h1 - s1 * energy[occ_mask]
    mo1 = hs.copy()
    mo1[:, vir_mask, :] = -hs[:, vir_mask, :] * e_ai
    mo1[:, occ_mask, :] = -0.5 * s1[:, occ_mask, :]
    for _ in range(max_cycle):
        v1 = vind(mo1.ravel()).reshape(3, nmo, occ_mask.sum())
        new = mo1.copy()
        new[:, vir_mask, :] = -(hs[:, vir_mask, :] + v1[:, vir_mask, :]) * e_ai
        new[:, occ_mask, :] = -0.5 * s1[:, occ_mask, :]
        if np.max(np.abs(new - mo1)) < conv_tol:
            return new
        mo1 = new
    raise RuntimeError("PySCF NMR CPHF reference did not converge")


def _pyscf_shielding(mol, method: str, gauge: str, origin_bohr=(0.0, 0.0, 0.0)):
    xc, cx = METHODS[method]
    t0 = time.perf_counter()
    with PeakMemorySampler() as mem:
        if xc is None:
            mf = scf.RHF(mol).run(verbose=0)
            cls = nmr_rhf.NMR
        else:
            mf = dft.RKS(mol)
            mf.xc = xc
            mf.run(verbose=0)
            cls = nmr_rks.NMR
        nmrobj = cls(mf)
        gauge_orig = None if gauge == "giao" else np.asarray(origin_bohr, dtype=float)
        coeff = mf.mo_coeff
        occ_mask = mf.mo_occ > 0
        h1ao = nmr_rhf.make_h10(nmrobj.mol, mf.make_rdm1(), gauge_orig=gauge_orig)
        s1ao = nmr_rhf.make_s10(nmrobj.mol, gauge_orig=gauge_orig)
        h1 = np.einsum("pi,xpq,qj->xij", coeff, h1ao, coeff[:, occ_mask])
        s1 = np.einsum("pi,xpq,qj->xij", coeff, s1ao, coeff[:, occ_mask])
        mo_coupled = _coupled_mo1(mf, h1, s1)
        tensor = (nmrobj.dia(gauge_orig=gauge_orig) + nmrobj.para(mo10=mo_coupled)[0]) * UNIT
    wall = time.perf_counter() - t0
    iso = np.asarray([float(np.trace(tensor[i]) / 3.0) for i in range(mol.natm)])
    return {
        "tensor_components_ppm": tensor.tolist(),
        "isotropic_shielding_ppm": iso.tolist(),
        "wall_time_seconds": wall,
        "peak_memory_mb": mem.peak_mb,
        "exact_exchange_fraction": cx,
    }


def _openqp_runtime_error() -> str | None:
    try:
        import oqp  # noqa: F401 - import initializes the installed package root
    except Exception as exc:
        return f"PyOpenQP is not importable from the active Python environment: {type(exc).__name__}: {exc}"

    root = Path(os.environ.get("OPENQP_ROOT", ""))
    if not root:
        return "Installed PyOpenQP did not initialize OPENQP_ROOT"
    libdir = root / "lib"
    if not ((libdir / "liboqp.dylib").exists() or (libdir / "liboqp.so").exists()):
        return f"installed PyOpenQP root {root} does not contain lib/liboqp.dylib or lib/liboqp.so"
    return None


def _run_openqp(name: str, basis: str, method: str, gauge: str):
    runtime_error = _openqp_runtime_error()
    if runtime_error:
        raise RuntimeError(runtime_error)
    xc, _ = METHODS[method]
    functional = f"functional={method}\n" if xc else ""
    dftgrid = "[dftgrid]\nrad_type=becke\n\n" if xc else ""
    system = "\n".join(
        f"   {gto.charge(sym):d} {float(x): .10f} {float(y): .10f} {float(z): .10f}"
        for sym, (x, y, z) in gto.format_atom(GEOMS_ANG[name], unit=1)
    )
    text = (
        f"[input]\nsystem=\n{system}\ncharge=0\nruntype=energy\n{functional}basis={basis}\nmethod=hf\n\n"
        f"[guess]\ntype=huckel\n\n[scf]\nmultiplicity=1\ntype=rhf\n\n{dftgrid}"
        f"[properties]\nscf_prop=nmr\nnmr_gauge={gauge}\n"
    )
    with tempfile.TemporaryDirectory() as wd:
        inp = Path(wd) / f"{name}_{method}_{gauge}.inp"
        inp.write_text(text)
        env = dict(os.environ)
        t0 = time.perf_counter()
        if psutil is not None:
            proc_handle = psutil.Popen(
                [sys.executable, "-m", "oqp.pyoqp", str(inp)],
                cwd=wd,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            with PeakMemorySampler(proc_handle.pid) as mem:
                stdout, stderr = proc_handle.communicate()
            proc = subprocess.CompletedProcess(proc_handle.args, proc_handle.returncode, stdout, stderr)
            peak_memory_mb = mem.peak_mb
        else:
            with _child_memory_delta() as mem:
                proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(inp)], cwd=wd, env=env, capture_output=True, text=True)
            peak_memory_mb = mem["peak_memory_mb"]
        wall = time.perf_counter() - t0
        log = inp.with_suffix(".log")
        stdout_stderr = proc.stdout + proc.stderr + (log.read_text() if log.exists() else "")
        if "NotImplementedError" in stdout_stderr or "not yet implemented" in stdout_stderr:
            return {"status": "gated", "reason": "OpenQP native GIAO path is still NotImplemented", "wall_time_seconds": wall, "peak_memory_mb": peak_memory_mb}
        if proc.returncode != 0:
            return {"status": "failed", "returncode": proc.returncode, "tail": stdout_stderr[-2000:], "wall_time_seconds": wall, "peak_memory_mb": peak_memory_mb}
        rows = []
        for line in stdout_stderr.splitlines():
            m = re.match(r"\s*(\d+)\s+[\d.]+" + r"\s+(-?\d+\.\d+)" * 5 + r"\s*$", line)
            if m:
                vals = [float(m.group(j)) for j in range(2, 7)]
                rows.append(vals[4])
        return {"status": "ok", "isotropic_shielding_ppm": rows, "wall_time_seconds": wall, "peak_memory_mb": peak_memory_mb}


def _write_outputs(records, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    json_path = outdir / "nmr_giao_benchmark_results.json"
    csv_path = outdir / "nmr_giao_benchmark_results.csv"
    md_path = outdir / "nmr_giao_benchmark_results.md"
    json_path.write_text(json.dumps(records, indent=2) + "\n")
    fields = [
        "system",
        "basis",
        "method",
        "backend",
        "gauge",
        "status",
        "origin_dependence_ppm",
        "sigma_iso_ppm",
        "pyscf_reference_iso_ppm",
        "openqp_minus_pyscf_delta_ppm",
        "wall_time_seconds",
        "peak_memory_mb",
        "notes",
    ]
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for r in records["rows"]:
            writer.writerow({k: r.get(k, "") for k in fields})
    lines = [
        "# NMR CGO/GIAO benchmark matrix",
        "",
        f"Status: {records['status']}",
        "",
        "|System|Basis|Method|Backend|Gauge|Status|Origin dependence/ppm|Sigma iso/ppm|PySCF reference iso/ppm|OpenQP-PySCF delta/ppm|Wall/s|Peak memory/MiB|Notes|",
        "|---|---:|---|---|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in records["rows"]:
        lines.append("|{system}|{basis}|{method}|{backend}|{gauge}|{status}|{origin_dependence_ppm}|{sigma_iso_ppm}|{pyscf_reference_iso_ppm}|{openqp_minus_pyscf_delta_ppm}|{wall_time_seconds}|{peak_memory_mb}|{notes}|".format(**{k: r.get(k, "") for k in fields}))
    md_path.write_text("\n".join(lines) + "\n")
    return json_path, csv_path, md_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true", help="run every matrix case; default runs a quick subset")
    parser.add_argument("--run-openqp", action="store_true", help="also run OpenQP CGO/GIAO through the normal driver")
    parser.add_argument("--check-openqp-runtime", action="store_true", help="check the active pip-installed OpenQP runtime and exit")
    parser.add_argument("--outdir", default=str(OUTDIR))
    args = parser.parse_args()

    if args.run_openqp or args.check_openqp_runtime:
        runtime_error = _openqp_runtime_error()
        if runtime_error:
            requested = "--check-openqp-runtime" if args.check_openqp_runtime else "--run-openqp"
            raise SystemExit(f"{requested} requested but pip-installed OpenQP runtime is unavailable: {runtime_error}")
        if args.check_openqp_runtime:
            print("pip-installed OpenQP runtime is available")
            print(f"OPENQP_ROOT={os.environ.get('OPENQP_ROOT', '')}")
            return

    _ensure_pyscf_imports()

    matrix = json.loads(MATRIX.read_text())
    translations = matrix["gauge_origin_translations_bohr"]
    cases = []
    sections = matrix["small_regression_systems"] + (matrix["origin_dependence_systems"] if args.all else [])
    if not args.all:
        sections = [matrix["small_regression_systems"][0], matrix["origin_dependence_systems"][0]]
    for item in sections:
        for method in item["methods"]:
            cases.append((item["name"], item["basis"], method.lower()))

    rows = []
    for name, basis, method in cases:
        if method not in METHODS:
            continue
        pyscf_reference = {}
        for gauge in ("cgo", "giao"):
            vals = []
            tensors = None
            wall = 0.0
            peak_memory_mb = 0.0
            status = "ok"
            notes = "PySCF trusted reference"
            for trans in translations:
                try:
                    mol = _mol(name, basis, translation_bohr=trans if gauge == "giao" else (0.0, 0.0, 0.0))
                    result = _pyscf_shielding(mol, method, gauge, origin_bohr=trans)
                    vals.append(np.asarray(result["isotropic_shielding_ppm"]))
                    tensors = result["tensor_components_ppm"]
                    wall += result["wall_time_seconds"]
                    peak_memory_mb = max(peak_memory_mb, float(result.get("peak_memory_mb", 0.0)))
                except Exception as exc:
                    status = "failed"
                    notes = type(exc).__name__ + ": " + str(exc)[:160]
                    break
            if vals:
                arr = np.vstack(vals)
                max_delta = float(np.max(np.ptp(arr, axis=0)))
                sigma_iso_ppm = ";".join(f"{x:.6f}" for x in arr[0])
                pyscf_reference[gauge] = arr[0]
            else:
                max_delta = ""
                sigma_iso_ppm = ""
            rows.append({
                "system": name, "basis": basis, "method": method, "backend": "pyscf", "gauge": gauge,
                "status": status, "origin_dependence_ppm": f"{max_delta:.6f}" if isinstance(max_delta, float) else max_delta,
                "sigma_iso_ppm": sigma_iso_ppm,
                "pyscf_reference_iso_ppm": sigma_iso_ppm,
                "openqp_minus_pyscf_delta_ppm": "",
                "wall_time_seconds": f"{wall:.3f}", "peak_memory_mb": f"{peak_memory_mb:.1f}", "notes": notes,
                "tensor_components_ppm": tensors,
            })
        if args.run_openqp:
            for gauge in ("cgo", "giao"):
                r = _run_openqp(name, basis, method, gauge)
                iso = r.get("isotropic_shielding_ppm") or []
                ref = pyscf_reference.get(gauge)
                delta = ""
                if iso and ref is not None:
                    assert np is not None
                    delta = f"{float(np.max(np.abs(np.asarray(iso) - ref))):.6f}"
                rows.append({
                    "system": name, "basis": basis, "method": method, "backend": "openqp", "gauge": gauge,
                    "status": r["status"],
                    "origin_dependence_ppm": "pending" if gauge == "giao" and r["status"] == "ok" else "",
                    "sigma_iso_ppm": ";".join(f"{x:.6f}" for x in iso),
                    "pyscf_reference_iso_ppm": ";".join(f"{x:.6f}" for x in ref) if ref is not None else "",
                    "openqp_minus_pyscf_delta_ppm": delta,
                    "wall_time_seconds": f"{r.get('wall_time_seconds', 0.0):.3f}", "peak_memory_mb": f"{r.get('peak_memory_mb', 0.0):.1f}",
                    "notes": r.get("reason", r.get("tail", ""))[:160],
                })
    records = {
        "status": "reference_ready_openqp_giao_gated",
        "metric_definition": "Origin-dependence, shielding magnitude, and OpenQP-minus-PySCF agreement are separate fields; no row mixes these quantities in one generic metric column.",
        "matrix_file": str(MATRIX),
        "metrics": matrix["metrics"],
        "origin_translations_bohr": translations,
        "rows": rows,
    }
    paths = _write_outputs(records, Path(args.outdir))
    for p in paths:
        print(p)


if __name__ == "__main__":
    main()
