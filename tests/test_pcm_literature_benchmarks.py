"""ddX-enabled scientific validation of the PCM energy path.

This harness drives the **Fortran/OpenQP** PCM path (`pcm_enabled` ->
`add_pcm_reaction_field` -> `E%e_pcm` -> ddX C adapter) through real
single-point inputs and compares the computed solvation energy against the
reference targets in ``tests/data/pcm_literature_benchmarks.json``
(see ``docs/solvent_pcm_literature_benchmarks.md``).

Skip semantics (so the gate never produces false confidence):

* The whole module skips when ``lib/liboqp`` is not built.
* A PCM-on benchmark skips when OpenQP was built without ddX
  (``OQP_ENABLE_DDX``): the C adapter reports this and the run aborts cleanly.
* A benchmark whose ``reference_value`` is still ``null`` / ``status`` is not
  ``"verified"`` is skipped even when ddX is available, with a message pointing
  to the procedure for populating it. This avoids comparing against a fabricated
  number.

Python here only builds inputs and parses output. It does **not** compute the
reaction field, Fock contribution, or solvation energy.
"""

import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data" / "pcm_literature_benchmarks.json"


def _lib_suffix() -> str:
    return {"Windows": "dll", "Linux": "so", "Darwin": "dylib"}.get(
        platform.uname()[0], "so"
    )


LIBOQP = ROOT / "lib" / f"liboqp.{_lib_suffix()}"


def _load_benchmarks():
    data = json.loads(DATA.read_text(encoding="utf-8"))
    return data["benchmarks"]


def _load_trusted_regressions():
    data = json.loads(DATA.read_text(encoding="utf-8"))
    return data.get("trusted_reference_regressions", [])


def _find_ddx_adapter_smoke():
    """Locate the ddX-enabled `oqp_ddx_adapter_smoke` CTest binary, if built.

    The binary only exists in a `-DENABLE_DDX=ON` build, so this is the natural
    skip guard for the Tier-1 ddX-adapter regression. A direct path may be given
    via OQP_DDX_ADAPTER_SMOKE; otherwise OQP_DDX_BUILD_DIR is searched.
    """
    direct = os.environ.get("OQP_DDX_ADAPTER_SMOKE")
    if direct and Path(direct).is_file() and os.access(direct, os.X_OK):
        return Path(direct)
    build_dir = os.environ.get("OQP_DDX_BUILD_DIR")
    if build_dir:
        for cand in (
            Path(build_dir) / "oqp_ddx_adapter_smoke",
            Path(build_dir) / "tests" / "oqp_ddx_adapter_smoke",
            Path(build_dir) / "bin" / "oqp_ddx_adapter_smoke",
        ):
            if cand.is_file() and os.access(cand, os.X_OK):
                return cand
    return None


def _input_text(bench, *, pcm_on: bool) -> str:
    system = "\n".join(bench["system"])
    pcm = ""
    if pcm_on:
        pcm = (
            "\n[pcm]\n"
            "enabled=true\n"
            "backend=ddx\n"
            "mode=reference_scf\n"
            "model=ddpcm\n"
            f"epsilon={bench['epsilon']}\n"
        )
    else:
        pcm = "\n[pcm]\nenabled=false\n"
    return (
        "[input]\n"
        "system=\n"
        f"{system}\n"
        f"charge={bench['charge']}\n"
        "runtype=energy\n"
        f"basis={bench['basis']}\n"
        f"method={bench['method']}\n"
        "\n"
        "[guess]\n"
        "type=huckel\n"
        "\n"
        "[scf]\n"
        f"multiplicity={bench['multiplicity']}\n"
        f"type={bench['scf_type']}\n"
        f"{pcm}"
    )


def _run(text: str):
    with tempfile.TemporaryDirectory() as d:
        inp = Path(d) / "bench.inp"
        inp.write_text(text)
        env = dict(os.environ)
        env["OPENQP_ROOT"] = str(ROOT)
        env["PYTHONPATH"] = str(ROOT / "pyoqp") + os.pathsep + env.get("PYTHONPATH", "")
        env["OMP_NUM_THREADS"] = "1"
        proc = subprocess.run(
            [sys.executable, "-m", "oqp.pyoqp", str(inp)],
            env=env,
            capture_output=True,
            text=True,
            cwd=d,
        )
        log_path = inp.with_suffix(".log")
        log = log_path.read_text() if log_path.is_file() else ""
        return proc, "\n".join([log, proc.stdout, proc.stderr])


def _total_energy(log: str):
    m = re.findall(r"TOTAL energy\s*=\s*(-?\d+\.\d+)", log)
    return float(m[-1]) if m else None


def _pcm_energy(log: str):
    m = re.findall(r"PCM solvent energy.*?=\s*(-?\d+\.\d+)", log)
    return float(m[-1]) if m else None


def _ddx_unavailable(log: str) -> bool:
    return "OQP_ENABLE_DDX" in log or "built without" in log


@unittest.skipUnless(
    LIBOQP.is_file(), f"liboqp not built at {LIBOQP}; build OpenQP first"
)
class PcmLiteratureBenchmarks(unittest.TestCase):
    """Per-benchmark validation; methods are attached dynamically below."""


def _make_vacuum_test(bench):
    def test(self):
        proc, log = _run(_input_text(bench, pcm_on=False))
        self.assertEqual(proc.returncode, 0, log)
        e = _total_energy(log)
        self.assertIsNotNone(e, log)
        # PCM-off must reproduce the established vacuum total energy exactly.
        self.assertAlmostEqual(e, bench["vacuum_total_energy"], places=7)
        # And emit no PCM term.
        self.assertIsNone(_pcm_energy(log), log)

    return test


def _make_pcm_test(bench):
    def test(self):
        proc, log = _run(_input_text(bench, pcm_on=True))

        if _ddx_unavailable(log):
            self.skipTest(
                "OpenQP built without ddX (OQP_ENABLE_DDX); the Fortran PCM "
                "path is wired but not executable here."
            )

        # ddX is available: the Fortran PCM path must run and converge (Q4),
        # and report a finite solvation energy (Q1) and total (Q3).
        self.assertEqual(proc.returncode, 0, log)
        e_pcm = _pcm_energy(log)
        self.assertIsNotNone(e_pcm, "PCM energy term not reported:\n" + log)
        e_tot = _total_energy(log)
        self.assertIsNotNone(e_tot, log)

        ref = bench.get("reference_value")
        if ref is None or bench.get("status") != "verified":
            self.skipTest(
                f"benchmark '{bench['id']}' has no verified reference yet "
                f"(status={bench.get('status')!r}); populate "
                "tests/data/pcm_literature_benchmarks.json per "
                "docs/solvent_pcm_literature_benchmarks.md before this is a "
                "pass/fail gate. Observed e_pcm=" + repr(e_pcm)
            )

        # Real scientific gate: match the reference within the stated tolerance.
        tol = float(bench["tolerance"])
        self.assertLessEqual(
            abs(e_pcm - float(ref)),
            tol,
            f"{bench['id']}: e_pcm={e_pcm} vs reference={ref} "
            f"exceeds tolerance {tol}\n{log}",
        )

    return test


def _attach_tests():
    for bench in _load_benchmarks():
        bid = re.sub(r"[^0-9a-zA-Z_]", "_", bench["id"])
        setattr(PcmLiteratureBenchmarks, f"test_pcm_{bid}", _make_pcm_test(bench))
        if bench.get("vacuum_total_energy") is not None:
            setattr(
                PcmLiteratureBenchmarks,
                f"test_vacuum_{bid}",
                _make_vacuum_test(bench),
            )


_attach_tests()


class DdxAdapterRegressionGate(unittest.TestCase):
    """Tier-1: ddX-generated trusted-reference regressions via the C adapter.

    These run the Fortran/C `oqp_ddx_adapter_smoke` driver (real ddX) and compare
    its printed energy to the verified `trusted_reference_regression` value in the
    benchmark table. They validate the OpenQP<->ddX adapter numerics only; they do
    NOT exercise the QM SCF coupling (electronic psi, q_cav-into-Fock,
    phi_cav-from-density), which the Tier-2 QM benchmarks cover. Python only runs
    the binary and parses its output.
    """

    def _run_point_charge(self, bench):
        binary = _find_ddx_adapter_smoke()
        if binary is None:
            self.skipTest(
                "ddX-enabled `oqp_ddx_adapter_smoke` not built; set "
                "OQP_DDX_BUILD_DIR (or OQP_DDX_ADAPTER_SMOKE) to a "
                "cmake -DENABLE_DDX=ON build to run this regression."
            )
        proc = subprocess.run([str(binary)], capture_output=True, text=True)
        out = "\n".join([proc.stdout, proc.stderr])
        self.assertEqual(proc.returncode, 0, out)
        # The point-charge lifecycle prints `energy=<float>` (not the later
        # `explicit_energy=`/`reaction_field_` lines); the lookbehind avoids them.
        m = re.search(r"(?<![A-Za-z_])energy=(-?\d+\.\d+)", out)
        self.assertIsNotNone(m, "adapter smoke did not print energy=:\n" + out)
        energy = float(m.group(1))
        ref = float(bench["reference_value"])
        tol = float(bench["tolerance"])
        self.assertLessEqual(
            abs(energy - ref),
            tol,
            f"{bench['id']}: ddX esolv={energy} vs trusted reference={ref} "
            f"exceeds tolerance {tol}\n{out}",
        )


def _attach_regression_tests():
    for bench in _load_trusted_regressions():
        if bench.get("reference_value") is None or bench.get("status") != "verified":
            continue
        bid = re.sub(r"[^0-9a-zA-Z_]", "_", bench["id"])

        def test(self, bench=bench):
            self._run_point_charge(bench)

        setattr(DdxAdapterRegressionGate, f"test_regression_{bid}", test)


_attach_regression_tests()


if __name__ == "__main__":
    unittest.main()
