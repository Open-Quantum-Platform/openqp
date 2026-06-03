"""Live validation for native one-electron GIAO h10 debug kernels."""

import ctypes
import os
import re
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TOL = 5.0e-7
COMP = ("x", "y", "z")
TERM_NAMES = (
    "raw_irjxp",
    "-0.5_irjxp",
    "raw_igkin",
    "-igkin",
    "raw_ignuc",
    "ignuc_asym",
    "-ignuc_asym",
    "final_h10",
)
SYSTEMS = {
    "h2_sto3g": {
        "system": "   1   0.000000000   0.000000000  -0.370000000\n   1   0.000000000   0.000000000   0.370000000",
        "atom": "H 0 0 -0.370000000; H 0 0 0.370000000",
    },
    "he_sto3g": {
        "system": "   2   0.000000000   0.000000000   0.000000000",
        "atom": "He 0 0 0",
    },
    "h2o_sto3g": {
        "system": "   8   0.000000000   0.000000000  -0.041061554\n   1  -0.533194329   0.533194329  -0.614469223\n   1   0.533194329  -0.533194329  -0.614469223",
        "atom": "O 0 0 -0.041061554; H -0.533194329 0.533194329 -0.614469223; H 0.533194329 -0.533194329 -0.614469223",
    },
}

INPUT_TEMPLATE = """\
[input]
system=
{system}
charge=0
runtype=energy
basis=sto-3g
method=hf

[guess]
type=huckel

[scf]
multiplicity=1
type=rhf

[properties]
scf_prop=
"""


def _oqp_root():
    root = Path(os.environ.get("OPENQP_ROOT", str(ROOT)))
    libdir = root / "lib"
    lib = libdir / "liboqp.dylib"
    if not lib.exists():
        lib = libdir / "liboqp.so"
    if not lib.exists():
        return None
    if not (root / "include" / "oqp.h").exists():
        return None
    try:
        getattr(ctypes.CDLL(str(lib)), "nmr_giao_h10_debug")
    except Exception:
        return None
    return str(root)


class NMRGIAOH10LiveTests(unittest.TestCase):
    def setUp(self):
        root = _oqp_root()
        if root is None:
            self.skipTest("OpenQP shared library/header not built with nmr_giao_h10_debug")
        self.root = root
        try:
            import pyscf  # noqa: F401
        except Exception as exc:  # pragma: no cover - optional dependency gate
            self.skipTest(f"PySCF unavailable for h10 oracle: {exc}")

    def _run_debug(self, workdir, case_name, case):
        inp = Path(workdir) / f"{case_name}.inp"
        log = Path(workdir) / f"{case_name}.log"
        inp.write_text(INPUT_TEMPLATE.format(system=case["system"]))
        script = Path(workdir) / f"run_{case_name}_debug.py"
        script.write_text(textwrap.dedent(f"""
            import oqp
            from oqp.pyoqp import Runner
            runner = Runner(project={case_name!r}, input_file={str(inp)!r}, log={str(log)!r}, silent=1, usempi=False)
            runner.run()
            oqp.nmr_giao_h10_debug(runner.mol)
        """))
        env = dict(os.environ)
        env["OPENQP_ROOT"] = self.root
        # Use this checkout's Python wrapper while allowing it to load the
        # selected compiled OpenQP runtime from OPENQP_ROOT.  This keeps live
        # tests from accidentally exercising an older pip-installed wrapper
        # when OPENQP_ROOT points at site-packages/oqp.
        root_parent = str(Path(self.root).parent)
        source_parent = str(ROOT / "pyoqp")
        env["PYTHONPATH"] = os.pathsep.join(
            item for item in (source_parent, root_parent, env.get("PYTHONPATH", "")) if item
        )
        proc = subprocess.run([sys.executable, str(script)], cwd=workdir, env=env,
                              capture_output=True, text=True, timeout=120)
        if proc.returncode != 0:
            self.fail(proc.stdout[-2000:] + proc.stderr[-4000:])
        if not log.exists():
            self.fail("h10 debug run did not create a log")
        return log.read_text()

    @staticmethod
    def _parse_h10(text):
        nbf_match = re.search(r"^GIAO_H10_DEBUG_NBF\s+(\d+)$", text, re.MULTILINE)
        if not nbf_match:
            raise AssertionError("GIAO_H10_DEBUG_NBF marker not found")
        nbf = int(nbf_match.group(1))
        packed = np.zeros((3, nbf, nbf))
        pat = re.compile(r"^GIAO_H10_DEBUG_PACKED\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.eE+\-]+)$", re.MULTILINE)
        count = 0
        for comp, i, j, val in pat.findall(text):
            c = int(comp) - 1
            ii = int(i) - 1
            jj = int(j) - 1
            v = float(val)
            packed[c, ii, jj] = v
            packed[c, jj, ii] = -v
            count += 1
        expected_count = 3 * nbf * (nbf + 1) // 2
        if count != expected_count:
            raise AssertionError(f"expected {expected_count} h10 packed records, found {count}")

        terms = {}
        full_pat = re.compile(
            r"^GIAO_H10_DEBUG_FULL\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.eE+\-]+)$",
            re.MULTILINE,
        )
        for term_idx, term_name, comp, i, j, val in full_pat.findall(text):
            name = term_name
            arr = terms.setdefault(name, np.zeros((3, nbf, nbf)))
            arr[int(comp) - 1, int(i) - 1, int(j) - 1] = float(val)
        return packed, terms

    @staticmethod
    def _pyscf_terms(case):
        from pyscf import gto
        mol = gto.M(atom=case["atom"], basis="sto-3g", unit="Angstrom", verbose=0)
        raw_irjxp = mol.intor("int1e_giao_irjxp", 3)
        raw_igkin = mol.intor("int1e_igkin", 3)
        raw_ignuc = mol.intor_asymmetric("int1e_ignuc", 3)
        ignuc_asym = raw_ignuc
        return {
            "raw_irjxp": raw_irjxp,
            "-0.5_irjxp": -0.5 * raw_irjxp,
            "raw_igkin": raw_igkin,
            "-igkin": -raw_igkin,
            "raw_ignuc": raw_ignuc,
            "ignuc_asym": ignuc_asym,
            "-ignuc_asym": -ignuc_asym,
            "final_h10": -0.5 * raw_irjxp - raw_igkin - ignuc_asym,
        }

    @staticmethod
    def _diag_line(term, native, ref):
        delta = native - ref
        comp, mu, nu = np.unravel_index(np.argmax(np.abs(delta)), delta.shape)
        ov = native[comp, mu, nu]
        pv = ref[comp, mu, nu]
        ratio = "nan" if abs(pv) < 1.0e-14 else f"{ov / pv:.6g}"
        sign = "same" if ov * pv > 0 else "opposite" if ov * pv < 0 else "zero"
        return (
            f"{term:14s} comp={COMP[comp]} maxOQP={np.max(np.abs(native)):.6e} "
            f"maxPySCF={np.max(np.abs(ref)):.6e} maxDelta={np.max(np.abs(delta)):.6e} "
            f"argmax=({mu + 1},{nu + 1}) oqp={ov:.12e} pyscf={pv:.12e} "
            f"ratio={ratio} sign={sign}"
        )

    def _case_diagnostics(self, case_name, native_terms, ref_terms):
        lines = [f"term-level h10 diagnostics for {case_name}"]
        first_wrong = None
        for term in TERM_NAMES:
            if term not in native_terms:
                lines.append(f"{term:14s} missing from OpenQP debug output")
                if first_wrong is None:
                    first_wrong = term
                continue
            line = self._diag_line(term, native_terms[term], ref_terms[term])
            lines.append(line)
            if first_wrong is None and np.max(np.abs(native_terms[term] - ref_terms[term])) > TOL:
                first_wrong = term
        lines.append(f"first_wrong_term={first_wrong}")
        return "\n".join(lines), first_wrong

    def test_native_h10_matches_pyscf_oneelectron_oracle(self):
        failures = []
        for case_name, case in SYSTEMS.items():
            with tempfile.TemporaryDirectory() as wd:
                text = self._run_debug(wd, case_name, case)
            packed, native_terms = self._parse_h10(text)
            ref_terms = self._pyscf_terms(case)
            self.assertEqual(packed.shape, ref_terms["final_h10"].shape)
            diagnostics, _first_wrong = self._case_diagnostics(case_name, native_terms, ref_terms)
            self.assertLess(np.max(np.abs(packed + packed.transpose(0, 2, 1))), 1.0e-10)
            final_err = np.max(np.abs(packed - ref_terms["final_h10"]))
            if final_err > TOL:
                failures.append(diagnostics + f"\npacked_final_delta={final_err:.6e}")
        if failures:
            self.fail("\n\n".join(failures))


if __name__ == "__main__":
    unittest.main()
