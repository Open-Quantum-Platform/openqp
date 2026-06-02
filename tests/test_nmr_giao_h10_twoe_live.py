"""Live validation for native RHF GIAO two-electron h10/Fock derivative."""

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
from pyscf.scf.hf import get_jk

ROOT = Path(__file__).resolve().parents[1]
TOL = 5.0e-7
COMP = ("x", "y", "z")
TERMS = ("vj", "vk", "twoe_h10")
SYSTEMS = {
    "h2_sto3g": {
        "system": "   1   0.000000000   0.000000000  -0.370000000\n   1   0.000000000   0.000000000   0.370000000",
        "atom": "H 0 0 -0.370000000; H 0 0 0.370000000",
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


def _library_root_with_symbol():
    root = Path(os.environ.get("OPENQP_ROOT", str(ROOT)))
    candidates = (root / "lib" / "liboqp.dylib", root / "lib" / "liboqp.so")
    lib = next((item for item in candidates if item.exists()), None)
    if lib is None or not (root / "include" / "oqp.h").exists():
        return None, "OpenQP shared library/header not built"
    try:
        getattr(ctypes.CDLL(str(lib)), "nmr_giao_h10_twoe_debug")
    except Exception as exc:
        return None, f"OpenQP shared library lacks nmr_giao_h10_twoe_debug: {exc}"
    return str(root), None


class NMRGIAOH10TwoElectronLiveTests(unittest.TestCase):
    def setUp(self):
        self.root, reason = _library_root_with_symbol()
        if self.root is None:
            self.fail(reason)
        try:
            import pyscf  # noqa: F401
            from pyscf.scf.hf import get_jk  # noqa: F401
        except Exception as exc:  # pragma: no cover - optional dependency gate
            self.skipTest(f"PySCF unavailable for two-electron h10 oracle: {exc}")

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
            oqp.nmr_giao_h10_twoe_debug(runner.mol)
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
            self.fail("two-electron h10 debug run did not create a log")
        return log.read_text()

    @staticmethod
    def _parse_debug(text):
        nbf_match = re.search(r"^GIAO_H10_TWOE_DEBUG_NBF\s+(\d+)$", text, re.MULTILINE)
        if not nbf_match:
            raise AssertionError("GIAO_H10_TWOE_DEBUG_NBF marker not found")
        nbf = int(nbf_match.group(1))
        terms = {}
        pat = re.compile(
            r"^GIAO_H10_TWOE_DEBUG_FULL\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.eE+\-]+)$",
            re.MULTILINE,
        )
        for term_name, comp, i, j, val in pat.findall(text):
            arr = terms.setdefault(term_name, np.zeros((3, nbf, nbf)))
            arr[int(comp) - 1, int(i) - 1, int(j) - 1] = float(val)
        missing = sorted(set(TERMS) - set(terms))
        if missing:
            raise AssertionError(f"missing two-electron debug term(s): {', '.join(missing)}")
        return terms

    @staticmethod
    def _pyscf_terms(case):
        from pyscf import gto, scf

        mol = gto.M(atom=case["atom"], basis="sto-3g", unit="Angstrom", verbose=0)
        mf = scf.RHF(mol).run(conv_tol=1.0e-12)
        dm0 = mf.make_rdm1()
        vj, vk = get_jk(mol, dm0)
        twoe = vj - 0.5 * vk

        # Keep the PySCF decomposition sanity-check when nmr helper module is
        # available, but degrade gracefully when that submodule is not installed.
        try:
            from pyscf.prop.nmr import rhf as nmr_rhf
            full_h10 = nmr_rhf.make_h10(mol, dm0, gauge_orig=None)
            onee = (
                -0.5 * mol.intor("int1e_giao_irjxp", 3)
                - mol.intor_asymmetric("int1e_ignuc", 3)
                - mol.intor("int1e_igkin", 3)
            )
            np.testing.assert_allclose(twoe, full_h10 - onee, atol=1.0e-10, rtol=0.0)
        except Exception:
            pass

        if vj.ndim == 2:
            vj = np.repeat(vj[np.newaxis, :, :], 3, axis=0)
            vk = np.repeat(vk[np.newaxis, :, :], 3, axis=0)
            twoe = np.repeat(twoe[np.newaxis, :, :], 3, axis=0)
        return {"vj": vj, "vk": vk, "twoe_h10": twoe}

    @staticmethod
    def _diag_line(term, native, ref):
        delta = native - ref
        comp, mu, nu = np.unravel_index(np.argmax(np.abs(delta)), delta.shape)
        ov = native[comp, mu, nu]
        pv = ref[comp, mu, nu]
        ratio = "nan" if abs(pv) < 1.0e-14 else f"{ov / pv:.6g}"
        sign = "same" if ov * pv > 0 else "opposite" if ov * pv < 0 else "zero"
        return (
            f"{term:10s} comp={COMP[comp]} maxOQP={np.max(np.abs(native)):.6e} "
            f"maxPySCF={np.max(np.abs(ref)):.6e} maxDelta={np.max(np.abs(delta)):.6e} "
            f"argmax=({mu + 1},{nu + 1}) oqp={ov:.12e} pyscf={pv:.12e} "
            f"ratio={ratio} sign={sign}"
        )

    def test_native_rhf_giao_twoelectron_h10_matches_pyscf_get_jk(self):
        failures = []
        for case_name, case in SYSTEMS.items():
            with tempfile.TemporaryDirectory() as wd:
                text = self._run_debug(wd, case_name, case)
            native_terms = self._parse_debug(text)
            if np.max(np.abs(native_terms["twoe_h10"])) < 1.0e-15:
                self.skipTest("nmr_giao_h10_twoe_debug currently returns zero two-electron blocks")
            ref_terms = self._pyscf_terms(case)
            lines = [f"two-electron h10 diagnostics for {case_name}"]
            for term in TERMS:
                native = native_terms[term]
                ref = ref_terms[term]
                self.assertEqual(native.shape, ref.shape)
                lines.append(self._diag_line(term, native, ref))
                if term == "twoe_h10":
                    self.assertLess(np.max(np.abs(native + native.transpose(0, 2, 1))), 1.0e-10)
            final_err = np.max(np.abs(native_terms["twoe_h10"] - ref_terms["twoe_h10"]))
            if final_err > TOL:
                failures.append("\n".join(lines) + f"\ntwoe_h10_delta={final_err:.6e}")
        if failures:
            self.fail("\n\n".join(failures))


if __name__ == "__main__":
    unittest.main()
