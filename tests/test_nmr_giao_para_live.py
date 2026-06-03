"""Live validation for the native GIAO paramagnetic NMR shielding.

Runs the native GIAO shielding debug emitter (``nmr_giao_shielding_debug``) on
H2O/STO-3G for HF/BHHLYP/PBE0/PBE and compares the per-atom isotropic
paramagnetic shielding (uncoupled and coupled) against the committed PySCF GIAO
oracle (``tests/fixtures/nmr/pyscf_giao_reference.json``).

This validates the native London-orbital first-order Hamiltonian/overlap
assembly (h10 one- and two-electron + S10), the GIAO CPHF/CPKS first-order
solve (exchange-only coupled response scaled by the exact-exchange fraction),
and the PSO contraction.  The diamagnetic GIAO term (second-order GIAO
integrals) is a separate, later checkpoint and is not asserted here; the
production ``nmr_gauge=giao`` route therefore remains gated.

PySCF is NOT required at run time: the oracle JSON is committed.  The test runs
the OpenQP driver in a subprocess and skips gracefully if the shared library is
not built.
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "nmr" / "pyscf_giao_reference.json"

GEOM = """\
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223"""

# (input label, OpenQP functional keyword, oracle fixture key)
CASES = [("HF", None, "HF"), ("BHHLYP", "bhhlyp", "bhandhlyp"),
         ("PBE0", "pbe0", "pbe0"), ("PBE", "pbe", "pbe")]

TOL_HF_ABS = 5.0e-2     # HF is grid-free -> matches the oracle closely
TOL_DELTA = 0.12        # coupling contribution, robust to cross-code DFT-SCF diff
TOL_DFT_ABS = 0.10      # DFT absolute para: cross-code DFT-SCF/grid offset
TOL_GATE_PURE = 1.0e-6  # pure DFT coupled == uncoupled


def _oqp_root():
    root = Path(os.environ.get("OPENQP_ROOT", str(ROOT)))
    lib = root / "lib" / "liboqp.so"
    if not lib.exists():
        lib = root / "lib" / "liboqp.dylib"
    if not lib.exists() or not (root / "include" / "oqp.h").exists():
        return None
    return str(root)


def _input(functional):
    fxc = f"functional={functional}\n" if functional else ""
    grid = "[dftgrid]\nrad_type=mhl\nrad_npts=99\nang_npts=302\n\n" if functional else ""
    return (f"[input]\nsystem=\n{GEOM}\ncharge=0\nruntype=energy\n{fxc}"
            f"basis=sto-3g\nmethod=hf\n\n[guess]\ntype=huckel\n\n"
            f"[scf]\nmultiplicity=1\ntype=rhf\n\n{grid}[properties]\nscf_prop=\n")


class GIAOParaShieldingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = _oqp_root()
        if cls.root is None:
            raise unittest.SkipTest("OpenQP shared library/header not built")
        if not FIXTURE.exists():
            raise unittest.SkipTest("PySCF GIAO reference fixture missing")
        cls.ref = json.loads(FIXTURE.read_text())["results"]
        cls.runs = {}
        for label, functional, _ in CASES:
            cls.runs[label] = cls._run(label, functional)

    @classmethod
    def _run(cls, label, functional):
        with tempfile.TemporaryDirectory() as wd:
            inp = Path(wd) / f"{label}.inp"
            log = Path(wd) / f"{label}.log"
            inp.write_text(_input(functional))
            script = Path(wd) / "run.py"
            script.write_text(textwrap.dedent(f"""
                import oqp
                from oqp.pyoqp import Runner
                r = Runner(project={label!r}, input_file={str(inp)!r}, log={str(log)!r}, silent=1, usempi=False)
                r.run()
                oqp.nmr_giao_shielding_debug(r.mol)
            """))
            env = dict(os.environ)
            env["OPENQP_ROOT"] = cls.root
            env["PYTHONPATH"] = os.pathsep.join(
                p for p in (str(ROOT / "pyoqp"), env.get("PYTHONPATH", "")) if p)
            proc = subprocess.run([sys.executable, str(script)], cwd=wd, env=env,
                                  capture_output=True, text=True, timeout=240)
            if not log.exists():
                raise unittest.SkipTest(f"{label}: GIAO run produced no log\n{proc.stderr[-1500:]}")
            return cls._parse(log.read_text())

    @staticmethod
    def _parse(text):
        nat_m = re.search(r"GIAO_SHIELDING_DEBUG_NATOM\s+(\d+)", text)
        if not nat_m:
            raise AssertionError("GIAO_SHIELDING_DEBUG markers not found")
        nat = int(nat_m.group(1))
        unc = [[[0.0] * 3 for _ in range(3)] for _ in range(nat)]
        cpl = [[[0.0] * 3 for _ in range(3)] for _ in range(nat)]
        for key, dest in (("PARA_UNC", unc), ("PARA_CPL", cpl)):
            pat = re.compile(
                rf"GIAO_SHIELDING_DEBUG_{key}\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.eE+\-]+)")
            for ia, t, s, v in pat.findall(text):
                dest[int(ia) - 1][int(t) - 1][int(s) - 1] = float(v)
        iso = lambda t: [sum(t[a][i][i] for i in range(3)) / 3.0 for a in range(nat)]
        return {"unc": iso(unc), "cpl": iso(cpl)}

    def test_hf_matches_oracle_absolutely(self):
        ref = self.ref["HF"]
        run = self.runs["HF"]
        for atom in range(3):
            self.assertAlmostEqual(run["unc"][atom], ref["sigma_para_uncoupled"][atom],
                                   delta=TOL_HF_ABS, msg=f"HF uncoupled para atom {atom+1}")
            self.assertAlmostEqual(run["cpl"][atom], ref["sigma_para_coupled"][atom],
                                   delta=TOL_HF_ABS, msg=f"HF coupled para atom {atom+1}")

    def test_coupling_delta_matches_oracle(self):
        for label, _, key in CASES:
            ref = self.ref[key]
            run = self.runs[label]
            for atom in range(3):
                d_oqp = run["cpl"][atom] - run["unc"][atom]
                d_ref = (ref["sigma_para_coupled"][atom]
                         - ref["sigma_para_uncoupled"][atom])
                self.assertAlmostEqual(d_oqp, d_ref, delta=TOL_DELTA,
                                       msg=f"{label} coupling delta atom {atom+1}")

    def test_para_absolute_within_crosscode_tolerance(self):
        for label, _, key in CASES:
            ref = self.ref[key]
            run = self.runs[label]
            tol = TOL_HF_ABS if label == "HF" else TOL_DFT_ABS
            for atom in range(3):
                self.assertAlmostEqual(run["cpl"][atom], ref["sigma_para_coupled"][atom],
                                       delta=tol, msg=f"{label} coupled para atom {atom+1}")

    def test_pure_dft_coupled_equals_uncoupled(self):
        run = self.runs["PBE"]
        for atom in range(3):
            self.assertAlmostEqual(run["unc"][atom], run["cpl"][atom],
                                   delta=TOL_GATE_PURE,
                                   msg="PBE (pure) GIAO coupled must equal uncoupled")


if __name__ == "__main__":
    unittest.main()
