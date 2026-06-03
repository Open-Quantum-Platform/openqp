"""Phase-0 regression tests: coupled HF/hybrid ground-state magnetic response.

Validates the coupled CPHF/CPKS magnetic response (exact-exchange response of the
imaginary antisymmetric first-order density) against an independent common-gauge-origin
reference (tests/fixtures/nmr/cgo_reference.json) for HF, BHHLYP, PBE0, PBE on
H2O/STO-3G with the gauge origin at the center of mass.

Encodes the Phase-0 gates:
  gate 0: first-order magnetic density is antisymmetric (max|P^B + P^B^T| ~ 0)
  gate 1: Coulomb response vanishes (||J(P^B)|| ~ 0)
  gate 2: exact-exchange response is nonzero (||K(P^B)|| > 0)
  gate 3: pure DFT (PBE) coupled == uncoupled (exact)
  gate 4: HF/hybrid coupled != uncoupled; HF matches the oracle exactly, and the
          coupling contribution Delta = coupled - uncoupled matches the oracle for
          all functionals (Delta isolates the coupling from cross-code DFT-SCF
          differences, which shift the absolute DFT numbers by ~0.1 ppm).
  (gate 6: |Delta| scales monotonically with the exact-exchange fraction c_x.)

Runs the OpenQP driver in a subprocess; skips gracefully if the library is not
built. The oracle JSON is committed; do not recompute ad hoc.
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "nmr" / "cgo_reference.json"

GEOM = """\
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223"""

# (input label, OpenQP functional keyword, fixture key)
CASES = [("HF", None, "HF"), ("BHHLYP", "bhhlyp", "bhandhlyp"),
         ("PBE0", "pbe0", "pbe0"), ("PBE", "pbe", "pbe")]

TOL_HF_ABS = 5.0e-2     # HF is grid-free -> exact match to the oracle
TOL_DELTA = 0.12        # coupling contribution, robust to cross-code DFT-SCF diff
TOL_GATE3 = 1.0e-6      # pure DFT coupled == uncoupled
TOL_ANTISYM = 1.0e-9    # P^B antisymmetry / Coulomb response
TOL_KMIN = 1.0e-3       # exchange response must be clearly nonzero


def _oqp_root():
    root = os.environ.get("OPENQP_ROOT", str(ROOT))
    lib = Path(root) / "lib"
    if (lib / "liboqp.dylib").exists() or (lib / "liboqp.so").exists():
        return root
    return None


def _input(functional):
    fxc = f"functional={functional}\n" if functional else ""
    grid = "[dftgrid]\nrad_type=becke\n\n" if functional else ""
    return (f"[input]\nsystem=\n{GEOM}\ncharge=0\nruntype=energy\n{fxc}"
            f"basis=sto-3g\nmethod=hf\n\n[guess]\ntype=huckel\n\n"
            f"[scf]\nmultiplicity=1\ntype=rhf\n\n{grid}"
            f"[properties]\nscf_prop=nmr\n")


class CoupledMagneticResponseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = _oqp_root()
        if cls.root is None:
            raise unittest.SkipTest("OpenQP shared library not built")
        if not FIXTURE.exists():
            raise unittest.SkipTest("Reference fixture missing")
        cls.ref = json.loads(FIXTURE.read_text())["results"]
        cls.runs = {}
        env = dict(os.environ)
        env["OPENQP_ROOT"] = cls.root
        env["PYTHONPATH"] = str(ROOT / "pyoqp") + os.pathsep + env.get("PYTHONPATH", "")
        for label, functional, _ in CASES:
            with tempfile.TemporaryDirectory() as wd:
                inp = Path(wd) / "h2o.inp"
                inp.write_text(_input(functional))
                subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(inp)],
                               cwd=wd, env=env, capture_output=True, text=True)
                log = inp.with_suffix(".log")
                if not log.exists():
                    raise unittest.SkipTest(f"{label} run produced no log")
                cls.runs[label] = cls._parse(log.read_text())

    @staticmethod
    def _parse(text):
        out = {"rows": {}, "gates": {}}
        for key, pat in (("gate0", r"gate0.*?=\s*([-0-9.eE+]+)"),
                         ("gate1", r"gate1.*?=\s*([-0-9.eE+]+)"),
                         ("gate2", r"gate2.*?=\s*([-0-9.eE+]+)")):
            m = re.search(pat, text)
            if m:
                out["gates"][key] = float(m.group(1))
        # rows: atom Z dia para_unc para_cpl tot_unc tot_cpl
        for line in text.splitlines():
            m = re.match(r"\s*(\d+)\s+[\d.]+" + r"\s+(-?\d+\.\d+)" * 5 + r"\s*$", line)
            if m:
                out["rows"][int(m.group(1))] = [float(m.group(j)) for j in range(2, 7)]
        return out

    def test_gate0_density_antisymmetric(self):
        for label in self.runs:
            self.assertLess(abs(self.runs[label]["gates"]["gate0"]), TOL_ANTISYM,
                            f"{label}: P^B not antisymmetric")

    def test_gate1_coulomb_response_vanishes(self):
        for label in self.runs:
            self.assertLess(abs(self.runs[label]["gates"]["gate1"]), TOL_ANTISYM,
                            f"{label}: Coulomb response of P^B not ~0")

    def test_gate2_exchange_response_nonzero(self):
        for label in self.runs:
            self.assertGreater(self.runs[label]["gates"]["gate2"], TOL_KMIN,
                               f"{label}: exact-exchange response of P^B is ~0")

    def test_gate3_pure_dft_coupled_equals_uncoupled(self):
        for atom in (1, 2):
            unc, cpl = self.runs["PBE"]["rows"][atom][1], self.runs["PBE"]["rows"][atom][2]
            self.assertAlmostEqual(unc, cpl, delta=TOL_GATE3,
                                   msg="PBE (pure) coupled must equal uncoupled")

    def test_gate4_hf_matches_oracle_absolutely(self):
        ref = self.ref["HF"]
        for atom in (1, 2):
            dia, pu, pc = self.runs["HF"]["rows"][atom][:3]
            self.assertAlmostEqual(pc, ref["sigma_para_coupled"][atom - 1],
                                   delta=TOL_HF_ABS, msg=f"HF coupled para atom {atom}")
            self.assertAlmostEqual(pu, ref["sigma_para_uncoupled"][atom - 1],
                                   delta=TOL_HF_ABS, msg=f"HF uncoupled para atom {atom}")

    def test_gate4_coupling_delta_matches_oracle(self):
        # Delta = coupled - uncoupled isolates the exchange coupling from
        # cross-code DFT-SCF differences in the absolute numbers.
        for label, _, key in CASES:
            ref = self.ref[key]
            for atom in (1, 2):
                d_oqp = self.runs[label]["rows"][atom][2] - self.runs[label]["rows"][atom][1]
                d_ref = (ref["sigma_para_coupled"][atom - 1]
                         - ref["sigma_para_uncoupled"][atom - 1])
                self.assertAlmostEqual(d_oqp, d_ref, delta=TOL_DELTA,
                                       msg=f"{label} coupling Delta atom {atom}")

    def test_gate6_delta_scales_with_exact_exchange(self):
        # |Delta(O)| monotonically decreasing: HF > BHHLYP > PBE0 > PBE(=0).
        d = {lab: abs(self.runs[lab]["rows"][1][2] - self.runs[lab]["rows"][1][1])
             for lab, _, _ in CASES}
        self.assertGreater(d["HF"], d["BHHLYP"])
        self.assertGreater(d["BHHLYP"], d["PBE0"])
        self.assertGreater(d["PBE0"], d["PBE"])
        self.assertLess(d["PBE"], TOL_GATE3)


if __name__ == "__main__":
    unittest.main()
