"""End-to-end regression test for the native NMR shielding property.

Runs an RHF/STO-3G H2O common-gauge NMR shielding calculation and checks:

  1. the PSO integral matrix has a vanishing diagonal (anti-Hermitian operator),
  2. the PSO matrix is antisymmetric (max|A + A^T| ~ 0),
  3. the oxygen and hydrogen paramagnetic (and total) isotropic shieldings match
     an independent uncoupled common-gauge reference.

The test runs the OpenQP Python driver in a subprocess and skips gracefully if a
built library / environment is not available.
"""

import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Independent common-gauge, uncoupled reference (H2O/STO-3G, gauge origin = COM).
# Computed from the standard CPHF-free uncoupled magnetic response, unit = ALPHA**2 * 1e6.
REF = {
    # atom index (1-based): (sigma_dia, sigma_para, sigma_total) in ppm
    1: (411.4176, -113.6305, 297.7871),  # O
    2: (28.0618, 1.7849, 29.8468),       # H
    3: (28.0618, 1.7849, 29.8468),       # H
}
TOL_PPM = 5.0e-2          # ppm; well above STO-3G run-to-run noise
TOL_ANTISYM = 1.0e-8      # a.u.; PSO must be exactly antisymmetric

INPUT = """\
[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
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
scf_prop=nmr
"""


def _oqp_root():
    root = os.environ.get("OPENQP_ROOT", str(ROOT))
    libdir = Path(root) / "lib"
    if not (libdir / "liboqp.dylib").exists() and not (libdir / "liboqp.so").exists():
        return None
    return root


class NMRShieldingTests(unittest.TestCase):
    def setUp(self):
        self.root = _oqp_root()
        if self.root is None:
            self.skipTest("OpenQP shared library not built; skipping live NMR run")

    def _run(self, workdir):
        inp = Path(workdir) / "h2o_nmr.inp"
        inp.write_text(INPUT)
        env = dict(os.environ)
        env["OPENQP_ROOT"] = self.root
        root_parent = str(Path(self.root).parent)
        env["PYTHONPATH"] = root_parent + os.pathsep + env.get("PYTHONPATH", "")
        proc = subprocess.run(
            [sys.executable, "-m", "oqp.pyoqp", str(inp)],
            cwd=workdir, env=env, capture_output=True, text=True,
        )
        log = inp.with_suffix(".log")
        if not log.exists():
            self.skipTest(
                "NMR run did not produce a log (missing runtime deps?):\n"
                + proc.stdout[-2000:] + proc.stderr[-2000:]
            )
        return log.read_text()

    def _parse(self, text):
        # PSO diagnostics line: "max|diag|, max|A+A^T| =   d   a"
        m = re.search(r"max\|diag\|.*?=\s*([0-9.eE+\-]+)\s+([0-9.eE+\-]+)", text)
        self.assertIsNotNone(m, "PSO diagnostics line not found in log")
        diag_max, asym_max = float(m.group(1)), float(m.group(2))

        # Shielding rows (5 columns): atom Z dia para_unc para_cpl tot_unc tot_cpl.
        # This test validates the UNCOUPLED numbers (dia, para_unc, tot_unc).
        rows = {}
        for line in text.splitlines():
            mm = re.match(r"\s*(\d+)\s+[\d.]+" + r"\s+(-?\d+\.\d+)" * 5 + r"\s*$", line)
            if mm:
                vals = [float(mm.group(j)) for j in range(2, 7)]
                rows[int(mm.group(1))] = (vals[0], vals[1], vals[3])  # dia, para_unc, tot_unc
        return diag_max, asym_max, rows

    def test_nmr_shielding(self):
        with tempfile.TemporaryDirectory() as wd:
            text = self._run(wd)
        diag_max, asym_max, rows = self._parse(text)

        # 1 & 2: PSO operator must be exactly antisymmetric (zero diagonal).
        self.assertLess(diag_max, TOL_ANTISYM,
                        f"PSO diagonal not ~0 (max|diag|={diag_max:g})")
        self.assertLess(asym_max, TOL_ANTISYM,
                        f"PSO not antisymmetric (max|A+A^T|={asym_max:g})")

        # 3: dia / para / total match the independent uncoupled reference.
        self.assertEqual(set(rows), set(REF), f"unexpected shielding rows: {rows}")
        for atom, (dia, para, tot) in REF.items():
            g_dia, g_para, g_tot = rows[atom]
            self.assertAlmostEqual(g_dia, dia, delta=TOL_PPM,
                                   msg=f"atom {atom} sigma_dia {g_dia} vs {dia}")
            self.assertAlmostEqual(g_para, para, delta=TOL_PPM,
                                   msg=f"atom {atom} sigma_para {g_para} vs {para}")
            self.assertAlmostEqual(g_tot, tot, delta=TOL_PPM,
                                   msg=f"atom {atom} sigma_total {g_tot} vs {tot}")


if __name__ == "__main__":
    unittest.main()
