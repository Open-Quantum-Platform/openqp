"""STG1X = 0.85 HF + 0.15 B88 exchange + LYP correlation (JPCA 123, 7991).

The SCF energies are pinned to an independent code: ORCA 5.0.3 with the same
hybrid built explicitly (BLYP, %method ScalHFX 0.85, ScalDFX 0.15), cc-pVDZ
(spherical), VeryTightSCF, NoRI, DefGrid3.  Agreement is grid-level (< 1e-6
Eh).  The former definition (BHANDHLYP's 0.5 B88 + LYP under 0.85 HF, total
exchange 1.35) is ~3 Eh lower for H2O and fails by a wide margin.
"""
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ORCA_STG1X = {  # Eh, ORCA 5.0.3 (chc r630), see module docstring
    "h2o_rks": -76.371572666922,
    "ch2_uks": -39.136625262077,
}
SYSTEMS = {
    "h2o_rks": ("   O 0.000000 0.000000 0.117300\n   H 0.000000 0.757200 -0.469200\n"
                "   H 0.000000 -0.757200 -0.469200\n", 1, "rhf"),
    "ch2_uks": ("   C 0.000000 0.000000 0.102700\n   H 0.000000 0.995600 -0.308000\n"
                "   H 0.000000 -0.995600 -0.308000\n", 3, "uhf"),
}


def _runtime_available():
    try:
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


def _scf_energy(key):
    geom, mult, scf_type = SYSTEMS[key]
    text = (f"[input]\nsystem=\n{geom}charge=0\nruntype=energy\nbasis=cc-pvdz\nispher=true\n"
            "method=hf\nfunctional=stg1x\n\n[guess]\ntype=huckel\n\n"
            f"[scf]\ntype={scf_type}\nmultiplicity={mult}\nmaxit=200\nconv=1e-10\n\n"
            "[dftgrid]\nrad_npts=150\nang_npts=974\npruned=none\n")
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "job.inp").write_text(text)
        env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "2"))
        proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", "job.inp"], cwd=tmp, env=env,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=900)
        log_path = Path(tmp) / "job.log"
        log = log_path.read_text(errors="ignore") if log_path.exists() else ""
        found = re.findall(r"Final \S+ energy is\s+(-?[\d.]+)", log)
        return proc.returncode, (float(found[-1]) if found else None), log[-2000:]


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestStg1xDefinition(unittest.TestCase):
    def test_scf_energy_matches_orca_explicit_hybrid(self):
        for key, ref in ORCA_STG1X.items():
            with self.subTest(system=key):
                rc, energy, tail = _scf_energy(key)
                self.assertEqual(rc, 0, tail)
                self.assertIsNotNone(energy, tail)
                self.assertAlmostEqual(energy, ref, delta=5e-6)


if __name__ == "__main__":
    unittest.main()
