"""End-to-end PCM energy-path tests.

These drive the built OpenQP library through a real H2O RHF/6-31G* single point
and prove the first closed PCM reaction-field loop:

    D -> phi_cav -> ddX q -> V_pcm -> Fock/E_pcm -> D'

Two guarantees are checked:

* PCM-off (mandatory): the vacuum SCF path is preserved exactly. A run with no
  ``[pcm]`` section and a run with ``pcm.enabled=false`` both reproduce the
  established vacuum reference energy and are identical to each other, and no
  PCM energy term is emitted.

* PCM-on (skip-guarded): the reaction-field loop runs to completion with a
  finite, nonzero solvent contribution. This requires a ddX-enabled build
  (OQP_ENABLE_DDX); when ddX is unavailable the C adapter reports so and the
  test is skipped rather than failed.

The whole module is skipped when ``lib/liboqp`` has not been built.

NOTE: PCM-on only asserts that the loop closes and produces a finite nonzero
contribution. The q sign/scale and the phi_cav convention are PROVISIONAL and
are deliberately NOT checked against a trusted PCM reference at this gate.
"""

import os
import platform
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _lib_suffix() -> str:
    return {"Windows": "dll", "Linux": "so", "Darwin": "dylib"}.get(
        platform.uname()[0], "so"
    )


LIBOQP = ROOT / "lib" / f"liboqp.{_lib_suffix()}"

# Established vacuum RHF/6-31G* energy for this H2O geometry
# (examples/HF/H2O_RHF-HF_ENERGY.json).
VACUUM_RHF_6_31GS = -76.01074651512697

H2O_SYSTEM = (
    "   8   0.000000000   0.000000000  -0.041061554\n"
    "   1  -0.533194329   0.533194329  -0.614469223\n"
    "   1   0.533194329  -0.533194329  -0.614469223"
)


def _input_text(pcm_section: str = "") -> str:
    return (
        "[input]\n"
        "system=\n"
        f"{H2O_SYSTEM}\n"
        "charge=0\n"
        "runtype=energy\n"
        "basis=6-31g*\n"
        "method=hf\n"
        "\n"
        "[guess]\n"
        "type=huckel\n"
        "\n"
        "[scf]\n"
        "multiplicity=1\n"
        "type=rhf\n"
        f"{pcm_section}"
    )


PCM_ON_SECTION = (
    "\n[pcm]\n"
    "enabled=true\n"
    "backend=ddx\n"
    "mode=reference_scf\n"
    "model=ddpcm\n"
    "epsilon=78.3553\n"
)


@unittest.skipUnless(
    LIBOQP.is_file(), f"liboqp not built at {LIBOQP}; build OpenQP first"
)
class TestPcmEnergyPath(unittest.TestCase):
    def _run(self, text: str):
        with tempfile.TemporaryDirectory() as d:
            inp = Path(d) / "h2o.inp"
            inp.write_text(text)
            env = dict(os.environ)
            env["OPENQP_ROOT"] = str(ROOT)
            env["PYTHONPATH"] = (
                str(ROOT / "pyoqp") + os.pathsep + env.get("PYTHONPATH", "")
            )
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

    @staticmethod
    def _total_energy(log: str):
        m = re.findall(r"TOTAL energy\s*=\s*(-?\d+\.\d+)", log)
        return float(m[-1]) if m else None

    @staticmethod
    def _pcm_energy(log: str):
        m = re.findall(r"PCM solvent energy.*?=\s*(-?\d+\.\d+)", log)
        return float(m[-1]) if m else None

    def test_pcm_off_matches_vacuum_reference(self):
        # (a) no [pcm] section -> pure vacuum path.
        proc0, log0 = self._run(_input_text())
        self.assertEqual(proc0.returncode, 0, log0)
        e0 = self._total_energy(log0)
        self.assertIsNotNone(e0, log0)
        self.assertAlmostEqual(e0, VACUUM_RHF_6_31GS, places=7)

        # (b) explicit pcm.enabled=false must reproduce (a) exactly and emit no
        #     PCM energy term: the off-switch leaves the vacuum path untouched.
        proc1, log1 = self._run(_input_text("\n[pcm]\nenabled=false\n"))
        self.assertEqual(proc1.returncode, 0, log1)
        e1 = self._total_energy(log1)
        self.assertIsNotNone(e1, log1)
        self.assertEqual(e1, e0)
        self.assertIsNone(self._pcm_energy(log1))

    def test_pcm_on_closes_loop(self):
        proc, log = self._run(_input_text(PCM_ON_SECTION))

        if "OQP_ENABLE_DDX" in log or "built without" in log:
            self.skipTest(
                "OpenQP built without ddX (OQP_ENABLE_DDX); the PCM-on "
                "reaction-field path is wired but not executable here."
            )

        self.assertEqual(proc.returncode, 0, log)

        e_pcm = self._pcm_energy(log)
        self.assertIsNotNone(e_pcm, "PCM energy term was not reported:\n" + log)
        self.assertNotEqual(e_pcm, 0.0, "PCM energy term is zero; loop not closed.")

        e_tot = self._total_energy(log)
        self.assertIsNotNone(e_tot, log)
        # The solvent contribution must move the total off the vacuum value.
        self.assertNotAlmostEqual(e_tot, VACUUM_RHF_6_31GS, places=7)


if __name__ == "__main__":
    unittest.main()
