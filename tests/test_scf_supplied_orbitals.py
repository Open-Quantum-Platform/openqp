"""Second-order SCF starts directly from valid supplied orbitals."""

import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source"
KEPT = "Second-order converger starts from the supplied orbitals."

DECK = """[input]
system=
   8   0.000000000   0.000000000   0.000000000
   1   0.000000000  -0.757000000   0.587000000
   1   0.000000000   0.757000000   0.587000000
charge=0
runtype=energy
basis=6-31g
method=hf

[guess]
{guess}

[scf]
type={scf}
multiplicity={mult}
converger_type={converger}
trh_impl=native
conv={conv}
"""


def _runtime_available():
    try:
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


class TestGuessOriginSource(unittest.TestCase):
    def test_model_guesses_are_cold_and_supplied_guesses_are_recorded(self):
        for name in ("guess_huckel", "guess_hcore", "guess_minao", "guess_sap"):
            src = (SOURCE / "modules" / f"{name}.F90").read_text()
            self.assertIn("infos%control%guess = GUESS_COLD", src, name)
        self.assertIn(
            "infos%control%guess = GUESS_SUPPLIED",
            (SOURCE / "modules" / "guess_json.F90").read_text(),
        )

        guess_py = (ROOT / "pyoqp" / "oqp" / "library" / "guess.py").read_text()
        self.assertIn("supplied = alpha in ('reloaded', 'reused')", guess_py)
        self.assertIn(
            "control.guess = GUESS_SUPPLIED if supplied else GUESS_COLD",
            guess_py,
        )

    def test_kept_orbitals_use_current_fock_energies(self):
        src = (SOURCE / "scf.F90").read_text()
        self.assertRegex(
            src,
            r"keep_supplied = orthonormal_orbitals\(mo_a, smat_full, nbf\)",
        )
        self.assertRegex(
            src,
            r"if \(keep_supplied\) then[\s\S]{0,700}?"
            r"call fock_diagonal_energies\(pfock\(:,1\), mo_a, mo_energy_a, nbf\)",
        )
        self.assertIn(
            "if ((use_soscf .or. use_trah) .and. iter == 1 .and. .not. do_pfon .and. &",
            src,
        )
        self.assertRegex(
            src,
            r"if \(infos%mol_energy%SCF_converged\) then\s*\n"
            r"\s*infos%control%guess = GUESS_SUPPLIED\s*\n"
            r"\s*else\s*\n\s*infos%control%guess = GUESS_COLD",
        )


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestSuppliedOrbitalsRuntime(unittest.TestCase):
    def _run(self, tmp, name, guess, scf, mult, converger, conv):
        deck = Path(tmp) / f"{name}.inp"
        deck.write_text(
            DECK.format(
                guess=guess,
                scf=scf,
                mult=mult,
                converger=converger,
                conv=conv,
            )
        )
        env = dict(os.environ, OMP_NUM_THREADS="2")
        try:
            proc = subprocess.run(
                [sys.executable, "-m", "oqp.pyoqp", deck.name],
                cwd=tmp,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=240,
            )
        except subprocess.TimeoutExpired:
            self.fail(f"{name} did not finish within 240 s")

        output = proc.stdout.decode(errors="ignore")
        log = (Path(tmp) / f"{name}.log").read_text(errors="ignore")
        self.assertEqual(proc.returncode, 0, output[-2000:] + log[-3000:])
        energies = re.findall(
            r"Final \w+ energy is\s+(-?\d+\.\d+) after\s+\d+ iterations", log
        )
        self.assertTrue(energies, f"{name}: no final SCF energy in the log")
        self.assertNotIn("escalating to", log)
        self.assertIn("SCF convergence achieved", log)
        return log, float(energies[-1])

    def _check(self, scf, mult):
        with tempfile.TemporaryDirectory() as tmp:
            first_log, _ = self._run(
                tmp,
                "first",
                "type=huckel\nsave_mol=true",
                scf,
                mult,
                "soscf",
                "1e-3",
            )
            self.assertNotIn(KEPT, first_log)
            self.assertTrue((Path(tmp) / "first.json").exists())
            reload_guess = "type=json\nfile=first.json"

            soscf_log, e_soscf = self._run(
                tmp, "soscf", reload_guess, scf, mult, "soscf", "1e-8"
            )
            self.assertIn(KEPT, soscf_log)

            trah_log, e_trah = self._run(
                tmp, "trah", reload_guess, scf, mult, "trah", "1e-8"
            )
            self.assertIn(KEPT, trah_log)
            residuals = re.findall(
                r"final fresh-Fock residual\s*=\s*([0-9.E+-]+)", trah_log
            )
            self.assertTrue(residuals, "TRAH did not report its final fresh-Fock residual")
            self.assertLess(float(residuals[-1]), 1.0e-8)

            diis_log, e_diis = self._run(
                tmp, "diis", reload_guess, scf, mult, "diis", "1e-8"
            )
            self.assertNotIn(KEPT, diis_log)
            self.assertAlmostEqual(e_soscf, e_diis, delta=1e-8)
            if scf == "uhf":
                # UHF can have several stationary solutions. TRAH may leave the
                # supplied basin and converge to a lower one, but it must not
                # finish above the DIIS solution used as the reference here.
                self.assertLessEqual(e_trah, e_diis + 1e-8)
            else:
                self.assertAlmostEqual(e_trah, e_diis, delta=1e-8)

    def test_rhf_json_reload(self):
        self._check("rhf", 1)

    def test_rohf_triplet_json_reload(self):
        self._check("rohf", 3)

    def test_uhf_triplet_json_reload(self):
        self._check("uhf", 3)


if __name__ == "__main__":
    unittest.main()
