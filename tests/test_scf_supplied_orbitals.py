"""SOSCF/TRAH start from supplied orbitals without re-diagonalising the first Fock.

The shared converger takes a steepest-descent step on its first iteration, which
refilled supplied orbitals in the order of the (ROHF effective) orbital energies
and could swap the occupations of a converged solution.  Orbitals reloaded from a
JSON guess, or left by a converged SCF, now start a second-order converger as they
are; model guesses and DIIS keep the diagonalisation.
"""
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
conv={conv}
"""


def _runtime_available():
    try:
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


class TestGuessMarkerSource(unittest.TestCase):
    def test_model_guesses_mark_cold_and_reloads_mark_supplied(self):
        for name in ("guess_huckel", "guess_hcore", "guess_minao", "guess_sap"):
            src = (SOURCE / "modules" / f"{name}.F90").read_text()
            self.assertIn("infos%control%guess = GUESS_COLD", src, name)
        self.assertIn("infos%control%guess = GUESS_SUPPLIED",
                      (SOURCE / "modules" / "guess_json.F90").read_text())
        guess_py = (ROOT / "pyoqp" / "oqp" / "library" / "guess.py").read_text()
        self.assertIn("mol.data._data.control.guess = GUESS_SUPPLIED if alpha == 'reloaded' else GUESS_COLD",
                      guess_py)
        types = (SOURCE / "types.F90").read_text()
        self.assertIn("GUESS_COLD = 1, GUESS_SUPPLIED = 2", types)
        self.assertIn("GUESS_COLD, GUESS_SUPPLIED = 1, 2", guess_py)

    def test_kept_orbitals_take_energies_from_the_current_fock(self):
        src = (SOURCE / "scf.F90").read_text()
        self.assertRegex(src, r"keep_supplied = orthonormal_orbitals\(mo_a, smat_full, nbf\)")
        self.assertRegex(src, r"if \(keep_supplied\) then\s*\n(?:\s*!.*\n)*"
                              r"\s*call fock_diagonal_energies\(pfock\(:,1\), mo_a, mo_energy_a, nbf\)\s*\n"
                              r"\s*if \(scf_type == scf_uhf \.and\. nelec_b /= 0\) &\s*\n"
                              r"\s*call fock_diagonal_energies\(pfock\(:,2\), mo_b, mo_energy_b, nbf\)")
        # pFON places its Fermi level by orbital index, so it keeps the diagonalisation
        self.assertIn("if ((use_soscf .or. use_trah) .and. iter == 1 .and. .not. do_pfon .and. &", src)
        # a converged SCF hands its orbitals on; an unconverged one does not
        self.assertRegex(src, r"if \(infos%mol_energy%SCF_converged\) then\s*\n\s*infos%control%guess = GUESS_SUPPLIED\s*\n"
                              r"\s*else\s*\n\s*infos%control%guess = GUESS_COLD")


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestSuppliedOrbitalsRuntime(unittest.TestCase):
    def _run(self, tmp, name, guess, scf, mult, converger, conv):
        deck = Path(tmp) / f"{name}.inp"
        deck.write_text(DECK.format(guess=guess, scf=scf, mult=mult, converger=converger, conv=conv))
        env = dict(os.environ, OMP_NUM_THREADS="2")
        try:
            proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", deck.name], cwd=tmp, env=env,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=240)
        except subprocess.TimeoutExpired:
            self.fail(f"{name} did not finish within 240 s")
        self.assertEqual(proc.returncode, 0, proc.stdout.decode(errors="ignore")[-2000:])
        log = (Path(tmp) / f"{name}.log").read_text(errors="ignore")
        energies = re.findall(r"Final \w+ energy is\s+(-?\d+\.\d+) after\s+\d+ iterations", log)
        self.assertTrue(energies, f"{name}: no final SCF energy in the log")
        # the requested converger must finish on its own: an escalation would
        # report the energy of another solver's attempt
        self.assertNotIn("escalating to", log)
        self.assertIn("SCF convergence achieved", log)
        return log, float(energies[-1])

    def _check(self, scf, mult):
        with tempfile.TemporaryDirectory() as tmp:
            # loosely converged orbitals, so a reloaded SCF still iterates past
            # its first convergence test and reaches the orbital update
            log, _ = self._run(tmp, "first", "type=huckel\nsave_mol=true", scf, mult, "soscf", "1e-3")
            self.assertNotIn(KEPT, log)          # a model guess is diagonalised
            self.assertTrue((Path(tmp) / "first.json").exists())
            reload = "type=json\nfile=first.json"
            log, e_kept = self._run(tmp, "soscf", reload, scf, mult, "soscf", "1e-8")
            self.assertIn(KEPT, log)
            log, e_diis = self._run(tmp, "diis", reload, scf, mult, "diis", "1e-8")
            self.assertNotIn(KEPT, log)          # DIIS keeps the diagonalisation
            self.assertAlmostEqual(e_kept, e_diis, delta=1e-8)

    def test_rhf_json_reload(self):
        # RHF reloads bypass guess_json; the Python guess marks them supplied
        self._check("rhf", 1)

    def test_rohf_triplet_json_reload(self):
        self._check("rohf", 3)

    def test_uhf_triplet_json_reload(self):
        self._check("uhf", 3)


if __name__ == "__main__":
    unittest.main()
