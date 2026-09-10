"""[qmmm] forcefield_files must build the PDB-based QM molecule.

The active QM/MM drivers read only [qmmm] forcefield_files, but the QM
molecule of a deck with ``[input] system = file.pdb <1-based indices>`` was
built from the legacy [qmmm] forcefield (default AMBER-14), so a residue
defined by a custom XML failed before any driver ran:
``No template found for residue 1 (FOR)`` for formaldehyde in water.
"""
import os
import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples" / "QMMM"


class TestResolveForceFieldFiles(unittest.TestCase):
    def setUp(self):
        from oqp.utils.qmmm import resolve_forcefield_files
        self.resolve = resolve_forcefield_files

    def test_lists_and_builtins(self):
        self.assertEqual(self.resolve("amber14-all.xml amber14/tip3p.xml"), ["amber14-all.xml", "amber14/tip3p.xml"])
        self.assertEqual(self.resolve("amber14-all.xml,amber14/tip3p.xml"), ["amber14-all.xml", "amber14/tip3p.xml"])
        self.assertEqual(self.resolve(""), [])
        self.assertEqual(self.resolve(None), [])

    def test_deck_relative_names_and_spaces(self):
        with tempfile.TemporaryDirectory() as deck, tempfile.TemporaryDirectory() as cwd:
            (Path(deck) / "my ff.xml").write_text("<ForceField/>")
            (Path(deck) / "custom.xml").write_text("<ForceField/>")
            here = os.getcwd(); os.chdir(cwd)
            try:
                self.assertEqual(self.resolve("my ff.xml", deck), [str(Path(deck) / "my ff.xml")])
                self.assertEqual(self.resolve("custom.xml amber14/tip3p.xml", deck),
                                 [str(Path(deck) / "custom.xml"), "amber14/tip3p.xml"])
            finally:
                os.chdir(here)

    def test_the_handler_is_registered_after_forcefield(self):
        src = (ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py").read_text()
        i_ff = src.index('"forcefield": "set_qmmm_forcefield",')
        i_ffs = src.index('"forcefield_files": "set_qmmm_forcefield_files",')
        self.assertLess(i_ff, i_ffs)
        # [qmmm] is applied before [input] (whose system key builds the molecule)
        self.assertLess(src.index("    'qmmm': {"), src.index("    'input': {"))


def _runtime():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        import oqp  # noqa: F401
        from oqp.molecule import Molecule  # noqa: F401
        import openmm  # noqa: F401
        return (EXAMPLES / "formaldehyde.xml").exists()
    except Exception:
        return False


DECK = """[input]
system=formaldehyde_water.pdb 1 2 3 4
charge=0
runtype=energy
basis=sto-3g
method=hf
qmmm_flag=True
[scf]
type=rhf
multiplicity=1
[qmmm]
pdb_file=formaldehyde_water.pdb
forcefield_files=formaldehyde.xml tip3p.xml
qm_atoms=0-3
cutoff=NoCutoff
"""


@unittest.skipUnless(_runtime(), "compiled OpenQP runtime or OpenMM unavailable")
class TestPdbMoleculeWithCustomResidue(unittest.TestCase):
    def test_custom_residue_from_forcefield_files_next_to_the_deck(self):
        from oqp.molecule import Molecule
        from oqp.utils import qmmm
        with tempfile.TemporaryDirectory() as deck_dir, tempfile.TemporaryDirectory() as cwd:
            for name in ("formaldehyde_water.pdb", "formaldehyde.xml", "tip3p.xml"):
                shutil.copy(EXAMPLES / name, Path(deck_dir) / name)
            deck = Path(deck_dir) / "m.inp"
            deck.write_text(DECK)
            here = os.getcwd(); os.chdir(cwd)          # files resolvable only next to the deck
            try:
                mol = Molecule("m", str(deck), str(Path(cwd) / "m.log"), silent=1)
                mol.load_config(str(deck))
            finally:
                os.chdir(here)
        self.assertEqual(int(mol.data["natom"]), 4)
        self.assertEqual([os.path.basename(f) for f in qmmm.force_field], ["formaldehyde.xml", "tip3p.xml"])

    def test_without_forcefield_files_the_legacy_default_is_kept(self):
        from oqp.molecule import Molecule
        from oqp.utils import qmmm
        deck_text = DECK.replace("forcefield_files=formaldehyde.xml tip3p.xml\n", "")
        with tempfile.TemporaryDirectory() as deck_dir:
            shutil.copy(EXAMPLES / "formaldehyde_water.pdb", Path(deck_dir) / "formaldehyde_water.pdb")
            deck = Path(deck_dir) / "m.inp"
            deck.write_text(deck_text)
            mol = Molecule("m", str(deck), str(Path(deck_dir) / "m.log"), silent=1)
            with self.assertRaisesRegex(Exception, "No template found"):
                mol.load_config(str(deck))                # the control: AMBER alone cannot build FOR
        self.assertEqual(list(qmmm.force_field), ["amber14-all.xml", "amber14/tip3p.xml"])


if __name__ == "__main__":
    unittest.main()
