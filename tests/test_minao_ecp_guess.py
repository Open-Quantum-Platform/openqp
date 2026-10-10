"""MINAO initial guess for elements beyond Kr and for atoms with an ECP.

The MINAO table stopped at Z = 36, so IF/def2-SVP aborted with "element beyond
tabulated range" before any SCF.  The table now covers H-Xe, the range of
OpenQP's STO-3G basis.  An atom with an effective core potential also needs its
core orbitals removed from the tabulated all-electron density: the def2 iodine
ECP replaces 28 core electrons, so the atomic density must hold 25 electrons,
not 53.  The core orbitals are the doubly occupied atomic natural orbitals that
lie within the span of the 1s-3d minimal-basis functions.

Checked here, with integral symmetry disabled:
  * IF/def2-SVP (I beyond the old table, def2 ECP with 28 core electrons):
    the guess reports 25 iodine electrons and RHF reaches the Hueckel energy;
  * HBr/LANL2DZ (tabulated Br, ECP with 28 core electrons): 7 Br electrons,
    same energy as Hueckel;
  * all-electron IF/3-21G (iodine beyond the old table, no ECP): 53 electrons,
    same energy as Hueckel;
  * the data file covers Z = 1-54, each with Z electrons and the Cartesian
    STO-3G size of that element.

The runtime cases are skipped unless the compiled OpenQP runtime is importable.
"""

import os
import re
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "basis_sets" / "minao_sto3g.dat"
STO3G = ROOT / "basis_sets" / "sto-3g.basis"

INPUT_TMPL = """[input]
system=
{system}
charge=0
runtype=energy
method=hf
basis={basis}
ispher=true
[guess]
type={guess}
[scf]
type=rhf
multiplicity=1
conv=1.0e-10
maxit=200
[symmetry]
enabled=false
"""

IF_SYSTEM = """   53   0.000000000   0.000000000   0.000000000
    9   0.100000000   0.200000000   1.910000000"""
HBR_SYSTEM = """   35   0.000000000   0.000000000   0.000000000
    1   0.100000000   0.200000000   1.395000000"""

NCART = {"S": 1, "P": 3, "D": 6, "F": 10}


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


def _read_table():
    lines = [ln for ln in DATA.read_text().splitlines() if ln and not ln.startswith("#")]
    toks = " ".join(lines).split()
    zmax, pos, table = int(toks[0]), 1, {}
    for _ in range(zmax):
        z, nao, nelec = int(toks[pos]), int(toks[pos + 1]), float(toks[pos + 2])
        pos += 3 + nao * nao
        table[z] = (nao, nelec)
    return zmax, table


def _sto3g_sizes():
    body = STO3G.read_text().split("$DATA", 1)[1]
    sizes, z = [], None
    for raw in body.splitlines():
        line = raw.strip()
        if line.startswith("!") or line.startswith("$END"):
            continue
        if not line:
            z = None
            continue
        tok = line.split()
        if z is None:
            z = len(sizes) + 1
            sizes.append(0)
        elif tok[0].upper() in NCART and len(tok) == 2:
            sizes[-1] += NCART[tok[0].upper()]
    return sizes


class MinaoTable(unittest.TestCase):
    def test_table_covers_sto3g_range(self):
        zmax, table = _read_table()
        sizes = _sto3g_sizes()
        self.assertEqual(zmax, 54)
        self.assertEqual(len(sizes), zmax)
        for z in range(1, zmax + 1):
            nao, nelec = table[z]
            self.assertEqual(nao, sizes[z - 1], f"Z={z}")
            self.assertAlmostEqual(nelec, z, places=8, msg=f"Z={z}")


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class MinaoEcpGuess(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory(prefix="oqp_minao_ecp_")
        cls.workdir = Path(cls._tmp.name)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def _run(self, name, system, basis, guess):
        from oqp.pyoqp import Runner

        tag = f"{name}_{guess}"
        inp = self.workdir / f"{tag}.inp"
        log = self.workdir / f"{tag}.log"
        inp.write_text(INPUT_TMPL.format(system=system, basis=basis, guess=guess))
        runner = Runner(project=tag, input_file=str(inp), log=str(log), silent=1, usempi=False)
        runner.run()
        return float(runner.mol.mol_energy.energy), log.read_text()

    def _atom_electrons(self, text):
        """(Z, ECP core, electrons) rows printed by the MINAO guess."""
        block = text.split("Electrons in atomic density", 1)[1]
        rows = []
        for line in block.splitlines()[1:]:
            m = re.match(r"\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s*$", line)
            if not m:
                break
            rows.append((int(m.group(2)), int(m.group(3)), float(m.group(4))))
        return rows

    def _check(self, name, system, basis, expected):
        e_minao, text = self._run(name, system, basis, "minao")
        e_huckel, _ = self._run(name, system, basis, "huckel")
        self.assertEqual(len(self._atom_electrons(text)), len(expected))
        for (z, ncore, nel), (z_ref, ncore_ref) in zip(self._atom_electrons(text), expected):
            self.assertEqual((z, ncore), (z_ref, ncore_ref))
            self.assertAlmostEqual(nel, z - ncore, places=5)
        self.assertAlmostEqual(e_minao, e_huckel, delta=1.0e-8)

    def test_iodine_def2_ecp(self):
        self._check("if_def2", IF_SYSTEM, "def2-svp", [(53, 28), (9, 0)])

    def test_bromine_lanl2dz_ecp(self):
        self._check("hbr_lanl2", HBR_SYSTEM, "lanl2dz", [(35, 28), (1, 0)])

    def test_iodine_all_electron(self):
        self._check("if_321g", IF_SYSTEM, "3-21g", [(53, 0), (9, 0)])


if __name__ == "__main__":
    unittest.main()
