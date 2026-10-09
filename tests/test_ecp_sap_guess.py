"""SAP initial guess with an effective core potential: IF / def2-SVP.

The def2 iodine ECP removes 28 core electrons.  The SAP guess builds
F = T + V_SAP, and the SAP table stops at Z = 36, so iodine contributed no
potential at all (and the ECP was never added): the guess had no attraction to
the iodine nucleus.  RHF then converged to a
symmetry-broken saddle point 0.276 Eh above the ground state: the I 4p and 4d
pi/delta pairs were split by ~2.5e-3 Eh, and because rigid rotations of the
molecule landed on differently oriented broken solutions, the static
polarizability was not rotationally covariant (errors up to 0.26 a.u.).

Elements beyond the table now get the bare (ECP-screened) attraction
-Z_val/r, tabulated ECP atoms have Z_eff capped at Z_val, and the ECP is added
to the guess Fock.  Checked here, with integral symmetry disabled:
  * SAP and Hueckel guesses reach the same RHF energy;
  * the occupied spectrum of the linear molecule has exactly degenerate pi
    pairs (7 sigma + 5 pi pairs for 17 occupied orbitals);
  * alpha has the bond axis as an eigenvector with two equal perpendicular
    eigenvalues, and alpha(R x) = R alpha(x) R^T for a rigid rotation R.
The alpha tolerances (5e-4 a.u., alpha ~ 27) sit at the CPHF solver level:
HCl/def2-SVP and all-electron IF/3-21G show 2e-5 to 8e-5 asymmetry and
rotation error; the broken-symmetry solution gave 0.16.

HBr/LANL2DZ covers a tabulated ECP atom (Br, Z = 35, 28 core electrons),
where Z_eff is capped at Z_val: SAP and Hueckel give the same energy and the
pi pair of the valence shell is degenerate.

Skipped unless the compiled OpenQP runtime is importable.
"""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

HBR_INPUT = """[input]
system=
   35   0.000000000   0.000000000   0.000000000
    1   0.100000000   0.200000000   1.395000000
charge=0
runtype=energy
method=hf
basis=lanl2dz
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

INPUT_TMPL = """[input]
system=
{system}
charge=0
runtype=energy
method=hf
basis=def2-svp
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

# I at the origin, F off-axis so that no coordinate axis is the bond axis
COORDS = np.array([[0.0, 0.0, 0.0], [0.1, 0.2, 1.91]])  # Angstrom
ZNUC = (53, 9)
# 25 I valence (def2 ECP, 28 core e-) + 9 F electrons
NOCC = 17
# CPHF-solver level for alpha (a.u.)
ALPHA_TOL = 5.0e-4


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


def _rotation(axis, theta):
    a = np.asarray(axis, dtype=float)
    a /= np.linalg.norm(a)
    k = np.array([[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]])
    return np.eye(3) + np.sin(theta) * k + (1.0 - np.cos(theta)) * k @ k


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class EcpSapGuess(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory(prefix="oqp_ecp_sap_")
        cls.workdir = Path(cls._tmp.name)
        cls.rot = _rotation([1.0, 1.0, 1.0], 0.5)
        cls.ref = cls._run("sap", COORDS)
        cls.huckel = cls._run("huckel", COORDS)
        cls.rotated = cls._run("sap", COORDS @ cls.rot.T)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    @classmethod
    def _run(cls, guess, coords):
        import oqp
        from oqp.pyoqp import Runner

        tag = f"{guess}_{len(list(cls.workdir.iterdir()))}"
        system = "\n".join(f"   {z}  {x:.12f}  {y:.12f}  {w:.12f}"
                           for z, (x, y, w) in zip(ZNUC, coords))
        inp = cls.workdir / f"{tag}.inp"
        inp.write_text(INPUT_TMPL.format(system=system, guess=guess))
        runner = Runner(project=f"if_{tag}", input_file=str(inp),
                        log=str(cls.workdir / f"{tag}.log"), silent=1, usempi=False)
        runner.run()
        mol = runner.mol
        alpha = np.zeros((3, 3), dtype=np.float64)
        oqp.cphf_static_polarizability(mol, oqp.ffi.cast("double *", oqp.ffi.from_buffer(alpha)))
        return {
            "energy": float(mol.mol_energy.energy),
            "mo": np.array(mol.data["OQP::E_MO_A"], dtype=float)[:NOCC],
            "alpha": alpha,
        }

    def test_sap_reaches_ground_state(self):
        self.assertAlmostEqual(self.ref["energy"], self.huckel["energy"], delta=1.0e-8)

    def test_pi_pairs_degenerate(self):
        mo = np.sort(self.ref["mo"])
        self.assertEqual(mo.size, NOCC)
        gaps = np.diff(mo)
        # 17 occupied = 7 sigma + 5 doubly degenerate pi pairs
        self.assertEqual(int(np.sum(gaps < 1.0e-6)), 5)
        self.assertGreater(float(np.min(gaps[gaps >= 1.0e-6])), 1.0e-4)

    def test_alpha_cylindrical(self):
        axis = COORDS[1] - COORDS[0]
        axis /= np.linalg.norm(axis)
        alpha = self.ref["alpha"]
        self.assertLess(np.max(np.abs(alpha - alpha.T)), ALPHA_TOL)
        self.assertLess(np.linalg.norm(alpha @ axis - (axis @ alpha @ axis) * axis), ALPHA_TOL)
        w = np.linalg.eigvalsh(alpha)
        a_par = axis @ alpha @ axis
        perp = sorted(w, key=lambda v: abs(v - a_par))[1:]
        self.assertLess(abs(perp[0] - perp[1]), ALPHA_TOL)

    def test_alpha_rotationally_covariant(self):
        r = self.rot
        self.assertAlmostEqual(self.rotated["energy"], self.ref["energy"], delta=1.0e-8)
        err = np.max(np.abs(self.rotated["alpha"] - r @ self.ref["alpha"] @ r.T))
        self.assertLess(err, ALPHA_TOL)


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class EcpSapGuessTabulated(unittest.TestCase):
    def _run(self, guess, workdir):
        from oqp.pyoqp import Runner

        inp = Path(workdir) / f"hbr_{guess}.inp"
        inp.write_text(HBR_INPUT.format(guess=guess))
        runner = Runner(project=f"hbr_{guess}", input_file=str(inp),
                        log=str(Path(workdir) / f"hbr_{guess}.log"), silent=1, usempi=False)
        runner.run()
        # 7 Br valence + 1 H electrons: sigma, sigma, pi pair
        mo = np.sort(np.array(runner.mol.data["OQP::E_MO_A"], dtype=float)[:4])
        return float(runner.mol.mol_energy.energy), mo

    def test_sap_matches_huckel_with_capped_zeff(self):
        with tempfile.TemporaryDirectory(prefix="oqp_ecp_sap_hbr_") as tmp:
            e_sap, mo = self._run("sap", tmp)
            e_huckel, _ = self._run("huckel", tmp)
        self.assertAlmostEqual(e_sap, e_huckel, delta=1.0e-8)
        self.assertEqual(int(np.sum(np.diff(mo) < 1.0e-6)), 1)


if __name__ == "__main__":
    unittest.main()
