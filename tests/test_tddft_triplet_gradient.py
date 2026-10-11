"""Analytic closed-shell triplet RPA/TDA gradients (tdhf_z_vector, tdhf_gradient).

For a closed-shell reference the triplet response has no Coulomb term in its
transition-density parts and the f_aa - f_ab kernel; the gradient code builds
H+[X+Y] without Coulomb, uses the spin-resolved f_xc / g_xc contractions with
dxb = -dxa, and drops the (X+Y)(X+Y) Coulomb term of the two-electron
derivative.  Checks:

* singlet routed through the same unrestricted code (OQP_TDDFT_XC_UNRESTRICTED=1,
  dxb = +dxa) equals the restricted singlet gradient: pins the normalization
  of the spin-resolved paths the triplet uses;
* triplet gradient == central finite differences of the triplet energy
  (C1 H2O/6-31G*, HF, SVWN, PBE, B3LYP, CAM-B3LYP; TDA and RPA).
Meta-GGA triplet gradients are rejected by the input checker (their
spin-antisymmetric tau terms do not reproduce finite differences); the
singlet identity above still covers M06-2X.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BOHR = 0.52917721092
GEOM = [("O", 0.010000, -0.020000, 0.117300),
        ("H", 0.080000, 0.770200, -0.449200),
        ("H", -0.050000, -0.731200, -0.509200)]
STEP = 1.0e-3          # Angstrom
COMPONENTS = [(0, 2), (1, 0), (2, 1)]   # (atom, xyz) checked by finite differences


def _runtime_available():
    try:
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except (ImportError, OSError, RuntimeError):
        return False


def _deck(geom, functional, td_type, mult, runtype):
    system = "".join(f"   {a} {x:.8f} {y:.8f} {z:.8f}\n" for a, x, y, z in geom)
    s = (f"[input]\nsystem=\n{system}charge=0\nruntype={runtype}\nbasis=6-31g*\nmethod=tdhf\n")
    if functional:
        s += f"functional={functional}\n"
    s += ("\n[guess]\ntype=huckel\nsave_mol=true\n\n"
          "[scf]\ntype=rhf\nmultiplicity=1\nmaxit=200\nconv=1e-11\n\n"
          f"[tdhf]\ntype={td_type}\nnstate=3\nmaxit=200\nconv=1e-11\nzvconv=1e-12\nmultiplicity={mult}\n")
    if runtype == "grad":
        s += "\n[properties]\ngrad=1\n"
    if functional:
        s += "\n[dftgrid]\nrad_npts=96\nang_npts=302\n"
    return s


def _run(text, env_extra=None):
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "job.inp").write_text(text)
        env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "2"))
        env.pop("OQP_TDDFT_XC_UNRESTRICTED", None)
        env.update(env_extra or {})
        proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", "job.inp"], cwd=tmp, env=env,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=900)
        log = (Path(tmp) / "job.log").read_text(errors="ignore") if (Path(tmp) / "job.log").exists() else ""
        jpath = Path(tmp) / "job.json"
        res = json.loads(jpath.read_text()) if jpath.exists() else None
        return proc.returncode, res, log[-2500:]


def _gradient(functional, td_type, mult, env_extra=None):
    rc, res, log = _run(_deck(GEOM, functional, td_type, mult, "grad"), env_extra)
    return rc, (res or {}).get("grad"), log


def _state_energy(geom, functional, td_type, mult):
    rc, res, log = _run(_deck(geom, functional, td_type, mult, "energy"))
    return rc, ((res or {}).get("energies") or [None, None])[1], log


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestTripletGradient(unittest.TestCase):

    def test_unrestricted_path_reproduces_the_restricted_singlet_gradient(self):
        for func in ("b3lypv5", "pbe", "m06-2x"):
            for td in ("tda", "rpa"):
                with self.subTest(functional=func, td=td):
                    r0, g0, d0 = _gradient(func, td, 1)
                    r1, g1, d1 = _gradient(func, td, 1, {"OQP_TDDFT_XC_UNRESTRICTED": "1"})
                    self.assertEqual(r0, 0, d0)
                    self.assertEqual(r1, 0, d1)
                    self.assertIsNotNone(g0, d0)
                    for a, b in zip(g0, g1):
                        self.assertAlmostEqual(a, b, delta=1e-7)

    def test_triplet_gradient_equals_finite_differences(self):
        for func in (None, "svwn", "pbe", "b3lypv5", "cam-b3lyp"):
            for td in ("tda", "rpa"):
                with self.subTest(functional=func or "hf", td=td):
                    rc, grad, dlog = _gradient(func, td, 3)
                    self.assertEqual(rc, 0, dlog)
                    self.assertIsNotNone(grad, dlog)
                    tol = 2e-5 if func else 1e-5
                    for atom, xyz in COMPONENTS:
                        energies = []
                        for sign in (1, -1):
                            geom = [list(g) for g in GEOM]
                            geom[atom][1 + xyz] += sign * STEP
                            rc_e, e, elog = _state_energy(geom, func, td, 3)
                            self.assertEqual(rc_e, 0, elog)
                            self.assertIsNotNone(e, elog)
                            energies.append(e)
                        fd = (energies[0] - energies[1]) / (2 * STEP / BOHR)
                        self.assertAlmostEqual(grad[3 * atom + xyz], fd, delta=tol,
                                               msg=f"atom {atom} xyz {xyz}")


if __name__ == "__main__":
    unittest.main()
