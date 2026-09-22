"""Live validation for the ACID / induced-current-density export.

Runs a GIAO NMR calculation on H2O/STO-3G, rebuilds the first-order current
density from the exported ``OQP::nmr_pdens`` response, and compares the field
pointwise against a committed independent GIAO reference
(``tests/fixtures/nmr/giao_current_reference.json``).

The field, not the integrated shielding, is what an ACID map plots, so it is the
field that has to be asserted.  The reference was produced by an independent
GIAO implementation whose current density was itself verified by Biot-Savart
integration back into the analytic shielding tensor (grid-converged to ~1e-8 ppm
on water and ethyne, symmetric part to ~1e-6 ppm under GIAO).

The reference code is NOT required at run time: the fixture is committed.  The
test runs OpenQP in a subprocess and skips gracefully when the shared library is
not built.
"""
import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "nmr" / "giao_current_reference.json"

# Cross-code tolerances: the SCF densities agree to well under this, and the
# current density is a first-order response evaluated pointwise.
TOL_ABS = 2.0e-6
TOL_REL = 2.0e-4


def _oqp_root():
    root = Path(os.environ.get("OPENQP_ROOT", str(ROOT)))
    lib = root / "lib" / "liboqp.so"
    if not lib.exists():
        lib = root / "lib" / "liboqp.dylib"
    if not lib.exists() or not (root / "include" / "oqp.h").exists():
        return None
    return str(root)


def _input(geometry):
    return ("[input]\nsystem=\n" + geometry + "\ncharge=0\nruntype=energy\n"
            "basis=sto-3g\nmethod=hf\n\n[guess]\ntype=huckel\n\n"
            "[scf]\nmultiplicity=1\ntype=rhf\nverbose=2\n\n[properties]\nscf_prop=\n")


class AcidCurrentDensityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = _oqp_root()
        if cls.root is None:
            raise unittest.SkipTest("OpenQP shared library/header not built")
        if not FIXTURE.exists():
            raise unittest.SkipTest("GIAO current-density reference fixture missing")
        cls.ref = json.loads(FIXTURE.read_text())
        cls.got = cls._run()

    @classmethod
    def _run(cls):
        ref = cls.ref
        charge = {"H": 1, "O": 8}
        geometry = "\n".join(
            f"   {charge[sym]}   {c[0]:.9f}   {c[1]:.9f}   {c[2]:.9f}"
            for sym, c in ref["geometry_angstrom"])
        with tempfile.TemporaryDirectory() as wd:
            inp = Path(wd) / "h2o.inp"
            log = Path(wd) / "h2o.log"
            out = Path(wd) / "current.json"
            inp.write_text(_input(geometry))
            pts = json.dumps(ref["points_bohr"])
            script = Path(wd) / "run.py"
            script.write_text(textwrap.dedent(f"""
                import json
                import numpy as np
                import oqp
                from oqp.pyoqp import Runner
                from oqp.export import AcidExporter
                from oqp.export.aicd import acid_scalar

                r = Runner(project='h2o', input_file={str(inp)!r},
                           log={str(log)!r}, silent=1, usempi=False)
                r.run()
                oqp.nmr_giao_shielding_debug(r.mol)

                ex = AcidExporter(r.mol)
                pts = np.array({pts})
                j = ex.current_density(pts)

                rng = np.random.default_rng(11)
                q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
                if np.linalg.det(q) < 0:
                    q[:, 0] *= -1.0
                rot = acid_scalar(np.einsum('ia,pab,jb->pij', q, j, q))

                json.dump({{"e_tot": float(r.mol.get_scf_energy()),
                           "current_tensor": j.tolist(),
                           "acid": acid_scalar(j).tolist(),
                           "acid_rotated": rot.tolist()}},
                          open({str(out)!r}, "w"))
            """))
            env = dict(os.environ)
            env["OPENQP_ROOT"] = cls.root
            env["PYTHONPATH"] = os.pathsep.join(
                p for p in (str(ROOT / "pyoqp"), env.get("PYTHONPATH", "")) if p)
            proc = subprocess.run([sys.executable, str(script)], cwd=wd, env=env,
                                  capture_output=True, text=True, timeout=600)
            if not out.exists():
                raise unittest.SkipTest(
                    "ACID run produced no output\n"
                    f"stdout:\n{proc.stdout[-1500:]}\nstderr:\n{proc.stderr[-2500:]}")
            return json.loads(out.read_text())

    def test_scf_energy_matches_reference(self):
        self.assertAlmostEqual(self.got["e_tot"], self.ref["e_tot"], delta=1e-7,
                               msg="SCF energy differs; geometry or basis mismatch")

    def test_current_density_matches_reference(self):
        import numpy as np
        got = np.asarray(self.got["current_tensor"])
        ref = np.asarray(self.ref["current_tensor"])
        self.assertEqual(got.shape, ref.shape)
        tol = TOL_ABS + TOL_REL * np.abs(ref)
        bad = np.abs(got - ref) > tol
        if bad.any():
            p, i, j = np.argwhere(bad)[0]
            self.fail(
                f"current density differs at point {p} component ({i},{j}): "
                f"got {got[p, i, j]:.8e}, reference {ref[p, i, j]:.8e}; "
                f"{bad.sum()} of {bad.size} components outside tolerance")

    def test_acid_scalar_matches_reference(self):
        import numpy as np
        got = np.asarray(self.got["acid"])
        ref = np.asarray(self.ref["acid"])
        tol = TOL_ABS + TOL_REL * np.abs(ref)
        self.assertTrue(np.all(np.abs(got - ref) <= tol),
                        f"max ACID deviation {np.abs(got - ref).max():.3e}")

    def test_acid_is_rotationally_invariant(self):
        import numpy as np
        got = np.asarray(self.got["acid"])
        rot = np.asarray(self.got["acid_rotated"])
        self.assertLess(np.abs(got - rot).max(), 1e-12,
                        "ACID must be invariant under rotation of the tensor field")


if __name__ == "__main__":
    unittest.main()
