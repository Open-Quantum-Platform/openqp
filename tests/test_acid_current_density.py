"""Live validation for the ACID / induced-current-density export.

Runs GIAO NMR calculations, rebuilds the first-order current density from the
exported ``OQP::nmr_pdens`` response, and compares the field pointwise against a
committed independent GIAO reference
(``tests/fixtures/nmr/giao_current_reference.json``).

The field, not the integrated shielding, is what an ACID map plots, so it is the
field that has to be asserted.  The reference was produced by an independent
GIAO implementation whose current density was itself verified by Biot-Savart
integration back into the analytic shielding tensor (grid-converged to ~1e-8 ppm
on water and ethyne, symmetric part to ~1e-6 ppm under GIAO).

Three things beyond the numbers are covered here, all in one process:

* ``OQP::nmr_pdens`` is a persistent tagarray replaced on every GIAO call, so
  the response is recomputed at a larger basis and then back at the smaller one
  and checked for both the resized shape and a value that is not stale.
* the response is AO-indexed, so a spherical basis (cc-pVDZ) is exercised for
  shape and antisymmetry even though the grid evaluator is Cartesian-only.
* that Cartesian-only limit must be an explicit refusal, not a partial map.

The reference code is NOT required at run time: the fixture is committed.  The
test skips only when the shared library is not built; once OpenQP runs, a
failure is a failure.
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


def _input(geometry, basis):
    return ("[input]\nsystem=\n" + geometry + "\ncharge=0\nruntype=energy\n"
            f"basis={basis}\nmethod=hf\n\n[guess]\ntype=huckel\n\n"
            "[scf]\nmultiplicity=1\ntype=rhf\nverbose=2\n\n[properties]\nscf_prop=\n")


DRIVER = """
import json
import numpy as np
import oqp
from oqp.pyoqp import Runner
from oqp.export import AcidExporter
from oqp.export.aicd import acid_scalar

POINTS = np.array({points})
RESULT = {{}}


def giao(tag, inp, log):
    r = Runner(project=tag, input_file=inp, log=log, silent=1, usempi=False)
    r.run()
    oqp.nmr_giao_shielding_debug(r.mol)
    return r.mol


def response(mol):
    nbf = int(round(np.asarray(mol.data["OQP::VEC_MO_A"]).size ** 0.5))
    raw = np.array(mol.data["OQP::nmr_pdens"], copy=True).ravel(order="C")
    return nbf, raw.reshape((3, nbf, nbf), order="F")


# --- small basis, the numerical reference case ------------------------------
mol = giao("small", {inp_small!r}, {log_small!r})
nbf_small, p_small = response(mol)
exporter = AcidExporter(mol)
j = exporter.current_density(POINTS)
rng = np.random.default_rng(11)
q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
if np.linalg.det(q) < 0:
    q[:, 0] *= -1.0
RESULT["e_tot"] = float(mol.get_scf_energy())
RESULT["current_tensor"] = j.tolist()
RESULT["acid"] = acid_scalar(j).tolist()
RESULT["acid_rotated"] = acid_scalar(np.einsum("ia,pab,jb->pij", q, j, q)).tolist()
RESULT["nbf_small"] = nbf_small

# --- same process, larger basis: the buffer must grow -----------------------
mol_big = giao("big", {inp_big!r}, {log_big!r})
nbf_big, p_big = response(mol_big)
RESULT["nbf_big"] = nbf_big
RESULT["shape_big"] = list(p_big.shape)
RESULT["antisymmetry_big"] = float(np.abs(p_big + np.swapaxes(p_big, 1, 2)).max())

# --- same process, back to the small basis: it must shrink and be fresh -----
mol_again = giao("again", {inp_small!r}, {log_again!r})
nbf_again, p_again = response(mol_again)
RESULT["nbf_again"] = nbf_again
RESULT["shape_again"] = list(p_again.shape)
RESULT["stale_after_shrink"] = float(np.abs(p_again - p_small).max())
RESULT["scale_small"] = float(np.abs(p_small).max())

# --- spherical basis: AO-indexed response, Cartesian-only exporter ----------
mol_sph = giao("sph", {inp_sph!r}, {log_sph!r})
nbf_sph, p_sph = response(mol_sph)
RESULT["nbf_spherical"] = nbf_sph
RESULT["shape_spherical"] = list(p_sph.shape)
RESULT["antisymmetry_spherical"] = float(np.abs(p_sph + np.swapaxes(p_sph, 1, 2)).max())
RESULT["scale_spherical"] = float(np.abs(p_sph).max())
try:
    AcidExporter(mol_sph)
except NotImplementedError as error:
    RESULT["spherical_refusal"] = str(error)
except Exception as error:  # noqa: BLE001 - recorded so the test can report it
    RESULT["spherical_refusal"] = f"{{type(error).__name__}}: {{error}}"
else:
    RESULT["spherical_refusal"] = ""

json.dump(RESULT, open({out!r}, "w"))
"""


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
            work = Path(wd)
            paths = {}
            for tag, basis in (("small", ref["basis"]), ("big", "6-31G"),
                               ("sph", "cc-pvdz")):
                inp = work / f"{tag}.inp"
                inp.write_text(_input(geometry, basis))
                paths[tag] = str(inp)
            out = work / "acid.json"
            script = work / "run.py"
            script.write_text(DRIVER.format(
                points=json.dumps(ref["points_bohr"]),
                inp_small=paths["small"], log_small=str(work / "small.log"),
                inp_big=paths["big"], log_big=str(work / "big.log"),
                log_again=str(work / "again.log"),
                inp_sph=paths["sph"], log_sph=str(work / "sph.log"),
                out=str(out)))
            env = dict(os.environ)
            env["OPENQP_ROOT"] = cls.root
            env["PYTHONPATH"] = os.pathsep.join(
                p for p in (str(ROOT / "pyoqp"), env.get("PYTHONPATH", "")) if p)
            proc = subprocess.run([sys.executable, str(script)], cwd=wd, env=env,
                                  capture_output=True, text=True, timeout=1200)
            if proc.returncode != 0 or not out.exists():
                raise AssertionError(
                    f"the ACID driver failed (exit {proc.returncode})\n"
                    f"stdout:\n{proc.stdout[-2000:]}\nstderr:\n{proc.stderr[-4000:]}")
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

    def test_response_buffer_grows_and_shrinks_in_one_process(self):
        self.assertGreater(self.got["nbf_big"], self.got["nbf_small"],
                           "the larger basis must actually be larger")
        self.assertEqual(self.got["shape_big"],
                         [3, self.got["nbf_big"], self.got["nbf_big"]],
                         "OQP::nmr_pdens did not grow with the basis")
        self.assertEqual(self.got["nbf_again"], self.got["nbf_small"])
        self.assertEqual(self.got["shape_again"],
                         [3, self.got["nbf_small"], self.got["nbf_small"]],
                         "OQP::nmr_pdens did not shrink back with the basis")
        self.assertLess(self.got["stale_after_shrink"],
                        1e-10 * max(self.got["scale_small"], 1.0),
                        "the response after shrinking differs from the first run "
                        "at the same basis, so the buffer carried stale values")

    def test_response_is_antisymmetric_in_a_spherical_basis(self):
        self.assertEqual(self.got["shape_spherical"],
                         [3, self.got["nbf_spherical"], self.got["nbf_spherical"]])
        scale = max(self.got["scale_spherical"], 1.0)
        self.assertLess(self.got["antisymmetry_spherical"], 1e-10 * scale,
                        "the magnetic density response must be antisymmetric")
        self.assertLess(self.got["antisymmetry_big"],
                        1e-10 * max(self.got["scale_small"], 1.0))

    def test_acid_uses_the_published_normalisation(self):
        """Pin the scale: 0.05 a.u. has to mean what the literature means.

        Recomputed here straight from Herges and Geuenich (J. Phys. Chem. A 105,
        3214), in the form the GIMIC reference implementation uses, so a change
        of convention in ``acid_scalar`` fails instead of silently moving every
        isosurface.  Uses the committed tensors, so it needs no calculation.
        """
        import numpy as np
        t = np.asarray(self.ref["current_tensor"])
        diag = ((t[:, 0, 0] - t[:, 1, 1]) ** 2 + (t[:, 1, 1] - t[:, 2, 2]) ** 2
                + (t[:, 2, 2] - t[:, 0, 0]) ** 2)
        off = ((t[:, 0, 1] + t[:, 1, 0]) ** 2 + (t[:, 0, 2] + t[:, 2, 0]) ** 2
               + (t[:, 1, 2] + t[:, 2, 1]) ** 2)
        published = np.sqrt(diag / 3.0 + off / 2.0)
        # Check the fixture against its own tensors, so this pins the convention
        # to machine precision; OpenQP is held to the fixture separately, where
        # the cross-code tolerance belongs.
        self.assertLess(np.abs(np.asarray(self.ref["acid"]) - published).max(),
                        1e-12, "the reference fixture is not on the published scale")

    def test_spherical_basis_is_refused_explicitly(self):
        message = self.got["spherical_refusal"]
        self.assertTrue(message,
                        "a spherical basis must be refused, not silently mapped")
        self.assertIn("Cartesian", message)


if __name__ == "__main__":
    unittest.main()
