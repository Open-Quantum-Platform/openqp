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


def basis_nbf(m):
    return int(m.data.get_basis()["nbf"])


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

# --- the workflow path, not just the exporter API --------------------------
# scf_prop=acid dispatches here, and the shipped example asserts nothing about
# what comes out of it -- it passes as long as nothing raises.  So check the
# files themselves: four cubes, a grid that matches the header, and a scalar
# field that is a norm and therefore cannot be negative.
import os as _os
from oqp.library.runfunc import write_acid_cubes

mol.config["properties"]["acid_spacing"] = "0.6"
mol.config["properties"]["acid_padding"] = "2.0"
cube_paths = write_acid_cubes(mol)
RESULT["cube_names"] = [_os.path.basename(q) for q in cube_paths]
RESULT["cubes_exist"] = [_os.path.exists(q) for q in cube_paths]


def read_cube(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    natom = int(lines[2].split()[0])
    dims = [int(lines[3 + k].split()[0]) for k in range(3)]
    body = " ".join(lines[6 + natom:]).split()
    return natom, dims, np.array(body, dtype=float)


natom_c, dims_c, vals_c = read_cube(cube_paths[0])
RESULT["cube_natom"] = natom_c
RESULT["cube_npts_match"] = bool(vals_c.size == dims_c[0] * dims_c[1] * dims_c[2])
RESULT["cube_all_finite"] = bool(np.all(np.isfinite(vals_c)))
RESULT["cube_acid_min"] = float(vals_c.min())
RESULT["cube_acid_max"] = float(vals_c.max())
RESULT["cube_vector_sizes"] = [int(read_cube(q)[2].size) for q in cube_paths[1:]]

# --- the SAME molecule and tagarray, at a different basis size -------------
# A fresh Runner would build a fresh OQPData and never replace anything, so the
# basis is changed on the molecule that already owns OQP::nmr_pdens.
import oqp.library
from oqp.library.single_point import SinglePoint


def rebasis(m, basis):
    m.config["input"]["basis"] = basis
    oqp.library.set_basis(m)
    SinglePoint(m).energy()
    oqp.nmr_giao_shielding_debug(m)
    return response(m)


nbf_big, p_big = rebasis(mol, "6-31G")
RESULT["nbf_big"] = nbf_big
RESULT["shape_big"] = list(p_big.shape)
RESULT["antisymmetry_big"] = float(np.abs(p_big + np.swapaxes(p_big, 1, 2)).max())

nbf_again, p_again = rebasis(mol, {basis_small!r})
RESULT["nbf_again"] = nbf_again
RESULT["shape_again"] = list(p_again.shape)
RESULT["stale_after_shrink"] = float(np.abs(p_again - p_small).max())
RESULT["scale_small"] = float(np.abs(p_small).max())

# --- the response is persistent, so a molecule that moved must be refused ---
# The geometry-optimisation shape of the defect: OQP::nmr_pdens survives
# update_system() untouched and, at an unchanged basis size, still loads with
# the right shape -- so without a provenance stamp the exporter silently
# combines the previous geometry's response with this geometry's coordinates.
AcidExporter(mol)  # current right here; everything below deliberately is not
home = np.asarray(mol.get_system(), dtype=float).copy()
away = home.copy()
away[0] += 0.05
mol.update_system(away)
try:
    AcidExporter(mol)
except ValueError as error:
    RESULT["moved_refusal"] = str(error)
else:
    RESULT["moved_refusal"] = ""

# A new SCF replaces the orbitals the response was built from, and nothing in
# a stamp can identify those, so scf_driver drops the response instead.  Put
# the coordinates back first, so the geometry check cannot be what answers.
mol.update_system(home)
SinglePoint(mol).energy()
try:
    AcidExporter(mol)
except ValueError as error:
    RESULT["restale_refusal"] = str(error)
else:
    RESULT["restale_refusal"] = ""

# The two steps above deliberately left the molecule stale; re-converge and
# re-run the GIAO so the basis checks below start from a response that really
# is current, rather than from one the density check would reject anyway.
SinglePoint(mol).energy()
oqp.nmr_giao_shielding_debug(mol)
AcidExporter(mol)

# --- the AO ordering is part of the response's identity --------------------
# Hold the density and the geometry fixed and permute the stamp's shell block,
# which is what the exporter sees when a same-size basis comes back with its
# shells reordered.  nbf is unchanged, and no invariant of the density could
# have helped either -- trace and sum of squares both survive D -> P D P^T --
# so only the shell block can catch this.
snap = np.array(mol.data["OQP::nmr_pdens_ref"], copy=True).ravel()
buf = np.asarray(mol.data["OQP::nmr_pdens_ref"]).ravel()
nat = len(mol.get_atoms())
head = 8 + 4 * nat          # header, charge/electron counts, Z, coords
nsh = int(round(snap[2]))
nprim = int(round(snap[3]))
blk = buf[head:head + 3 * nsh].reshape(nsh, 3).copy()
RESULT["permuted_rows_differ"] = bool(np.any(blk[0] != blk[2]))
blk[[0, 2]] = blk[[2, 0]]
buf[head:head + 3 * nsh] = blk.ravel()
try:
    AcidExporter(mol)
except ValueError as error:
    RESULT["permuted_refusal"] = str(error)
else:
    RESULT["permuted_refusal"] = ""
buf[:] = snap
AcidExporter(mol)  # the restore has to leave a usable response behind

# Nuclear and charge identity: the coordinates and the basis are untouched, so
# only these slots say which system the response belongs to.  Poking them is
# what the exporter sees when a composition or a charge is edited in place with
# no SCF in between.
buf[4] += 1.0
try:
    AcidExporter(mol)
except ValueError as error:
    RESULT["charge_refusal"] = str(error)
else:
    RESULT["charge_refusal"] = ""
buf[:] = snap
buf[8] += 1.0
try:
    AcidExporter(mol)
except ValueError as error:
    RESULT["z_refusal"] = str(error)
else:
    RESULT["z_refusal"] = ""
buf[:] = snap

# The sentinel the driver writes when it allocates the buffer and clears only
# once the response is complete.  Nothing else reaches this branch: it stands
# for a GIAO run that died between the two.
buf[0] = -1.0
try:
    AcidExporter(mol)
except ValueError as error:
    RESULT["incomplete_refusal"] = str(error)
else:
    RESULT["incomplete_refusal"] = ""
buf[:] = snap

# Primitives that share a per-shell summary: keep the centres, angular momenta
# and contraction counts, and move two exponents within one shell to their
# mean.  The sum is preserved exactly, so any per-shell reduction of the
# primitives is blind to this -- only storing them in full is not.
pbase = head + 3 * nsh
ncontr = buf[head:head + 3 * nsh].reshape(nsh, 3)[:, 2].astype(int)
off, target = 0, -1
for nc in ncontr:
    if nc >= 2 and buf[pbase + off] != buf[pbase + off + 1]:
        target = off
        break
    off += nc
RESULT["prim_target_found"] = bool(target >= 0)
if target >= 0:
    e0, e1 = float(buf[pbase + target]), float(buf[pbase + target + 1])
    mid = 0.5 * (e0 + e1)
    buf[pbase + target] = mid
    buf[pbase + target + 1] = mid
    RESULT["prim_sum_drift"] = float(abs((mid + mid) - (e0 + e1)))
    RESULT["prim_values_moved"] = float(abs(mid - e0))
    try:
        AcidExporter(mol)
    except ValueError as error:
        RESULT["prim_refusal"] = str(error)
    else:
        RESULT["prim_refusal"] = ""
    buf[:] = snap
    AcidExporter(mol)

# A real same-size basis swap: STO-6G has the same seven AOs as STO-3G here, so
# nbf cannot tell them apart.  No SCF here on purpose -- an SCF would erase the
# response and answer before the stamp ever got the chance.
mol.config["input"]["basis"] = "sto-6g"
oqp.library.set_basis(mol)
RESULT["rebasis_nbf"] = int(basis_nbf(mol))
try:
    AcidExporter(mol)
except ValueError as error:
    RESULT["rebasis_refusal"] = str(error)
else:
    RESULT["rebasis_refusal"] = ""

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
            for tag, basis in (("small", ref["basis"]), ("sph", "cc-pvdz")):
                inp = work / f"{tag}.inp"
                inp.write_text(_input(geometry, basis))
                paths[tag] = str(inp)
            out = work / "acid.json"
            script = work / "run.py"
            script.write_text(DRIVER.format(
                points=json.dumps(ref["points_bohr"]),
                inp_small=paths["small"], log_small=str(work / "small.log"),
                basis_small=ref["basis"],
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

    def test_a_moved_geometry_invalidates_the_response(self):
        """A response left over from the previous geometry must not be plotted."""
        msg = self.got["moved_refusal"]
        self.assertTrue(
            msg, "the exporter accepted a response from a geometry that has "
                 "since moved, and would have written a map mixing the two")
        self.assertIn("geometry moved", msg)

    def test_a_new_scf_drops_the_response(self):
        """Same atoms, same basis, new orbitals: the response must not survive."""
        msg = self.got["restale_refusal"]
        self.assertTrue(
            msg, "the exporter accepted a response built from orbitals that "
                 "the new SCF has since replaced")
        # scf_driver erases it, so the export finds nothing rather than judging
        # a fingerprint of a state it cannot identify.
        self.assertIn("absent", msg)

    def test_the_workflow_writes_four_readable_cubes(self):
        """scf_prop=acid has to produce files, not just avoid raising."""
        self.assertEqual(len(self.got["cube_names"]), 4)
        stems = [n.rsplit("_", 1)[-1] for n in self.got["cube_names"]]
        self.assertEqual(stems, ["acid.cube", "jx.cube", "jy.cube", "jz.cube"])
        self.assertTrue(all(self.got["cubes_exist"]), self.got["cube_names"])
        self.assertEqual(self.got["cube_natom"], 3)
        self.assertTrue(self.got["cube_npts_match"],
                        "the cube body does not fill the grid its header declares")
        self.assertTrue(self.got["cube_all_finite"])
        # ACID is sqrt(a sum of squares), so a negative value means the field
        # and the header disagree about layout.
        self.assertGreaterEqual(self.got["cube_acid_min"], 0.0)
        self.assertGreater(self.got["cube_acid_max"], 0.0)
        for size in self.got["cube_vector_sizes"]:
            self.assertEqual(size, self.got["cube_vector_sizes"][0])

    def test_a_different_system_invalidates_the_response(self):
        """Same coordinates and basis, different nuclei or charge."""
        charge = self.got["charge_refusal"]
        self.assertTrue(
            charge, "the exporter accepted a response built at a different "
                    "molecular charge")
        self.assertIn("charge", charge)
        z = self.got["z_refusal"]
        self.assertTrue(
            z, "the exporter accepted a response built for different nuclei, "
               "and would have written their atomic numbers into the cube "
               "header over another system's field")
        self.assertIn("Z=", z)

    def test_an_incomplete_response_is_refused(self):
        """The driver marks the buffer invalid until the response is finished."""
        msg = self.got["incomplete_refusal"]
        self.assertTrue(
            msg, "the exporter accepted a response the driver had marked "
                 "incomplete")
        self.assertIn("incomplete", msg)

    def test_a_permuted_ao_order_invalidates_the_response(self):
        """nbf and the density invariants cannot see a reordered basis."""
        self.assertTrue(self.got["permuted_rows_differ"],
                        "the permuted shells were identical, so the case was "
                        "never exercised")
        msg = self.got["permuted_refusal"]
        self.assertTrue(
            msg, "the exporter accepted a response indexed in a different AO "
                 "order than the basis it was about to be contracted with")
        self.assertIn("basis changed", msg)

    def test_primitives_sharing_a_summary_invalidate_the_response(self):
        """Exponents [1, 2] and [1.5, 1.5] share every per-shell summary."""
        self.assertTrue(self.got["prim_target_found"],
                        "no shell had two distinct primitives, so the case was "
                        "never exercised")
        self.assertEqual(self.got["prim_sum_drift"], 0.0,
                         "the substitution was supposed to preserve the "
                         "per-shell sum exactly")
        self.assertGreater(self.got["prim_values_moved"], 0.0)
        msg = self.got["prim_refusal"]
        self.assertTrue(
            msg, "the exporter accepted a response built on different AO "
                 "functions that happen to share a per-shell summary")
        self.assertIn("primitives changed", msg)

    def test_a_same_size_basis_swap_invalidates_the_response(self):
        """STO-6G has the same nbf as STO-3G, so only the shell block sees it."""
        self.assertEqual(self.got["rebasis_nbf"], self.got["nbf_small"],
                         "STO-6G was expected to have the same nbf as STO-3G, "
                         "so that nbf alone cannot catch the swap")
        msg = self.got["rebasis_refusal"]
        self.assertTrue(msg, "the exporter accepted a response built in a "
                             "different basis of the same size")
        # It has to be the basis block that answers.  No SCF ran, so the
        # response is still there to be judged; STO-6G differs from STO-3G in
        # contraction length, so which part of the block speaks first is an
        # implementation detail.
        self.assertIn("basis", msg)
        self.assertNotIn("absent", msg)

    def test_spherical_basis_is_refused_explicitly(self):
        message = self.got["spherical_refusal"]
        self.assertTrue(message,
                        "a spherical basis must be refused, not silently mapped")
        self.assertIn("Cartesian", message)


if __name__ == "__main__":
    unittest.main()
