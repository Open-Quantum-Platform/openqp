"""Validate the SG-1 pruned grid against its published definition.

SG-1 (P.M.W. Gill, B.G. Johnson, J.A. Pople, *Chem. Phys. Lett.* **209**
(1993) 506) is defined on a 50-point Murray-Handy-Laming / Euler-Maclaurin
radial grid, with the region boundaries alpha expressed in units of the
atomic radius.  Both halves of that statement have to hold for the grid to
*be* SG-1:

  * the radial size is 50, not whatever ``[dftgrid] rad_npts`` says;
  * the radial map is MHL, not whatever ``[dftgrid] rad_type`` says -- the
    same shell counts on a different map sit at different radii and would
    silently describe a different grid.

So SG-1 pins both, and the pin has to be **local to the grid being built**.
It is carried on the pruned-grid spec (``pruned%rad_grid_type``) and applied
inside ``dft_prepare_grid``; it must not be written back into ``infos``,
because ``infos`` is shared, long-lived user configuration and several code
paths build a temporary grid and then rebuild the production one (the MRSF
z-vector coarse-grid lever ``OQP_MRSF_ZV_COARSEGRID`` builds SG-1 by default
and then hands back to the production grid; the EKT-EA Fock rebuild
constructs another grid later still).  A write-back would leave every later
grid in the run on the MHL radial map and on the Gill Bragg-Slater radii
instead of the configured ones, with no diagnostic printed anywhere.

SG-1 is only defined for H-Ar.  Heavier atoms fall back to an unpruned
194-point sphere on the radial grid the user configured -- the pin above
covers H-Ar only.  Pinning heavy atoms to 50-point MHL as well cost
2.3e-4 Ha on a lone Kr atom and 1.7e-3 Ha on HBr against a converged grid.

These are behavioural checks: they build real grids through the compiled
runtime and compare energies and grid sizes.  Skipped unless the compiled
OpenQP runtime is importable.
"""

import contextlib
import io
import os
import re
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# rad_grid_type encoding, from OQPData._rad_grid_types / radial_grid_types.F90
RAD_MHL, RAD_TA, RAD_BECKE = 0, 2, 3

H2O = ("   8   0.000000000   0.000000000  -0.041061554\n"
       "   1  -0.533194329   0.533194329  -0.614469223\n"
       "   1   0.533194329  -0.533194329  -0.614469223")
# a lone atom above Ar: every atom is SG-1 type 4
KR = "  36   0.000000000   0.000000000   0.000000000"

INPUT_TMPL = """[input]
system=
{system}
charge=0
runtype=energy
basis=6-31g
functional=pbe
method=hf

[guess]
type=huckel
save_mol=false

[scf]
type=rhf
multiplicity=1
conv=1.0e-8

[dftgrid]
rad_npts={rad_npts}
ang_npts={ang_npts}
rad_type={rad_type}
pruned={pruned}
"""


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        import oqp.library  # noqa: F401
        from oqp.molecule import Molecule  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class SG1GridDefinition(unittest.TestCase):
    """SG-1 must be SG-1, and must not disturb anything built after it."""

    def _run(self, pruned, rad_npts=96, ang_npts=302, rad_type="ta",
             system=H2O):
        """Run one SCF; return (energy, grid points, rad_type before, after)."""
        import oqp
        import oqp.library
        from oqp.molecule import Molecule

        workdir = tempfile.mkdtemp()
        inp = os.path.join(workdir, "m.inp")
        log = os.path.join(workdir, "m.log")
        with open(inp, "w") as fh:
            fh.write(INPUT_TMPL.format(pruned=pruned, rad_npts=rad_npts,
                                       ang_npts=ang_npts, rad_type=rad_type,
                                       system=system))

        # banner/SCF chatter is noise here; only the log file is parsed
        with contextlib.redirect_stdout(io.StringIO()):
            mol = Molecule("m", inp, log)
            mol.load_config(inp)
            mol.data["OQP::log_filename"] = log
            oqp.oqp_banner(mol)
            oqp.library.set_basis(mol)
            before = int(mol.data._data.dft.rad_grid_type)
            oqp.library.ints_1e(mol)
            oqp.library.guess(mol)
            oqp.hf_energy(mol)
            after = int(mol.data._data.dft.rad_grid_type)
            energy = float(mol.data._data.mol_energy.energy)

        with open(log) as fh:
            sizes = re.findall(r"Molecular grid:\s+(\d+) points", fh.read())
        self.assertTrue(sizes, "no molecular grid size reported in the log")
        return energy, int(sizes[-1]), before, after

    def test_sg1_does_not_mutate_the_configured_radial_type(self):
        """The SG-1 pin is local: the user's rad_type survives the grid build.

        This is the regression guard.  If the override is written into
        ``infos%dft%rad_grid_type`` instead of the pruned spec, every grid
        built later in the same run silently switches to the MHL radial map
        and the Gill Bragg-Slater radii.
        """
        _, _, before, after = self._run("SG1", rad_type="ta")
        self.assertEqual(before, RAD_TA, "the deck should start on rad_type=ta")
        self.assertEqual(
            after, RAD_TA,
            "building an SG-1 grid overwrote the configured rad_type; the "
            "override must live on pruned%rad_grid_type, not in infos")

    def test_other_pruned_schemes_also_leave_the_radial_type_alone(self):
        """Control arm: SG-2 and the unpruned path have never pinned rad_type."""
        for pruned in ("SG2", "none"):
            with self.subTest(pruned=pruned):
                _, _, before, after = self._run(pruned, rad_type="ta")
                self.assertEqual(before, after)

    def test_sg1_ignores_the_user_radial_size_and_map(self):
        """SG-1 is pinned to its own definition, so these settings cannot matter.

        Two decks that disagree on *both* ``rad_npts`` and ``rad_type`` must
        produce the identical SG-1 grid.  This is what actually fails if the
        pin is dropped: unpinned, the first deck builds a 96-point TA grid
        and the second a 200-point Becke grid.
        """
        e1, n1, _, _ = self._run("SG1", rad_npts=96, rad_type="ta")
        e2, n2, _, _ = self._run("SG1", rad_npts=200, rad_type="becke")
        self.assertEqual(n1, n2,
                         "SG-1 grid size followed the user's rad_npts/rad_type")
        # Same grid, so the same energy to within summation order: the two
        # runs allocate different scratch sizes and reduce in a different
        # order, which is worth a few ULP.  Unpinned these decks differ in
        # the third decimal, so this tolerance is nowhere near the failure.
        self.assertAlmostEqual(
            e1, e2, places=10,
            msg="SG-1 energy followed the user's rad_npts/rad_type")

    def test_sg1_heavy_atoms_keep_the_configured_radial_grid(self):
        """Above Ar, SG-1 is undefined: keep the user's radial grid there.

        A lone Kr atom is all SG-1 type 4, so its SG-1 grid must be exactly
        the unpruned 194-point grid on the user's own rad_npts / rad_type --
        and must follow them.  Pinned to SG-1's 50-point MHL grid instead,
        SG-1 ignores both settings and this fails.
        """
        for rad_npts, rad_type in ((96, "ta"), (128, "becke")):
            with self.subTest(rad_npts=rad_npts, rad_type=rad_type):
                e_sg1, n_sg1, _, _ = self._run(
                    "SG1", rad_npts=rad_npts, rad_type=rad_type, system=KR)
                e_ref, n_ref, _, _ = self._run(
                    "none", rad_npts=rad_npts, ang_npts=194,
                    rad_type=rad_type, system=KR)
                self.assertEqual(
                    n_sg1, n_ref,
                    "SG-1 heavy atom is not on the configured radial grid")
                self.assertAlmostEqual(
                    e_sg1, e_ref, places=10,
                    msg="SG-1 heavy atom is not on the configured radial grid")

    def test_sg1_is_smaller_than_the_default_grid_it_is_pruned_from(self):
        """Sanity-check the pin is doing something, not silently inert.

        SG-1 on its own 50 radial shells is well under the shipped SG-2
        default.  An SG-1 that reported itself enabled while still running
        on ``rad_npts`` would fail *safe* -- same output shape, wrong grid
        -- so pin the order of magnitude here; the exact energy is pinned
        by the examples/ reference.
        """
        _, n_sg1, _, _ = self._run("SG1")
        _, n_sg2, _, _ = self._run("SG2")
        self.assertLess(n_sg1, 0.75 * n_sg2,
                        "SG-1 is not appreciably cheaper than SG-2; the "
                        "50-shell radial pin may not be taking effect")
        self.assertGreater(n_sg1, 0.10 * n_sg2)


if __name__ == "__main__":
    unittest.main()
