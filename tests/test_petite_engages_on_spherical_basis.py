"""The petite-list reduction must actually engage on a spherical basis.

This is the regression test for a bug that survived because it was invisible:
``stage_integral_symmetry_maps`` built its shell list without per-shell purity,
computed a Cartesian AO count, failed the ``n_ao != nbf`` guard and fell back to
the full C1 integral list -- on *every* spherical basis, which is cc-pVXZ,
aug-cc-pVXZ, def2 and the whole 6-311G family.

Two properties of that failure are worth stating, because they dictate the shape
of this test:

1. The fallback is a fail-safe, so the energy is **identical** either way. No
   assertion on energies, gradients or any physical quantity can detect it. The
   test must assert on the recorded status.
2. Cartesian and spherical shell sizes agree for l <= 1, so the AO total cannot
   distinguish "all shells pure" from OpenQP's real convention (only l >= 2 is
   spherical). A test that merely checks ``n_ao == nbf`` would pass on maps that
   apply the spherical component order to p shells and are silently wrong. So we
   also assert the per-shell purity flags directly.
"""

import os
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def _runtime_available():
    """True when the compiled OpenQP runtime can actually be driven.

    These are live calculations, so without this guard setUpClass raises
    during collection in a source-only or dependency-light checkout, turning
    an environment limitation into a test failure.
    """
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False

INPUT = """\
[input]
system=
   8   0.000000000   0.000000000   0.117300000
   1   0.000000000   0.757200000  -0.469200000
   1   0.000000000  -0.757200000  -0.469200000
charge=0
runtype=energy
basis={basis}
method=hf

[symmetry]
{symmetry}

[guess]
type=huckel

[scf]
type=rhf
multiplicity=1
conv=1.0e-9
"""


def _run(workdir, name, basis, symmetry):
    from oqp.pyoqp import Runner

    case = os.path.join(workdir, name)
    os.makedirs(case, exist_ok=True)
    path = os.path.join(case, f'{name}.inp')
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write(INPUT.format(basis=basis, symmetry=symmetry))
    runner = Runner(project=name, input_file=path,
                    log=os.path.join(case, f'{name}.log'), usempi=False)
    runner.run()
    return runner.mol


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class PetiteEngagesOnSphericalBasisTests(unittest.TestCase):
    """Water/C2v, cc-pVDZ (spherical d) -- the case that silently ran C1."""

    @classmethod
    def setUpClass(cls):
        cls.workdir = tempfile.mkdtemp(prefix='oqp_petite_spherical_')
        cls.c1 = _run(cls.workdir, 'water_c1', 'cc-pvdz', 'enabled=false')
        cls.sym = _run(cls.workdir, 'water_petite', 'cc-pvdz',
                       'enabled=true\nuse_integral_symmetry=true')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.workdir, ignore_errors=True)

    def test_reduction_is_active_not_silently_skipped(self):
        state = self.sym.symmetry_metadata.get('integral_symmetry', {})
        self.assertEqual(
            state.get('status'), 'active',
            'the petite reduction did not engage on a spherical basis; it fell '
            f'back to C1 with status {state.get("status")!r}. The energy is '
            'identical either way, so this assertion is the only thing that '
            'can catch it.')
        self.assertGreater(int(state.get('n_operations', 0)), 1)

    def test_energy_is_unchanged_by_the_reduction(self):
        e_c1 = float(np.asarray(self.c1.energies).ravel()[0])
        e_sym = float(np.asarray(self.sym.energies).ravel()[0])
        # Orbit weighting is algebraically exact for a totally symmetric
        # density; the residual is summation order only.
        self.assertLess(abs(e_sym - e_c1), 5.0e-9,
                        f'C1 {e_c1!r} vs petite {e_sym!r}')

    def test_shell_purity_follows_the_l_ge_2_rule(self):
        basis = self.sym.data.get_basis()
        spherical = basis.get('spherical')
        self.assertIsNotNone(
            spherical, 'per-shell purity is not exported; staging cannot know '
                       'the AO layout and will fall back to C1')
        angs = np.asarray(basis['angs'])
        spherical = np.asarray(spherical)
        # OpenQP transforms a shell only for l >= 2 (c2s_ncomp in
        # source/integrals/cart2sph.F90), so s and p stay Cartesian even here.
        self.assertTrue(np.all(spherical[angs < 2] == 0),
                        's or p shells reported as spherical')
        self.assertTrue(np.any(spherical[angs >= 2] == 1),
                        'cc-pVDZ must have spherical d shells')

    def test_cartesian_basis_reports_no_spherical_shells(self):
        """Guards the other direction: 6-31G* d shells are Cartesian."""
        mol = _run(self.workdir, 'water_cart', '6-31g*',
                   'enabled=true\nuse_integral_symmetry=true')
        basis = mol.data.get_basis()
        spherical = np.asarray(basis.get('spherical'))
        self.assertTrue(np.all(spherical == 0),
                        '6-31G* has Cartesian d shells; reporting them pure '
                        'would apply the wrong component order')
        self.assertEqual(
            mol.symmetry_metadata.get('integral_symmetry', {}).get('status'),
            'active')


if __name__ == '__main__':
    unittest.main()
