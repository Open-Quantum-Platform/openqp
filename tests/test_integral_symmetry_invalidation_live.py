"""Turning the petite reduction OFF must actually turn it off.

The staging path writes `OQP::sym_shell_map`, a frame-specific AO operator
(`sym_op_blocks` in the default input-frame route, or `sym_ao_target` and
`sym_ao_sign` in the standard frame), and sets `sym_petite_enable` to 1. Those
tags live in a store that outlives one step, and nothing on the abelian path or
in the Fortran layer removes them. So a reused molecule that ran once with
`use_integral_symmetry=true` and then runs again with it false would keep
applying the previous geometry's petite list and skeleton symmetrisation --
wrong integrals, silently, because the reduction is energy-identical when it IS
valid and nothing prints when it is not.

This is a LIVE test rather than a unit test on purpose. The unit tests assert
that `stage_integral_symmetry_maps()` clears when called; the thing that can
regress here is whether it gets CALLED at all, since `SinglePoint.reference()`
used to invoke it only when symmetry was already on -- exactly the branch a
"symmetry is now off" run does not take.
"""

import os
import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TAGS = (
    'OQP::sym_shell_map',
    'OQP::sym_ao_target',
    'OQP::sym_ao_sign',
    'OQP::sym_atom_weight',
    'OQP::sym_op_blocks',
    'OQP::sym_nonabelian',
    'OQP::sym_petite_enable',
)


def _runtime_available():
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
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g
functional=
method=hf

[guess]
type=huckel

[scf]
type=rhf
multiplicity=1
conv=1.0e-9

[symmetry]
enabled=true
use_integral_symmetry={use}
"""


def _has(mol, tag):
    try:
        mol.data[tag]
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class IntegralSymmetryInvalidationTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.workdir = tempfile.mkdtemp(prefix='oqp_sym_invalidation_')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.workdir, ignore_errors=True)

    def _run(self, name, use):
        from oqp.pyoqp import Runner
        case = os.path.join(self.workdir, name)
        os.makedirs(case, exist_ok=True)
        path = os.path.join(case, f'{name}.inp')
        with open(path, 'w', encoding='utf-8') as handle:
            handle.write(INPUT.format(use=use))
        runner = Runner(project=name, input_file=path,
                        log=os.path.join(case, f'{name}.log'), usempi=False)
        runner.run()
        return runner.mol

    def test_the_reduction_stages_its_tags_when_enabled(self):
        """Control. Without this the OFF assertion below proves nothing --
        every tag would be absent whether or not invalidation works."""
        mol = self._run('sym_on', 'true')
        self.assertEqual(
            mol.symmetry_metadata.get('integral_symmetry', {}).get('status'),
            'active')
        for tag in ('OQP::sym_shell_map', 'OQP::sym_op_blocks',
                    'OQP::sym_nonabelian', 'OQP::sym_petite_enable'):
            self.assertTrue(_has(mol, tag), f'{tag} was not staged when enabled')

    def test_no_tags_survive_a_run_with_the_reduction_disabled(self):
        mol = self._run('sym_off', 'false')
        for tag in TAGS:
            self.assertFalse(_has(mol, tag),
                             f'{tag} present although the reduction is off')
        self.assertNotEqual(
            mol.symmetry_metadata.get('integral_symmetry', {}).get('status'),
            'active')

    def test_reusing_one_molecule_across_the_switch_clears_the_tags(self):
        """The actual regression: stage once, then re-run the SAME molecule
        with the flag off. `reference()` must still reach the invalidation."""
        from oqp.library.single_point import SinglePoint
        mol = self._run('sym_reuse', 'true')
        self.assertTrue(_has(mol, 'OQP::sym_petite_enable'))

        mol.symmetry_metadata['use_integral_symmetry'] = False
        mol.config['symmetry']['use_integral_symmetry'] = 'False'
        SinglePoint(mol).reference(do_init_scf=False)

        for tag in TAGS:
            self.assertFalse(_has(mol, tag),
                             f'{tag} survived into a run with the reduction off')


    def test_reload_config_route_still_clears(self):
        """The route the previous test MISSED.

        `load_config()` calls `initialize_symmetry_metadata()`, which replaces
        the metadata dict and erases `integral_symmetry.status`, while
        `apply_config()` leaves the native tags alone. A "was it staged
        before?" check keyed on metadata therefore answers False in exactly
        this case, and the staged maps survive with `sym_petite_enable=1`.
        The earlier test mutated the metadata in place and so never exercised
        it.
        """
        from oqp.library.single_point import SinglePoint
        mol = self._run('sym_reload', 'true')
        self.assertTrue(_has(mol, 'OQP::sym_petite_enable'))

        case = os.path.join(self.workdir, 'sym_reload')
        off = os.path.join(case, 'off.inp')
        with open(off, 'w', encoding='utf-8') as handle:
            handle.write(INPUT.format(use='false'))
        mol.load_config(off)
        # Precondition: the reload really did erase the status, so the fix
        # cannot be passing for the wrong reason.
        self.assertNotEqual(
            mol.symmetry_metadata.get('integral_symmetry', {}).get('status'),
            'active')
        self.assertTrue(mol.has_staged_integral_symmetry(),
                        'tags should still be live before reference() runs')

        SinglePoint(mol).reference(do_init_scf=False)
        for tag in TAGS:
            self.assertFalse(_has(mol, tag),
                             f'{tag} survived a load_config() that disabled the reduction')


    def test_tags_are_gone_before_the_initial_scf_runs(self):
        """Ordering, not just end state.

        With `init_scf != 'no'`, `_init_convergence()` runs a full
        two-electron SCF -- possibly in a DIFFERENT basis (`init_basis`) --
        before the staging call further down. Clearing only at staging time
        hands that initial SCF the previous run's maps with
        `sym_petite_enable=1`: a corrupted starting solution, and potentially a
        target SCF converged into a different basin.

        Asserting the final state cannot catch this, because the tags are gone
        by the end either way. So capture what was live at the moment
        `_init_convergence()` was entered.
        """
        from oqp.library.single_point import SinglePoint
        mol = self._run('sym_order', 'true')
        self.assertTrue(_has(mol, 'OQP::sym_petite_enable'))

        case = os.path.join(self.workdir, 'sym_order')
        off = os.path.join(case, 'off.inp')
        with open(off, 'w', encoding='utf-8') as handle:
            handle.write(INPUT.format(use='false'))
        mol.load_config(off)
        self.assertTrue(mol.has_staged_integral_symmetry(),
                        'precondition: tags live before reference()')

        sp = SinglePoint(mol)
        seen = {}
        original = sp._init_convergence

        def spy():
            seen['staged_at_entry'] = mol.has_staged_integral_symmetry()
            return original()

        sp._init_convergence = spy
        # The deck is RHF; 'no' would skip _init_convergence entirely and
        # the spy would never fire, which the assertion below catches.
        if sp.init_scf == 'no':
            sp.init_scf = 'rhf'
        sp.reference(do_init_scf=True)

        self.assertIn('staged_at_entry', seen,
                      '_init_convergence did not run; the test proves nothing')
        self.assertFalse(
            seen['staged_at_entry'],
            'stale petite maps were still live when the initial SCF started')


if __name__ == '__main__':
    unittest.main()
