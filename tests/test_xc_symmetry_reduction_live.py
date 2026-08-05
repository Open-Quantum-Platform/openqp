"""The XC symmetry reduction must be reachable from the path the SCF uses.

`symAtomWeight` lets run_xc integrate only one atom per symmetry orbit and
scale the surviving slices by the orbit size. It landed live in PR #184, and
PR #202 then added a density-driven XC branch that orphaned it:

  - the only setter sat in `dmatd_blk`, reachable only via `calc_dft_xc`;
  - `calc_dft_xc` has one call site, in the `else` of `if (use_density_xc)`;
  - `use_density_xc` is `present(dens_in)`, and every `calc_fock` call site in
    the tree passes a density.

So the branch was unreachable and the reduction never executed anywhere --
while the log still reported the feature as enabled. Benzene/6-31G*/BHHLYP,
D2h, 12 SCF iterations: XC build 0.311 s -> 0.108 s once revived, with the
total energy unchanged at -232.1002160904.

These tests pin the structure rather than the timing, so they stay meaningful
without a build.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SCF_ADDONS = ROOT / 'source/scf_addons.F90'
GRIDINT_ENERGY = ROOT / 'source/dftlib/dft_gridint_energy.F90'


def subroutine_body(path, name):
    """Text of `subroutine <name>` up to its `end subroutine`."""
    source = path.read_text()
    start = re.search(r'^\s*subroutine\s+%s\b' % name, source,
                      flags=re.IGNORECASE | re.MULTILINE)
    assert start, f'{name} not found in {path.name}'
    end = re.search(r'^\s*end\s+subroutine\s+%s\b' % name, source[start.end():],
                    flags=re.IGNORECASE | re.MULTILINE)
    assert end, f'end of {name} not found'
    return source[start.start():start.end() + end.end()]


class TheLivePathCanReachTheReduction(unittest.TestCase):
    def test_the_density_branch_gates_on_symmetry(self):
        """This branch is the one the SCF actually takes."""
        body = subroutine_body(SCF_ADDONS, 'calc_jk_xc')
        density_branch = body[body.index('if (use_density_xc) then'):]
        self.assertIn('get_sym_atom_weight', density_branch,
                      'the live XC path must consult the symmetry gate')

    def test_the_weight_reaches_the_grid_driver(self):
        body = subroutine_body(GRIDINT_ENERGY, 'dmatd_density_blk')
        self.assertIn('xc_opts%symAtomWeight => sym_atom_weight', body,
                      'dmatd_density_blk must set the weight run_xc consumes')

    def test_the_gate_is_shared_by_both_branches(self):
        """One helper, so the two branches cannot drift apart again."""
        source = SCF_ADDONS.read_text()
        self.assertEqual(source.count('subroutine get_sym_atom_weight'), 2,
                         'expected exactly one definition (and its end)')
        self.assertEqual(source.count('call get_sym_atom_weight'), 2,
                         'both XC branches must use the shared gate')


class ReducedGridsAreProjected(unittest.TestCase):
    def test_the_density_path_symmetrizes_its_skeleton(self):
        """Without this the SCF converges to a wrong energy.

        run_xc scales surviving slices by the orbit size, so E_xc and the
        electron count come out exact, but the XC MATRIX is a skeleton and
        must be projected onto the totally symmetric component.
        """
        body = subroutine_body(SCF_ADDONS, 'calc_dft_xc_density')
        reduced = body[body.index('if (present(sym_atom_weight)) then'):]
        cut = reduced.index('else')
        self.assertIn('symmetrize_skeleton_fock', reduced[:cut],
                      'the reduced branch must project its skeleton')

    def test_the_unreduced_branch_does_not_project(self):
        body = subroutine_body(SCF_ADDONS, 'calc_dft_xc_density')
        after_else = body[body.rindex('else'):]
        self.assertNotIn('symmetrize_skeleton_fock', after_else)


class PhiCacheCannotServeReducedWeights(unittest.TestCase):
    """The one way this change could produce a silently wrong number.

    run_xc bakes symw into the grid weights *before* they are stored, replay
    deliberately does not re-apply it, and the cache validity signature carries
    no symmetry key. The petite gate is toggled off mid-run at a fixed geometry
    by the stability/TRAH stage, so a slice cached under the reduction would be
    replayed unreduced.
    """

    def test_the_cache_is_disabled_whenever_a_weight_is_set(self):
        body = subroutine_body(GRIDINT_ENERGY, 'dmatd_density_blk')
        enable = body.index('xc_opts%use_phi_cache = (infos%control%xc_phi_cache')
        guard = body.index(
            'if (associated(xc_opts%symAtomWeight)) xc_opts%use_phi_cache = .false.')
        self.assertGreater(guard, enable,
                           'the guard must come after the enable, or it is '
                           'overwritten')

    def test_the_guard_precedes_the_grid_call(self):
        body = subroutine_body(GRIDINT_ENERGY, 'dmatd_density_blk')
        guard = body.index('xc_opts%use_phi_cache = .false.')
        run = body.index('call run_xc(')
        self.assertLess(guard, run)


class ResponseKernelsAreDeliberatelyNotPlumbed(unittest.TestCase):
    """Extending this to fxc/gxc/giao would be a correctness bug, not plumbing.

    Those kernels act on trial vectors that are not totally symmetric, so an
    orbit-reduced grid does not reproduce the full-grid result and there is no
    totally-symmetric projector to repair it.
    """

    def test_no_response_grid_path_sets_the_weight(self):
        for name in ('dft_gridint_fxc.F90', 'dft_gridint_gxc.F90',
                     'dft_gridint_giao.F90'):
            with self.subTest(file=name):
                text = (ROOT / 'source/dftlib' / name).read_text()
                self.assertNotIn('symAtomWeight', text)


if __name__ == '__main__':
    unittest.main()
