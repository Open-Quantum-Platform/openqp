"""Static wiring checks for the rotated-axis direct pure-spherical ERI path."""
import os
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def _read(rel):
    with open(os.path.join(ROOT, rel), encoding='utf-8') as handle:
        return handle.read()


class RotaxisDirectPureWiring(unittest.TestCase):
    def test_pure_submodule_exists(self):
        src = _read('source/integrals/int_rotaxis_pure.F90')
        self.assertIn('submodule (int2e_rotaxis) int2e_rotaxis_pure', src)
        self.assertIn('module subroutine genr22_pure', src)
        # fused transform must come from the shared projection tables
        self.assertIn('int2_init_shell_projection', src)

    def test_parent_module_declares_interface(self):
        src = _read('source/integrals/int_rotaxis.F90')
        self.assertIn('public genr22_pure', src)
        self.assertIn('module subroutine genr22_pure', src)
        # the Cartesian path must keep its unrolled rotation
        self.assertIn('call r30s1d(jtype, grotspd, prot)', src)

    def test_driver_branches_to_direct_pure(self):
        src = _read('source/integrals/int2.F90')
        # both rotaxis call sites (shellquartet + ints_exchange) use the
        # direct-pure entry for quartets containing harmonic-flagged shells
        self.assertEqual(src.count('call genr22_pure('), 4)
        # the post-hoc projection stays only on the Cartesian branch
        self.assertEqual(src.count('call genr22_reduce_pure('), 2)


if __name__ == '__main__':
    unittest.main()
