"""SinglePoint.excitation must drop the previous run's imaginary-root status
before the TD solve.

A reused Molecule keeps ``td_imaginary_states`` from an earlier mixed-root RPA
run.  If the next TD solve raises or ends in a Davidson breakdown, the stale
status would make ``Molecule.get_results()`` hide energies and publish
``td_imaginary`` for a run that has none.
"""
import types
import unittest


def _single_point_module():
    try:
        from oqp.library import single_point
        return single_point
    except ImportError:
        return None


SP = _single_point_module()


@unittest.skipUnless(SP is not None, "oqp.library.single_point not importable")
class TestImaginaryStatusReset(unittest.TestCase):
    def _excitation_with(self, tddft):
        mol = types.SimpleNamespace(
            td_imaginary_states=(1,),             # left over from a previous run
            symmetry_metadata=None,
            config={'tdhf': {'type': 'rpa'}, 'input': {}},
            mol_energy=types.SimpleNamespace(Davidson_converged=False),
        )
        sp = SP.SinglePoint.__new__(SP.SinglePoint)
        sp.method = 'tdhf'
        sp.exception = True
        sp.mol = mol
        sp.tddft = tddft
        return sp, mol

    def test_breakdown_clears_previous_imaginary_states(self):
        sp, mol = self._excitation_with(lambda: None)
        original = SP.dump_log
        SP.dump_log = lambda *a, **k: None
        try:
            with self.assertRaises(SP.TDnotConverged):
                sp.excitation([-78.0])
        finally:
            SP.dump_log = original
        self.assertEqual(mol.td_imaginary_states, ())

    def test_exception_in_solver_clears_previous_imaginary_states(self):
        def boom():
            raise RuntimeError("solver failed")
        sp, mol = self._excitation_with(boom)
        with self.assertRaises(RuntimeError):
            sp.excitation([-78.0])
        self.assertEqual(mol.td_imaginary_states, ())


if __name__ == "__main__":
    unittest.main()
