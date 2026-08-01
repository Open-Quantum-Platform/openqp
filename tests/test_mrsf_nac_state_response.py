"""Focused regressions for the analytic MRSF NAC state response."""

import importlib.util
import os
from pathlib import Path
from types import SimpleNamespace
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "nac_analytic_state_response",
    ROOT / "pyoqp" / "oqp" / "library" / "nac_analytic.py",
)
assert SPEC is not None and SPEC.loader is not None
NAC_ANALYTIC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(NAC_ANALYTIC)


class _Data(dict):
    def __init__(self, bvec, cutoff):
        super().__init__({'OQP::td_bvec_mo': np.array(bvec, copy=True)})
        self._data = SimpleNamespace(
            control=SimpleNamespace(int2e_cutoff=cutoff)
        )


class _Mol:
    def __init__(self, bvec, cutoff):
        self.data = _Data(bvec, cutoff)


class MRSFNACStateResponseTests(unittest.TestCase):
    def test_orthonormal_response_and_gap_guard_are_resident_fortran(self):
        driver = (
            ROOT / "source" / "modules" / "mrsf_nac_driver.F90"
        ).read_text()
        self.assertIn(
            "gap = energies_saved(pair_j(batch_pair))", driver
        )
        self.assertIn("energies_saved(pair_i(batch_pair))", driver)
        self.assertIn("energies_saved = energies", driver)
        self.assertIn(
            "gap_floor = 128.0_dp*epsilon(1.0_dp)*energy_scale", driver
        )
        self.assertIn(".not. ieee_is_finite(gap)", driver)
        self.assertIn("abs(gap) <= gap_floor", driver)
        self.assertIn("bvec_saved(:,pair_i(batch_pair))/gap", driver)
        self.assertIn(
            "wpair_ytil(redundant_index,wpair_index) = 0.0_dp", driver
        )
        production = (
            ROOT / "pyoqp" / "oqp" / "library" / "nac_analytic.py"
        ).read_text()
        self.assertNotIn("_orthonormal_pair_response", production)
        self.assertNotIn("for I in range(nstate)", production)

    def test_mutated_state_is_restored_after_an_exception(self):
        names = NAC_ANALYTIC._NAC_MUTATED_ENV
        old_environment = {
            name: os.environ.get(name, NAC_ANALYTIC._MISSING)
            for name in names
        }

        def restore_test_environment():
            for name, value in old_environment.items():
                if value is NAC_ANALYTIC._MISSING:
                    os.environ.pop(name, None)
                else:
                    os.environ[name] = value

        self.addCleanup(restore_test_environment)
        os.environ['NAC_DUMP_ROHF_RESPONSE'] = 'previous-rohf-setting'

        original_bvec = np.arange(6.0).reshape(2, 3)
        mol = _Mol(original_bvec, cutoff=1.0e-12)

        @NAC_ANALYTIC._with_temporary_nac_state
        def fail_after_mutation(active_mol):
            active_mol.data['OQP::td_bvec_mo'] = np.full((2, 3), -9.0)
            active_mol.data._data.control.int2e_cutoff = 1.0e-20
            os.environ['NAC_DUMP_ROHF_RESPONSE'] = '1'
            raise LookupError('intentional regression-test failure')

        with self.assertRaisesRegex(LookupError, 'intentional'):
            fail_after_mutation(mol)

        np.testing.assert_array_equal(
            mol.data['OQP::td_bvec_mo'], original_bvec
        )
        self.assertEqual(mol.data._data.control.int2e_cutoff, 1.0e-12)
        self.assertEqual(
            os.environ['NAC_DUMP_ROHF_RESPONSE'], 'previous-rohf-setting'
        )


if __name__ == "__main__":
    unittest.main()
