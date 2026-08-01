"""Regressions for resident ordered-pair MRSF NAC final assembly."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"
DRIVER = ROOT / "source" / "modules" / "mrsf_nac_driver.F90"
HEADER = ROOT / "include" / "oqp.h"
PRODUCTION = ROOT / "pyoqp" / "oqp" / "library" / "nac_analytic.py"
SPEC = importlib.util.spec_from_file_location("nac_analytic_pair_assembly", PRODUCTION)
assert SPEC is not None and SPEC.loader is not None
NAC_ANALYTIC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(NAC_ANALYTIC)


class MRSFNACPairAssemblyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SOURCE.read_text()
        cls.driver = DRIVER.read_text()
        cls.header = HEADER.read_text()
        cls.production = PRODUCTION.read_text()

    def test_resident_c_api_is_declared_and_used(self):
        for name in (
            "mrsf_nac_pair_accumulator_init",
            "mrsf_nac_pair_accumulate",
            "mrsf_nac_pair_finalize",
        ):
            self.assertIn(f'bind(C, name="{name}")', self.source)
            self.assertIn(f"void {name}(struct oqp_handle_t *inf", self.header)
        self.assertIn("call mrsf_nac_pair_accumulator_init(infos)", self.driver)
        self.assertIn("call mrsf_nac_pair_finalize(infos)", self.driver)
        self.assertIn("call mrsf_nac_pair_accumulate_antisym(", self.driver)
        self.assertIn("oqp.mrsf_nac_lagrangian(mol)", self.production)
        self.assertEqual(self.production.count("oqp.mrsf_nac_lagrangian(mol)"), 1)

    def test_fortran_owns_accumulation_antisymmetry_and_gap_scaling(self):
        accumulate = self.source.split(
            "subroutine mrsf_nac_pair_accumulate(infos", 1
        )[1].split("end subroutine mrsf_nac_pair_accumulate", 1)[0]
        finalize = self.source.split(
            "subroutine mrsf_nac_pair_finalize(infos)", 1
        )[1].split("end subroutine mrsf_nac_pair_finalize", 1)[0]
        self.assertIn("t1 = amp(coord,istate,jstate) + esum(cart,atom)", accumulate)
        self.assertIn("z_response = z_hf(cart,atom) + z_xc(cart,atom)", accumulate)
        self.assertIn("OQP_td_energies", finalize)
        self.assertIn("0.5_dp", finalize)
        self.assertIn("energies(jstate) - energies(istate)", finalize)
        self.assertIn("gap*dcv(coord,istate,jstate)", finalize)

        self.assertNotIn("dp = np.zeros", self.production)
        self.assertNotIn("dpa = 0.5", self.production)
        self.assertNotIn("nacv[I, J] = gap", self.production)

    def test_antisymmetric_accumulator_preserves_final_layout(self):
        accumulate = self.source.split(
            "subroutine mrsf_nac_pair_accumulate_antisym(infos", 1
        )[1].split("end subroutine mrsf_nac_pair_accumulate_antisym", 1)[0]
        self.assertIn(
            "value = nonz_antisym(coord) + z_hf(cart,atom) + z_xc(cart,atom)",
            accumulate,
        )
        self.assertIn("dp_ordered(coord,istate,jstate) = value", accumulate)
        self.assertIn("dp_ordered(coord,jstate,istate) = -value", accumulate)

    def test_fortran_pair_layout_is_exposed_without_numeric_reassembly(self):
        nstate = 3
        natom = 2
        ncoord = 3 * natom
        flat = np.empty(ncoord * nstate * nstate)
        for jstate in range(nstate):
            for istate in range(nstate):
                start = (istate + jstate * nstate) * ncoord
                flat[start:start + ncoord] = 100 * istate + 10 * jstate + np.arange(ncoord)

        exposed = NAC_ANALYTIC._resident_pair_cartesian(flat, nstate, natom)
        self.assertEqual(exposed.shape, (nstate, nstate, natom, 3))
        for istate in range(nstate):
            for jstate in range(nstate):
                expected = (
                    100 * istate + 10 * jstate + np.arange(ncoord)
                ).reshape(natom, 3)
                np.testing.assert_array_equal(exposed[istate, jstate], expected)

    def test_pair_layout_rejects_wrong_resident_size(self):
        with self.assertRaisesRegex(RuntimeError, "inconsistent size"):
            NAC_ANALYTIC._resident_pair_cartesian(np.zeros(7), 2, 1)


if __name__ == "__main__":
    unittest.main()
