"""Source-level guards for native TRAH final-state handoff.

Native TRAH optimizes orbitals inside ``trah_native_run`` and the SCF driver
then reads orbitals, Fock matrices, and orbital energies through
``scf_conv_trah_result``. These guards make sure the native final state is
published through that result buffer before post-SCF consumers such as MRSF run.
"""

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCF_CONVERGER = ROOT / "source" / "scf_converger.F90"
TRAH_NATIVE = ROOT / "source" / "trah_converger.F90"


class NativeTRAHResultHandoffTests(unittest.TestCase):
    def test_trah_result_exposes_final_mo_energies(self):
        src = SCF_CONVERGER.read_text()

        result_type = re.search(
            r"type,\s*extends\(scf_conv_result\)\s*::\s*scf_conv_trah_result"
            r".*?end type scf_conv_trah_result",
            src,
            re.S | re.I,
        )
        self.assertIsNotNone(result_type, "missing scf_conv_trah_result type")
        block = result_type.group(0)

        self.assertIn("get_mo_e_a   => conv_result_trah_get_mo_e_a", block)
        self.assertIn("get_mo_e_b   => conv_result_trah_get_mo_e_b", block)
        self.assertIn("subroutine conv_result_trah_get_mo_e_a", src)
        self.assertIn("subroutine conv_result_trah_get_mo_e_b", src)
        self.assertIn("vector = self%dat%buffer(self%dat%slot)%mo_e_a", src)
        self.assertIn("vector = self%dat%buffer(self%dat%slot)%mo_e_b", src)

    def test_native_trah_publishes_final_orbitals_focks_density_and_energies(self):
        src = TRAH_NATIVE.read_text()

        run = re.search(
            r"subroutine trah_native_run\(.*?end subroutine trah_native_run",
            src,
            re.S | re.I,
        )
        self.assertIsNotNone(run, "missing trah_native_run")
        block = run.group(0)

        self.assertIn("call compute_native_mo_energies", block)
        self.assertIn("conv%dat%buffer(conv%dat%slot)%focks = conv%fock_ao", block)
        self.assertIn("conv%dat%buffer(conv%dat%slot)%densities = conv%dens", block)
        self.assertIn("conv%dat%buffer(conv%dat%slot)%energy = e0", block)
        self.assertIn("conv%dat%buffer(conv%dat%slot)%mo_a = conv%mo_a", block)
        self.assertIn("conv%dat%buffer(conv%dat%slot)%mo_e_a = mo_e_a", block)
        self.assertIn("conv%dat%buffer(conv%dat%slot)%mo_b = conv%mo_b", block)
        self.assertIn("conv%dat%buffer(conv%dat%slot)%mo_e_b = mo_e_b", block)
        self.assertNotIn("conv%dat%put", block)

        # The publication must happen after the final fresh Fock rebuild; otherwise
        # MRSF can receive the pre-TRAH reference while SCF reports the TRAH energy.
        self.assertLess(
            block.rindex("call build_fock_grad"),
            block.index("conv%dat%buffer(conv%dat%slot)%focks = conv%fock_ao"),
        )

    def test_native_mo_energy_helper_uses_final_fock_projected_into_final_mos(self):
        src = TRAH_NATIVE.read_text()
        helper = re.search(
            r"subroutine compute_native_mo_energies\(.*?end subroutine compute_native_mo_energies",
            src,
            re.S | re.I,
        )
        self.assertIsNotNone(helper, "missing native MO-energy helper")
        block = helper.group(0)

        self.assertIn("call unpack_matrix(fock, work_1)", block)
        self.assertIn("mo_coeffs", block)
        self.assertIn("mo_energies(i) = work_1(i, i)", block)


if __name__ == "__main__":
    unittest.main()
