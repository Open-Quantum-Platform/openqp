import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENERGY = ROOT / "source" / "modules" / "tdhf_mrsf_energy.F90"
LIB = ROOT / "source" / "tdhf_mrsf_lib.F90"
SINGLE_POINT = ROOT / "pyoqp" / "oqp" / "library" / "single_point.py"
OQPDATA = ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
INPUT_CHECKER = ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py"


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text.lower())


class UMRSFEnergyRegressionTests(unittest.TestCase):
    def test_umrsf_mixed_exchange_channels_follow_mrsf_permutation_pattern(self):
        source = compact(LIB.read_text())
        self.assertIn(
            "f3(:nf,9:10,i,k)=f3(:nf,9:10,i,k)-xval*d3(:nf,9:10,j,l)",
            source,
        )
        self.assertIn(
            "f3(:nf,9:10,k,i)=f3(:nf,9:10,k,i)-xval*d3(:nf,9:10,l,j)",
            source,
        )
        self.assertIn(
            "f3(:nf,9:10,i,l)=f3(:nf,9:10,i,l)-xval*d3(:nf,9:10,j,k)",
            source,
        )
        self.assertIn(
            "f3(:nf,9:10,l,i)=f3(:nf,9:10,l,i)-xval*d3(:nf,9:10,k,j)",
            source,
        )

    def test_umrsf_flag_is_scoped_to_umrsf_entry_point(self):
        source = compact(ENERGY.read_text())
        self.assertIn("subroutinetdhf_mrsf_energy_c", source)
        self.assertIn("inf%tddft%umrsf=.false.", source)
        self.assertIn("logical::previous_umrsf", source)
        self.assertIn("previous_umrsf=inf%tddft%umrsf", source)
        self.assertIn("inf%tddft%umrsf=previous_umrsf", source)

    def test_umrsf_jacobi_rotation_intent_and_diagonal_are_consistent(self):
        lib = LIB.read_text().lower()
        energy = compact(ENERGY.read_text())
        self.assertRegex(
            lib,
            r"real\(kind=dp\),\s*intent\(inout\),\s*dimension\(:,:\)\s*::\s*mo_a,\s*mo_b",
        )
        self.assertIn("mo_energy_work_a", energy)
        self.assertIn("mo_energy_work_a(i)=fa(i,i)", energy)
        self.assertIn("mo_energy_work_b(i)=fb(i,i)", energy)
        self.assertIn("callmrinivec(infos,mo_energy_work_a,mo_energy_work_b", energy)

    def test_umrsf_is_registered_but_gradient_and_optimizers_are_blocked_cleanly(self):
        oqpdata = compact(OQPDATA.read_text())
        checker = INPUT_CHECKER.read_text().lower()
        single = SINGLE_POINT.read_text().lower()

        self.assertIn("'umrsf'", oqpdata)
        self.assertIn("self._data.tddft.umrsf=td_type=='umrsf'", oqpdata)
        self.assertIn("td_type == \"umrsf\" and runtype == \"grad\"", checker)
        self.assertIn("td_type == \"umrsf\" and runtype in", checker)
        self.assertIn("umrsf-tddft gradients are not implemented", single)

    def test_umrsf_energy_does_not_use_mrsf_transition_density_output_path(self):
        source = compact(ENERGY.read_text())
        self.assertIn("if(umrsf)then", source)
        self.assertIn("trden=0.0_dp", source)
        self.assertIn("else", source)
        self.assertIn("callget_mrsf_transition_density", source)

    def test_spin_pair_scaling_avoids_hfscale_division_by_zero(self):
        source = compact(ENERGY.read_text())
        self.assertIn("if(abs(infos%tddft%hfscale)>epsilon(1.0_dp))then", source)
        self.assertIn("spc_scale_coco", source)
        self.assertIn("spc_scale_ovov", source)
        self.assertIn("spc_scale_coov", source)


if __name__ == "__main__":
    unittest.main()
