import importlib.util
import re
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENERGY = ROOT / "source" / "modules" / "tdhf_mrsf_energy.F90"
LIB = ROOT / "source" / "tdhf_mrsf_lib.F90"
SINGLE_POINT = ROOT / "pyoqp" / "oqp" / "library" / "single_point.py"
OQPDATA = ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
INPUT_CHECKER = ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py"

# Every UMRSF runtype other than "energy" eventually drives a gradient,
# Hessian, or Z-vector, none of which are implemented for UMRSF. ("thermo"
# is also gradient-driven but is not in the checker's recognized runtype set,
# so it is rejected earlier as an unknown runtype rather than by this guard.)
UMRSF_BLOCKED_RUNTYPES = (
    "grad", "prop", "data", "hess", "nac", "nacme",
    "optimize", "meci", "mecp", "mep", "ts", "irc", "neb",
)


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text.lower())


def _load_input_checker():
    """Import input_checker.py in isolation with a minimal MPI stub."""
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        size = 1
        use_mpi = False

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils

    spec = importlib.util.spec_from_file_location(
        "input_checker_umrsf_under_test", INPUT_CHECKER
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _umrsf_config(runtype):
    return {
        "input": {
            "runtype": runtype,
            "method": "tdhf",
            "basis": "6-31g",
            "system": "\nO 0.0 0.0 0.0\nH 0.0 0.0 0.95\nH 0.9 0.0 -0.3",
        },
        "guess": {},
        "scf": {"type": "uhf", "multiplicity": 3},
        "tdhf": {"type": "umrsf", "nstate": 2},
        "properties": {"grad": [1]},
        "optimize": {"lib": "geometric", "istate": 1, "jstate": 2},
        "nac": {"states": [[1, 2]]},
        "neb": {"product": "", "nimage": 3},
    }


def _umrsf_guard_errors(report):
    return [
        diag
        for diag in report.errors
        if diag.path == "tdhf.type" and "only supports runtype=energy" in diag.message.lower()
    ]


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

    def test_umrsf_is_registered_as_a_tdhf_type(self):
        oqpdata = compact(OQPDATA.read_text())
        single = SINGLE_POINT.read_text().lower()

        self.assertIn("'umrsf'", oqpdata)
        self.assertIn("self._data.tddft.umrsf=td_type=='umrsf'", oqpdata)
        self.assertIn("umrsf-tddft gradients are not implemented", single)

    def test_umrsf_energy_runtype_is_not_blocked(self):
        checker = _load_input_checker()
        report = checker.check_input_values(
            _umrsf_config("energy"), raise_error=False, emit=False
        )
        self.assertEqual(
            _umrsf_guard_errors(report),
            [],
            "UMRSF energy must not be rejected by the runtype guard:\n"
            + report.to_text(),
        )

    def test_umrsf_non_energy_runtypes_are_blocked_at_the_single_choke_point(self):
        checker = _load_input_checker()
        for runtype in UMRSF_BLOCKED_RUNTYPES:
            with self.subTest(runtype=runtype):
                report = checker.check_input_values(
                    _umrsf_config(runtype), raise_error=False, emit=False
                )
                guard_errors = _umrsf_guard_errors(report)
                self.assertEqual(
                    len(guard_errors),
                    1,
                    f"runtype={runtype} should raise exactly one UMRSF guard "
                    f"error, got {len(guard_errors)}:\n" + report.to_text(),
                )
                self.assertIn(runtype, guard_errors[0].value)

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
