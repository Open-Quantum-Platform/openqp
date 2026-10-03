"""DKH2 examples: decoupling in the decontracted basis ([scf] scal_rel_decontract).

The four HBr examples under examples/DKH share one geometry and one DKH2 RHF
setup and differ only in the basis (contracted vs uncontracted 6-31G*) and in
the decoupling basis (scal_rel_decontract=1, the default, vs 0, the legacy
contracted-basis route).  Their committed references encode two properties of
the decontracted route that do not depend on the platform:

* with an already uncontracted basis the two routes must coincide (the
  contraction matrix is the identity up to primitive order);
* with a contracted basis the decontracted route must lie above the
  uncontracted-basis energy (the contracted functions span a subspace of the
  primitives for the same Hamiltonian), while the legacy route falls below it.
"""
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DKH_EXAMPLES = ROOT / "examples" / "DKH"

CONTRACTED = "HBr_RHF_DKH2_ENERGY"
CONTRACTED_LEGACY = "HBr_RHF_DKH2_LEGACY_ENERGY"
UNCONTRACTED = "HBr_RHF_DKH2_UNCONTRACTED_ENERGY"
UNCONTRACTED_LEGACY = "HBr_RHF_DKH2_UNCONTRACTED_LEGACY_ENERGY"
ALL = (CONTRACTED, CONTRACTED_LEGACY, UNCONTRACTED, UNCONTRACTED_LEGACY)


def _reference_energy(name):
    with open(DKH_EXAMPLES / f"{name}.json", encoding="utf-8") as handle:
        return float(json.load(handle)["energy"])


class TestDKHExamples(unittest.TestCase):

    def _run_example(self, name):
        from oqp.utils.oqp_tester import OQPTester
        from oqp.utils.mpi_utils import MPIManager
        tester = OQPTester(output_dir="openqp_dkh_test_tmp", omp_threads=2,
                           mpi_manager=MPIManager())
        tester.run(str(DKH_EXAMPLES / f"{name}.inp"))
        result = tester.results[0]
        self.assertEqual(result["status"], "PASSED",
                         f"{name} regression failed:\n{result['message']}")

    def test_example_files_present(self):
        for name in ALL:
            for suffix in (".inp", ".json"):
                self.assertTrue((DKH_EXAMPLES / f"{name}{suffix}").is_file(),
                                f"missing: {name}{suffix}")
        self.assertTrue((DKH_EXAMPLES / "x2c-tzvpall_uncontracted_HBr.json").is_file())

    def test_reference_uncontracted_routes_coincide(self):
        # C is the identity for an uncontracted basis: both routes must agree to
        # numerical noise (the gates inside dk_scalar hold to 1e-12).
        self.assertAlmostEqual(_reference_energy(UNCONTRACTED),
                               _reference_energy(UNCONTRACTED_LEGACY), delta=1e-8)

    def test_reference_contracted_route_is_variational(self):
        e_unc = _reference_energy(UNCONTRACTED)
        e_new = _reference_energy(CONTRACTED)
        e_old = _reference_energy(CONTRACTED_LEGACY)
        # the decontracted route lies above the uncontracted-basis energy by the
        # contraction error; the legacy route lies below it (non-variational)
        self.assertGreater(e_new, e_unc)
        self.assertLess(e_new - e_unc, 0.5)
        self.assertLess(e_old, e_unc)

    def test_contracted_default(self):
        self._run_example(CONTRACTED)

    def test_contracted_legacy(self):
        self._run_example(CONTRACTED_LEGACY)

    def test_uncontracted_default(self):
        self._run_example(UNCONTRACTED)

    def test_uncontracted_legacy(self):
        self._run_example(UNCONTRACTED_LEGACY)


if __name__ == "__main__":
    unittest.main()
