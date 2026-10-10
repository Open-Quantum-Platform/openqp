import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = (
    ROOT
    / "tests"
    / "fixtures"
    / "higher_root_gradient"
    / "tddft_bhhlyp_631gstar_gradient_matrix.json"
)


class HigherRootGradientRegressionTests(unittest.TestCase):
    def test_tddft_fixture_documents_higher_root_gradient_fd_agreement(self):
        """Stored patched-branch FD diagnostic covers TDDFT roots 1-3."""
        data = json.loads(FIXTURE.read_text())

        self.assertEqual("6-31g*", data["basis"])
        self.assertEqual("bhhlyp", data["functional"])
        self.assertEqual(4, data["nstate"])
        self.assertEqual(9, len(data["cases"]))

        tolerance = data["tolerance_ha_per_bohr"]
        failing = []
        for case in data["cases"]:
            if case["max_abs_diff"] > tolerance or case["trah_analytic"]:
                failing.append(
                    f"{case['molecule']} root {case['root']}: "
                    f"max_abs_diff={case['max_abs_diff']:.9f}, "
                    f"rms_diff={case['rms_diff']:.9f}, "
                    f"cosine={case['cosine_similarity']:.9f}, "
                    f"{case['worst_component']} "
                    f"ratio={case['worst_component_ratio_fd_over_analytic']}"
                )

        self.assertEqual(
            [],
            failing,
            "TDDFT analytic gradients must agree with central finite differences "
            f"within {tolerance} and avoid TRAH; mismatches: " + "; ".join(failing),
        )

    def test_tdhf_gradient_uses_target_state_transition_vectors(self):
        """TDDFT gradient must not flatten all roots and silently reuse root-1 X/Y."""
        source = (ROOT / "source" / "modules" / "tdhf_gradient.F90").read_text()

        self.assertRegex(
            source,
            r"xpy\s*\(:\s*,\s*:\s*\).*xmy\s*\(:\s*,\s*:\s*\)",
            "OQP_td_xpy/OQP_td_xmy are stored as [lexc, nstates]; the gradient "
            "reader must keep them rank-2 so target_state can select the requested root.",
        )
        self.assertRegex(
            source,
            r"xpy\s*\(:\s*,\s*infos%tddft%target_state\s*\)",
            "TDDFT gradient must build transition-density terms from the requested "
            "target_state, not from the first flat column of OQP_td_xpy.",
        )
        self.assertRegex(
            source,
            r"xmy\s*\(:\s*,\s*infos%tddft%target_state\s*\)",
            "TDDFT gradient must build transition-density terms from the requested "
            "target_state, not from the first flat column of OQP_td_xmy.",
        )


if __name__ == "__main__":
    unittest.main()
