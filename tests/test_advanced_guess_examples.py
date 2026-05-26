import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples" / "other"


class TestAdvancedGuessExamples(unittest.TestCase):
    def test_sad_and_sap_examples_are_available_to_openqp_test_runner(self):
        expected = {
            "h2o_rhf_3-21g_sad.inp",
            "h2o_rhf_3-21g_sad.json",
            "h2o_rhf_3-21g_sap.inp",
            "h2o_rhf_3-21g_sap.json",
        }

        missing = sorted(name for name in expected if not (EXAMPLES / name).is_file())

        self.assertEqual(missing, [])

    def test_sad_and_sap_examples_use_advanced_guess_keywords(self):
        sad_input = (EXAMPLES / "h2o_rhf_3-21g_sad.inp").read_text()
        sap_input = (EXAMPLES / "h2o_rhf_3-21g_sap.inp").read_text()

        self.assertIn("[guess]", sad_input)
        self.assertIn("type=sad", sad_input)
        self.assertIn("[guess]", sap_input)
        self.assertIn("type=sap", sap_input)

    def test_cli_keeps_short_test_alias_for_example_runner(self):
        pyoqp_source = (ROOT / "pyoqp" / "oqp" / "pyoqp.py").read_text()

        self.assertIn("--test", pyoqp_source)
        self.assertIn("--run_tests", pyoqp_source)


if __name__ == "__main__":
    unittest.main()
