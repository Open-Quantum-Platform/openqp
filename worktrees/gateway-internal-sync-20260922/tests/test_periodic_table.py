import unittest

from pyoqp.oqp.periodic_table import ELEMENTS_LONG_NAME, ELEMENTS_NAME, MASSES, SYMBOL_MAP


class TestPeriodicTableCoverage(unittest.TestCase):
    def test_all_official_elements_have_names_masses_and_symbol_map_entries(self):
        self.assertGreaterEqual(len(ELEMENTS_NAME), 119)
        self.assertGreaterEqual(len(ELEMENTS_LONG_NAME), 119)
        self.assertGreaterEqual(len(MASSES), 119)

        for atomic_number in range(1, 119):
            symbol = ELEMENTS_NAME[atomic_number]
            long_name = ELEMENTS_LONG_NAME[atomic_number]

            self.assertIsInstance(symbol, str, atomic_number)
            self.assertTrue(symbol.strip(), atomic_number)
            self.assertIsInstance(long_name, str, atomic_number)
            self.assertTrue(long_name.strip(), atomic_number)
            self.assertGreater(MASSES[atomic_number], 0.0, atomic_number)

            self.assertEqual(SYMBOL_MAP[atomic_number], atomic_number)
            self.assertEqual(SYMBOL_MAP[str(atomic_number)], atomic_number)
            self.assertEqual(SYMBOL_MAP[symbol.strip()], atomic_number)
            self.assertEqual(SYMBOL_MAP[symbol.strip().lower()], atomic_number)
            self.assertEqual(SYMBOL_MAP[long_name.upper()], atomic_number)
            self.assertEqual(SYMBOL_MAP[long_name.lower()], atomic_number)


if __name__ == "__main__":
    unittest.main()
