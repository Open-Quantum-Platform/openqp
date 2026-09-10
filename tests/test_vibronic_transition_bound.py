"""The vibronic spectrum must bound its transition count before evaluating it."""

import unittest
from unittest import mock

import numpy as np

import oqp.library.vibronic as vibronic


class VibronicTransitionBound(unittest.TestCase):
    def test_state_product_is_bounded_before_any_overlap(self):
        nmode = 10
        model = vibronic.HarmonicVibronicModel.create(
            [1000.0] * nmode,
            [950.0] * nmode,
            np.eye(nmode),
            np.zeros(nmode),
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="synthetic ten-mode bound test",
        )
        # Nine quanta over ten modes is 92,378 states on each side -- each set
        # is admissible under max_states, but the product is ~8.5e9 transitions.
        with mock.patch.object(
            vibronic.HarmonicOverlapEngine,
            "overlap",
            side_effect=AssertionError("an overlap was evaluated"),
        ):
            with self.assertRaisesRegex(ValueError, "max_transitions"):
                vibronic.harmonic_vibronic_spectrum(
                    model,
                    electronic_origin_cm1=20000.0,
                    origin_kind="zero_zero",
                    max_final_quanta=9,
                    temperature_kelvin=300.0,
                    max_initial_quanta=9,
                )

    def test_small_spectra_are_unaffected(self):
        model = vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [900.0],
            [[1.0]],
            [0.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="synthetic one-mode test",
        )
        result = vibronic.harmonic_vibronic_spectrum(
            model,
            electronic_origin_cm1=20000.0,
            origin_kind="zero_zero",
            max_final_quanta=4,
            normalization="sum",
            minimum_franck_condon_completeness=0.0,
        )
        self.assertEqual(len(result.lines), 5)

    def test_nonpositive_bound_is_rejected(self):
        model = vibronic.HarmonicVibronicModel.create(
            [1000.0], [900.0], [[1.0]], [0.0],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="synthetic one-mode test",
        )
        with self.assertRaisesRegex(ValueError, "max_transitions must be positive"):
            vibronic.harmonic_vibronic_spectrum(
                model, electronic_origin_cm1=20000.0, origin_kind="zero_zero",
                max_final_quanta=1, max_transitions=0,
            )


if __name__ == "__main__":
    unittest.main()
