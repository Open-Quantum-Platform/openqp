"""MRSF property finite-difference thresholds must be finite.

Each safety gate is a comparison, and comparisons with NaN are false, so a NaN
threshold would otherwise pass validation and silently disable its gate.
"""

import math
import types
import unittest

import numpy as np

import oqp.library.mrsf_spectroscopy_fd as fd

REQUEST = dict(
    electronic_method="MRSF-TDDFT",
    electronic_state="MRSF S1",
    state_index=1,
    normal_modes=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
    displacement=1.0e-3,
    coordinate_phase_convention="synthetic fixed modes",
)


class PropertyThresholdsMustBeFinite(unittest.TestCase):
    def test_finite_request_is_accepted(self):
        fd.MRSFPropertyFDRequest.create(**REQUEST)

    def test_non_finite_request_thresholds_are_rejected(self):
        for name in (
            "minimum_tracking_margin",
            "fd_relative_tolerance",
            "fd_absolute_tolerance",
            "sos_tail_relative_tolerance",
            "sos_minimum_gap_hartree",
        ):
            for bad in (math.nan, math.inf):
                with self.subTest(name=name, value=bad):
                    with self.assertRaisesRegex(ValueError, "finite"):
                        fd.MRSFPropertyFDRequest.create(**REQUEST, **{name: bad})

    def test_truncated_sos_rejects_non_finite_controls(self):
        states = types.SimpleNamespace(energies=np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]))
        for controls in (
            {"tail_relative_tolerance": math.nan, "minimum_gap_hartree": 1.0e-5},
            {"tail_relative_tolerance": 0.05, "minimum_gap_hartree": math.nan},
        ):
            with self.subTest(**{k: repr(v) for k, v in controls.items()}):
                with self.assertRaisesRegex(ValueError, "finite"):
                    fd.truncated_sos_polarizability(states, 0, tail_states=2, **controls)


if __name__ == "__main__":
    unittest.main()
