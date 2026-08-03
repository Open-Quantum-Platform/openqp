"""Tests for the optional native NAMD hop trace observer."""

import numpy as np

from tools.diagnostics.trace_namd_hop import build_trace_row


def test_two_state_trace_uses_nan_for_three_state_only_columns():
    row = build_trace_row(
        istep=2,
        dt_fs=0.5,
        active=1,
        hopped=False,
        last_hop_random=0.25,
        coef=np.array([np.sqrt(0.75), 0.5j]),
        params=np.array([0.0, 0.0, 0.0, 0.25]),
        tdc=np.array([[0.0, -0.1], [0.1, 0.0]]),
        cmhp=np.array([[0.0, 0.2], [0.3, 0.0]]),
    )

    assert np.isclose(row["pop_1"], 0.75)
    assert np.isclose(row["pop_2"], 0.25)
    assert row["tdc_12_au"] == -0.1
    for name in ("pop_3", "tdc_13_au", "tdc_23_au", "p_31", "p_32"):
        assert np.isnan(row[name])
