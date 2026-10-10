"""Tests for the optional native NAMD hop trace observer."""

import csv
from types import SimpleNamespace

import numpy as np

from tools.diagnostics.trace_namd_hop import build_trace_row, install_trace


def test_two_state_trace_uses_nan_for_three_state_only_columns():
    row = build_trace_row(
        istep=2,
        dt_fs=0.5,
        active=1,
        hopped=False,
        last_hop_random=0.375,
        coef=np.array([np.sqrt(0.75), 0.5j]),
        params=np.array([0.0, 0.0, 0.0, 0.25]),
        tdc=np.array([[0.0, -0.1], [0.1, 0.0]]),
        cmhp=np.array([[0.0, 0.2], [0.3, 0.0]]),
        time_fs=0.625,
    )

    assert row["t_fs"] == 0.625
    assert row["random"] == 0.375
    assert np.isclose(row["pop_1"], 0.75)
    assert np.isclose(row["pop_2"], 0.25)
    assert row["tdc_12_au"] == -0.1
    assert row["p_12"] == 0.2
    assert row["p_21"] == 0.3
    for name in ("pop_3", "tdc_13_au", "tdc_23_au", "p_31", "p_32"):
        assert np.isnan(row[name])


def test_variant_loggers_write_trace_through_common_observer(tmp_path):
    class FakeNAMD:
        def _prepare_hop_step(self, _istep):
            return True

        def _log_step(self, istep, r, hopped=False,
                      transition_energy_jump=np.nan):
            return None

    class FakeQMMM(FakeNAMD):
        def _log_qmmm(self, istep, epot, hopped=False,
                      transition_energy_jump=np.nan):
            return None

    class FakeSOC(FakeNAMD):
        def _log_soc(self, istep, e_pure, mult, state, w, hopped):
            return None

    class FakeMCH(FakeNAMD):
        def _log_mch(self, istep, e_pure, mult, state, hopped):
            return None

    class FakeSOCQMMM(FakeNAMD):
        def _log_soc_qmmm(self, istep, epot, mult, state, w, hopped):
            return None

    class FakeMCHQMMM(FakeNAMD):
        def _log_mch_qmmm(self, istep, epot, mult, state, hopped):
            return None

    module = SimpleNamespace(
        NAMD=FakeNAMD,
        NAMD_QMMM=FakeQMMM,
        NAMD_SOC=FakeSOC,
        NAMD_SOC_MCH=FakeMCH,
        NAMD_SOC_QMMM=FakeSOCQMMM,
        NAMD_SOC_MCH_QMMM=FakeMCHQMMM,
    )
    output = tmp_path / "variant-hop.csv"
    install_trace(output, skip_first_hop=False, namd_module=module)

    def driver(cls):
        instance = cls()
        instance.coef = np.array([1.0 + 0.0j, 0.0 + 0.0j])
        instance.dt_fs = 0.5
        instance.active = 1
        instance._last_hop_random = 0.25
        instance.dt_adaptive = True
        instance._t_fs = 0.125 * instance.active
        instance.mol = SimpleNamespace(data={
            "OQP::namd_params": np.array([0.0, 0.0, 0.0, 0.25]),
            "OQP::namd_results": np.arange(4.0),
            "OQP::td_states_overlap": np.eye(2),
        })
        instance._compute_tdc = lambda overlap: overlap - overlap.T
        return instance

    driver(FakeQMMM)._log_qmmm(1, -1.0, hopped=True)
    soc_driver = driver(FakeSOC)
    soc_driver.coef = np.array([0.0 + 0.0j, 0.0 + 0.0j, 1.0 + 0.0j])
    soc_driver.active = 3
    soc_driver._last_hop_probabilities = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.4, 0.5, 0.0],
    ])
    soc_driver._t_fs = 0.125
    soc_driver._log_soc(2, -1.0, 1, 1, 1.0, True)
    driver(FakeMCH)._log_mch(3, -1.0, 1, 1, True)
    driver(FakeSOCQMMM)._log_soc_qmmm(4, -1.0, 1, 1, 1.0, True)
    driver(FakeMCHQMMM)._log_mch_qmmm(5, -1.0, 1, 1, True)

    with output.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert [row["step"] for row in rows] == ["1", "2", "3", "4", "5"]
    assert all(row["hopped"] == "1" for row in rows)
    assert all(row["random"] == "0.25" for row in rows)
    assert all(row["t_fs"] == "0.125" for row in rows)
    assert rows[1]["p_31"] == "0.4"
    assert rows[1]["p_32"] == "0.5"
