"""Common MD controls must retain explicit provenance across QM/MM lowering."""

from datetime import date
import importlib.util
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "pyoqp" / "oqp" / "utils" / "md_controls.py"
SPEC = importlib.util.spec_from_file_location("_md_controls_direct", MODULE_PATH)
md_controls = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = md_controls
SPEC.loader.exec_module(md_controls)


def test_schema_defaults_do_not_replace_explicit_qmmm_controls():
    qmmm = {
        "n_steps": 2,
        "trajectory_file": "alanine.pdb",
        "ensemble": "nvt",
        "temperature": 310.0,
    }
    materialized_md_defaults = {
        "common_controls": True,
        "common_control_keys": "",
        "nstep": 100,
        "trajectory_file": "",
        "thermostat": "off",
        "thermostat_temperature": 300.0,
    }
    assert md_controls.merge_explicit_md_controls(
        qmmm, materialized_md_defaults) == qmmm


def test_only_named_md_controls_override_qmmm_controls():
    merged = md_controls.merge_explicit_md_controls(
        {
            "n_steps": 2,
            "trajectory_file": "legacy.pdb",
            "ensemble": "nvt",
            "temperature": 310.0,
        },
        {
            "common_control_keys": "nstep,trajectory_file",
            "nstep": 8,
            "trajectory_file": "public.pdb",
            "thermostat": "off",
            "thermostat_temperature": 300.0,
        },
    )
    assert merged == {
        "n_steps": 8,
        "trajectory_file": "public.pdb",
        "ensemble": "nvt",
        "temperature": 310.0,
    }


def test_explicit_temperature_friction_and_thermostat_are_translated():
    merged = md_controls.merge_explicit_md_controls(
        {"ensemble": "nve", "temperature": 280.0, "friction": 0.5},
        {
            "common_control_keys": "temperature friction thermostat",
            "thermostat_temperature": 325.0,
            "thermostat_friction": 2.0,
            "thermostat": "langevin",
        },
    )
    assert merged == {
        "ensemble": "nvt",
        "temperature": 325.0,
        "friction": 2.0,
    }


def test_explicit_legacy_temperature_and_friction_are_translated():
    merged = md_controls.merge_explicit_md_controls(
        {"ensemble": "nve", "temperature": 280.0, "friction": 0.5},
        {
            "common_control_keys": (
                "ensemble init_temp thermostat_friction"
            ),
            "ensemble": "nvt",
            "init_temp": 315.0,
            "thermostat_friction": 1.5,
            "thermostat_temperature": 300.0,
        },
    )
    assert merged == {
        "ensemble": "nvt",
        "temperature": 280.0,
        "initial_temperature": 315.0,
        "friction": 1.5,
    }


def test_initial_and_thermostat_temperatures_remain_separate():
    merged = md_controls.merge_explicit_md_controls(
        {"ensemble": "nve", "temperature": 280.0, "friction": 0.5},
        {
            "common_control_keys": (
                "init_temp thermostat_temperature thermostat"
            ),
            "init_temp": 50.0,
            "thermostat_temperature": 300.0,
            "thermostat": "langevin",
        },
    )
    assert merged == {
        "ensemble": "nvt",
        "temperature": 300.0,
        "initial_temperature": 50.0,
        "friction": 0.5,
    }


def test_openmm_seed_is_reproducible_and_stream_separated():
    today = date(2026, 9, 29)
    first = md_controls.openmm_random_seed(0, 4, today=today)
    assert first == md_controls.openmm_random_seed(0, 4, today=today)
    assert first != md_controls.openmm_random_seed(0, 5, today=today)
    assert 1 <= first <= 2147483646
    with pytest.raises(ValueError, match="non-negative"):
        md_controls.openmm_random_seed(1, -1)


def test_qmmm_driver_consumes_velocity_and_seed_controls():
    source = (
        ROOT / "pyoqp" / "oqp" / "library" / "qmmm_md.py"
    ).read_text(encoding="utf-8")
    assert "merge_explicit_md_controls(qmmm_cfg, md_cfg)" in source
    assert "integrator.setRandomNumberSeed(seed)" in source
    assert "continuation_seed(seed, self._start_state[\"step\"])" in source
    assert "self._set_initial_velocities()" in source
    assert "context.setVelocitiesToTemperature(" in source
    assert "self.initial_temperature" in source
    assert "AU_VELOCITY_TO_NM_PER_PS" in source


def test_output_names_fall_back_when_the_materialised_default_is_empty():
    # The [qmmm] output-name defaults are empty strings, so a materialised
    # config carries the key with an empty value: get(key, fallback) returns
    # "" and the driver opened '' (FileNotFoundError at the first frame).
    materialized = {"trajectory_file": "", "log_file": "   ", "energy_file": ""}
    assert md_controls.resolve_output_name(
        materialized, "trajectory_file", "qmmm_trajectory.pdb"
    ) == "qmmm_trajectory.pdb"
    assert md_controls.resolve_output_name(
        materialized, "log_file", "qmmm_trajectory.dat"
    ) == "qmmm_trajectory.dat"
    assert md_controls.resolve_output_name(
        materialized, "energy_file", "total_energy.npz"
    ) == "total_energy.npz"

    # absent behaves like empty, and an explicit name always wins
    assert md_controls.resolve_output_name({}, "trajectory_file", "d.pdb") == "d.pdb"
    assert md_controls.resolve_output_name(
        {"trajectory_file": " ala-gsmd.pdb "}, "trajectory_file", "d.pdb"
    ) == "ala-gsmd.pdb"


def test_qmmm_md_schema_output_defaults_are_empty_so_the_fallback_matters():
    schema = (
        ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
    ).read_text(encoding="utf-8")
    for key in ("trajectory_file", "log_file", "energy_file"):
        assert "'%s': {'type': str, 'default': ''}" % key in schema

    source = (
        ROOT / "pyoqp" / "oqp" / "library" / "qmmm_md.py"
    ).read_text(encoding="utf-8")
    for key in ("trajectory_file", "log_file", "energy_file"):
        assert 'resolve_output_name(\n            qmmm_cfg, "%s"' % key in source
        assert 'qmmm_cfg.get("%s",' % key not in source


def test_explicit_ensemble_wins_over_legacy_thermostat_provenance():
    """job.workflow.md(thermostat="langevin") followed by md(ensemble="npt")
    leaves both keys in the provenance list; the thermostat branch then
    overwrote npt with nvt and the barostat silently never ran."""
    md = {"common_control_keys": "ensemble,thermostat",
          "ensemble": "npt", "thermostat": "langevin"}
    assert md_controls.merge_explicit_md_controls({}, md)["ensemble"] == "npt"
    md = {"common_control_keys": "ensemble,thermostat",
          "ensemble": "nve", "thermostat": "off"}
    assert md_controls.merge_explicit_md_controls({}, md)["ensemble"] == "nve"
    # the legacy thermostat alone still selects the ensemble
    md = {"common_control_keys": "thermostat", "thermostat": "langevin"}
    assert md_controls.merge_explicit_md_controls({}, md)["ensemble"] == "nvt"


def test_explicit_trajectory_interval_reaches_the_qmmm_reporter():
    """md(trajectory_interval=N) was accepted for QM/MM MD and then ignored:
    the driver's name for it is report_interval and nothing mapped the two."""
    md = {"common_control_keys": "trajectory_interval", "trajectory_interval": 10}
    assert md_controls.merge_explicit_md_controls({}, md)["report_interval"] == 10
    # only when named: the [md] schema default must not override [qmmm]
    md = {"common_control_keys": "nstep", "trajectory_interval": 1, "nstep": 5}
    merged = md_controls.merge_explicit_md_controls({"report_interval": 4}, md)
    assert merged["report_interval"] == 4


def test_sectioned_decks_supply_their_own_md_provenance(tmp_path):
    """A sectioned .inp deck has no common_control_keys marker, so a plain
    ``[md] nstep = 2`` looked like a schema default and QM/MM MD ignored it
    (running its own 1000 steps with Maxwell velocities)."""
    keys = md_controls.explicit_md_control_keys
    # the marker, when present, is authoritative (concise input, Python API)
    assert keys({"common_control_keys": "dt,nstep", "velocity": "zero"},
                written_keys={"velocity"}) == {"dt", "nstep"}
    # no marker: what the deck wrote is what the user set
    assert keys({"nstep": 2, "velocity": "zero"},
                written_keys={"nstep", "velocity", "common_controls"}) == {"nstep", "velocity"}
    # no marker and no deck information: nothing is assumed (materialised config)
    assert keys({"nstep": 100, "dt": 0.5}) == set()

    deck = tmp_path / "run.inp"
    deck.write_text("[input]\nruntype = md\n\n[md]\nnstep = 2\nVelocity = zero\n")
    assert md_controls.sectioned_md_keys(deck) == {"nstep", "velocity"}
    deck.write_text("[input]\nruntype = md\n")
    assert md_controls.sectioned_md_keys(deck) == set()
    concise = tmp_path / "run.oqp"
    concise.write_text('rks/bhhlyp/6-31g md(nstep=2)\ngeom="w.pdb"\n')
    assert md_controls.sectioned_md_keys(concise) is None
    assert md_controls.sectioned_md_keys(tmp_path / "missing.inp") is None
    assert md_controls.sectioned_md_keys(None) is None

    merged = md_controls.merge_explicit_md_controls(
        {"n_steps": 1000}, {"common_control_keys": "nstep", "nstep": 2})
    assert merged["n_steps"] == 2


def test_continuation_seed_gives_each_restart_point_its_own_stream():
    """Reseeding a restarted Langevin integrator with the user's seed replays
    the random sequence of the first segment in every segment."""
    f = md_controls.continuation_seed
    seeds = {f(7, step) for step in range(0, 4000, 100)}
    assert len(seeds) == 40                              # distinct per restart point
    assert 7 not in seeds                                # never the original stream
    assert f(7, 200) == f(7, 200)                        # reproducible
    assert f(7, 200) != f(8, 200)                        # still depends on the seed
    assert all(1 <= value <= 2147483646 for value in seeds)
