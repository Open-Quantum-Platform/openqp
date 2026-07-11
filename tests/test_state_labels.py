"""Tests for the physical-state labels used in user-facing logs."""

import importlib.util
from pathlib import Path


MODULE = Path(__file__).parents[1] / "pyoqp" / "oqp" / "utils" / "state_labels.py"
SPEC = importlib.util.spec_from_file_location("oqp_state_labels", MODULE)
state_labels = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(state_labels)


def mrsf_config(runtype="energy", target_mult=1):
    return {
        "input": {
            "method": "tdhf",
            "runtype": runtype,
            "functional": "bhhlyp",
            "basis": "6-31g*",
        },
        "scf": {"type": "rohf", "multiplicity": 3},
        "tdhf": {"type": "mrsf", "multiplicity": target_mult, "nstate": 5},
        "properties": {"grad": [1]},
        "optimize": {"istate": 1, "jstate": 2, "kstate": 3,
                     "imult": 1, "jmult": 3},
        "hess": {"state": 1},
        "nac": {"states": [[1, 2]]},
        "md": {"active": 1, "init_state": ""},
    }


def test_singlet_mrsf_root_one_is_s0_and_reference_is_not_s0():
    config = mrsf_config()

    assert state_labels.public_state_label(config, 0) == \
        "triplet ROHF reference (internal)"
    assert state_labels.public_state_label(config, 1) == "S0"
    assert state_labels.public_state_label(config, 2) == "S1"


def test_each_mrsf_spin_manifold_is_zero_based():
    triplet = mrsf_config(target_mult=3)
    quintet = mrsf_config(target_mult=5)

    assert state_labels.public_state_label(triplet, 1) == "T0"
    assert state_labels.public_state_label(triplet, 2) == "T1"
    assert state_labels.public_state_label(quintet, 2) == "Q1"


def test_opt_and_meci_summaries_use_physical_labels():
    opt = mrsf_config(runtype="optimize")
    meci = mrsf_config(runtype="meci")

    assert state_labels.requested_states(opt) == "S0"
    assert state_labels.requested_states(meci) == "S0/S1"


def test_mecp_can_label_two_spin_manifolds():
    config = mrsf_config(runtype="mecp")
    config["optimize"]["jstate"] = 1

    assert state_labels.requested_states(config) == "S0/T0"


def test_soc_summary_reports_counts_for_both_zero_based_manifolds():
    config = mrsf_config(runtype="soc")
    config["tdhf"].update(nstate=3, nstate_s=2, nstate_t=4)

    assert state_labels.requested_states(config) == "S0-S1 and T0-T3"
    text = state_labels.format_calculation_request(config)
    assert "Target spin:" not in text


def test_start_summary_separates_target_from_internal_reference():
    config = mrsf_config(runtype="optimize")

    text = state_labels.format_calculation_request(
        config, source="water.oqp", resolved="water.resolved.oqp")

    assert "Method:                      MRSF-TDDFT" in text
    assert "Physical target state(s):    S0" in text
    assert "Target spin:                 singlet" in text
    assert "Reference:                   triplet ROHF (internal working reference)" in text
    assert "Resolved input:              water.resolved.oqp" in text


def test_non_mrsf_state_labels_remain_neutral():
    config = {
        "input": {"method": "hf", "runtype": "grad", "basis": "def2-svp"},
        "tdhf": {"type": "rpa"},
    }

    assert state_labels.public_state_label(config, 0) == "state 0"
    assert state_labels.public_method_name(config) == "HF"


def test_mrsf_ekt_summary_labels_the_physical_parent_state():
    config = mrsf_config(runtype="ekt")
    config["tdhf"]["target"] = 1

    assert state_labels.requested_states(config) == "S0"


def test_tdhf_method_name_uses_functional_presence():
    config = mrsf_config()
    assert state_labels.public_method_name(config) == "MRSF-TDDFT"
    config["input"]["functional"] = ""
    assert state_labels.public_method_name(config) == "MRSF-TDHF"


def test_prop_and_branching_plane_summaries_keep_physical_meaning():
    prop = mrsf_config(runtype="prop")
    prop["properties"]["grad"] = [2]
    assert state_labels.requested_states(prop) == "S1"

    bp = mrsf_config(runtype="nac")
    bp["nac"]["bp"] = True
    text = state_labels.format_calculation_request(bp)
    assert "Calculation:                 Branching-plane analysis" in text
