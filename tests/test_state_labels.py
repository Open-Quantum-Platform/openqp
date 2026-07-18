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


def dftb_config(dftb_type="ground", *, td_type="tda", nstate=3, model=""):
    return {
        "input": {
            "method": "dftb",
            "runtype": "energy",
            # This is a legacy schema default and must never be advertised as
            # the basis of a minimal-basis DFTB calculation.
            "basis": "6-31g*",
            "functional": "",
        },
        "scf": {"type": "rhf", "multiplicity": 1},
        "tdhf": {"type": td_type, "multiplicity": 1, "nstate": nstate},
        "dftb": {
            "type": dftb_type,
            "backend": "auto",
            "parameter_path": "",
            "model": model,
        },
        "properties": {"grad": []},
        "optimize": {"istate": 0},
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


def test_soc_namd_active_is_labeled_as_a_spin_adiabatic_surface():
    config = mrsf_config(runtype="namd")
    config["md"].update(soc=True, active=5)

    assert state_labels.requested_states(config) == "spin-adiabatic surface 5"
    config["md"]["init_state"] = "T0"
    assert state_labels.requested_states(config) == \
        "T0 initialization; spin-adiabatic surface 5"


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


def test_dftb_aliases_share_one_canonical_method_display():
    expected = {
        "dftb": ("ground", "DFTB"),
        "noscc": ("ground_noscc", "DFTB0"),
        "tda": ("tddftb", "TD-DFTB (TDA)"),
        "sf-tddftb": ("sf", "SF-TDDFTB"),
        "mrsftddftb": ("mrsf", "MRSF-TDDFTB"),
    }
    for alias, (canonical, display) in expected.items():
        config = dftb_config(alias)
        assert state_labels.canonical_dftb_type(alias) == canonical
        assert state_labels.public_method_name(config) == display


def test_dftb_calculation_summary_hides_fake_gaussian_basis():
    config = dftb_config("tddftb", nstate=4)
    text = state_labels.format_calculation_request(config)

    assert "Method:                      TD-DFTB (TDA)" in text
    assert "Physical target state(s):    singlet S1-S4" in text
    assert "Reference:                   singlet restricted SCC-DFTB reference" in text
    assert "6-31g*" not in text
    assert "Basis:" not in text


def test_dftb_target_summary_reports_the_selected_gradient_state():
    config = dftb_config("mrsf", td_type="mrsf", nstate=4)
    config["input"]["runtype"] = "grad"
    config["properties"]["grad"] = [2]

    assert state_labels.dftb_target_description(config) == "S1"


def test_dftb_mrsf_labels_use_backend_target_multiplicity():
    config = dftb_config("mrsf", td_type="mrsf", nstate=3)
    config["dftb"]["target_multiplicity"] = 3
    # The legacy TDHF field may remain at its schema default; the native DFTB
    # adapter sends dftb.target_multiplicity to the backend.
    config["tdhf"]["multiplicity"] = 1

    assert state_labels.public_state_label(config, 1) == "T0"
    assert state_labels.requested_states(config) == "T0-T2"
    assert state_labels.dftb_target_description(config) == "T0-T2"
    text = state_labels.format_calculation_request(config)
    assert "Target spin:                 triplet" in text


def test_explicit_ground_dftb_ignores_stale_mrsf_response_defaults():
    config = dftb_config("ground", td_type="mrsf", nstate=4)

    assert not state_labels.is_mrsf(config)
    assert state_labels.public_state_label(config, 0) == "state 0"
    text = state_labels.format_calculation_request(config)
    assert "Method:                      DFTB" in text
    assert "Physical target state(s):    S0 ground state" in text
    assert "Target spin:" not in text
    assert "State labels:" not in text


def test_sf_reference_multiplicity_one_logs_effective_triplet_reference():
    config = dftb_config("sf", td_type="sf")
    config["dftb"]["reference_multiplicity"] = 1

    assert state_labels.dftb_reference_description(config).startswith("triplet")


def test_open_shell_ground_dftb_is_not_labeled_singlet():
    config = dftb_config("ground")
    config["dftb"]["reference_multiplicity"] = 2

    assert state_labels.dftb_reference_description(config).startswith("doublet")
    assert state_labels.dftb_target_description(config) == "D0 ground state"


def test_conventional_tddftb_log_omits_inactive_sf_and_cam_controls():
    config = dftb_config("tddftb", nstate=3)
    text = state_labels.format_dftb_settings(config, backend="native")

    assert "Response manifold:" in text
    assert "closed-shell singlet TDA" in text
    assert "Spin completion:" not in text
    assert "Gamma kernel / omega:" not in text
    assert "CAM alpha / beta:" not in text


def test_dftb_settings_are_grouped_and_report_resolved_paths():
    config = dftb_config("mrsf", td_type="mrsf", nstate=3)
    config["dftb"].update(
        scc_mixer="trah",
        scc_tolerance=1.0e-9,
        response_solver="davidson",
        zvector=True,
    )
    text = state_labels.format_dftb_settings(
        config,
        backend="native",
        parameter_path="/parameters/OB2W0PT3",
        library_path="/site-packages/openqp_dftb/libopenqp_dftb_c.dylib",
    )

    assert "   Calculation\n" in text
    assert "   SCC\n" in text
    assert "   Response\n" in text
    assert "   Gradient\n" in text
    assert "   Operator\n" in text
    assert "MRSF-TDDFTB" in text
    assert "/parameters/OB2W0PT3" in text
    assert "trah (requested)" in text
    assert "davidson (requested)" in text
    assert "Relaxed Z-vector:" in text


def test_zero_capability_mask_never_claims_new_native_features():
    config = dftb_config("mrsf", td_type="mrsf", nstate=3)
    config["dftb"]["state_to_state_spectrum"] = True

    text = state_labels.format_dftb_settings(
        config,
        backend="native",
        library_path="/wheel/libopenqp_dftb_c.dylib",
        abi_version=3,
        capabilities=0,
    )

    assert "C API ABI:" in text
    assert "none advertised" in text
    assert "structured progress tracing" in text
    assert "loaded library (C ABI 3) does not advertise all-pair state spectrum" in text
    assert "loaded library (C ABI 3) does not advertise final trust-region SCC recovery" in text


def test_advertised_capabilities_enable_native_feature_summary():
    config = dftb_config("mrsf", td_type="mrsf", nstate=3)
    text = state_labels.format_dftb_settings(
        config,
        backend="native",
        abi_version=3,
        capabilities=(
            state_labels.DFTB_CAP_STRUCTURED_TRACE
            | state_labels.DFTB_CAP_STATE_SPECTRUM
            | state_labels.DFTB_CAP_SCC_FINAL_TRUST
        ),
    )

    assert "structured trace, state spectrum, SCC final trust recovery" in text
    assert "Native progress trace:" in text
    assert "level 1" in text
    assert "yes (unrelaxed TDA/state interaction)" in text
    assert "final charge/spin trust-TRAH pass when eligible" in text


def test_probe_settings_do_not_advertise_native_trace_or_spectrum():
    config = dftb_config("mrsf", td_type="mrsf", nstate=3)
    text = state_labels.format_dftb_settings(
        config,
        backend="probe",
        executable="/tmp/openqp_dftb_probe",
    )

    assert "Native progress trace:" not in text
    assert "requested; unavailable with probe backend" in text
    assert "Failed-cycle recovery:" in text
    assert "unavailable with probe backend" in text


def test_dftb_preset_log_does_not_claim_schema_defaults_are_effective():
    config = dftb_config("mrsf", td_type="mrsf", model="dtcam-tb")
    text = state_labels.format_dftb_settings(config, backend="native")

    assert "dtcam-tb preset" in text
    assert "dtcam-tb (operator resolved by openqp-dftb backend)" in text
    assert "Mixer:" not in text
    assert "CAM alpha / beta:" not in text
