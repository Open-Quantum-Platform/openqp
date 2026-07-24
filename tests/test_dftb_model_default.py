"""Policy tests for the production [dftb] model default resolution."""

from oqp.utils.input_checker import apply_dftb_model_default


def _config(dftb=None, tdhf_type="mrsf", runtype="energy"):
    config = {
        "input": {"method": "dftb", "runtype": runtype},
        "tdhf": {"type": tdhf_type, "nstate": 3},
        "properties": {"grad": []},
        "optimize": {},
        "dftb": {"backend": "native", "type": "auto", "model": ""},
    }
    if dftb:
        config["dftb"].update(dftb)
    return config


def test_mrsf_defaults_to_dtcam_tb():
    config = _config(dftb={"type": "mrsf"})
    assert apply_dftb_model_default(config) == "dtcam-tb"
    assert config["dftb"]["model"] == "dtcam-tb"


def test_sf_defaults_to_dftb_plus():
    config = _config(dftb={"type": "sf"}, tdhf_type="sf")
    assert apply_dftb_model_default(config) == "dftb+"


def test_closed_shell_routes_stay_preset_free():
    # Singlet ground and closed-shell TD-DFTB: LC-DFTB2 is not implemented
    # for the restricted reference, so no preset default there.
    for dftb, td in ((({"type": "ground"}), "rpa"),
                     (({"type": "tddftb"}), "tda"),
                     (({"type": "ground_noscc"}), "rpa")):
        config = _config(dftb=dftb, tdhf_type=td)
        assert apply_dftb_model_default(config) == ""
        assert config["dftb"]["model"] == ""


def test_open_shell_ground_defaults_to_dftb_plus():
    config = _config(dftb={"type": "ground", "reference_multiplicity": 3},
                     tdhf_type="rpa")
    assert apply_dftb_model_default(config) == "dftb+"


def test_model_none_is_an_explicit_opt_out():
    config = _config(dftb={"type": "mrsf", "model": "none"})
    assert apply_dftb_model_default(config) == ""
    assert config["dftb"]["model"] == ""


def test_explicit_model_wins():
    config = _config(dftb={"type": "mrsf", "model": "dftb+"})
    assert apply_dftb_model_default(config) == "dftb+"
    assert config["dftb"]["model"] == "dftb+"


def test_tuned_locked_key_keeps_the_explicit_keys_route():
    # erf-tuning style inputs keep meaning what they said.
    config = _config(dftb={"type": "mrsf", "omega": 0.25, "lc_gamma": "erf"})
    assert apply_dftb_model_default(config) == ""
    assert config["dftb"]["model"] == ""


def test_probe_backend_gets_no_preset_default():
    config = _config(dftb={"type": "mrsf", "backend": "probe"})
    assert apply_dftb_model_default(config) == ""


def test_non_dftb_methods_are_untouched():
    config = _config(dftb={"type": "mrsf"})
    config["input"]["method"] = "tdhf"
    assert apply_dftb_model_default(config) == ""
    assert config["dftb"]["model"] == ""
