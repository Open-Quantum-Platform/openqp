"""Policy tests for the production [dftb] model default resolution."""

from oqp.utils.input_checker import (
    apply_dftb_model_default,
    expand_dftb_method_alias,
)


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
    # MRSF-TDDFTB defaults to the tuned DTCAM-TB operator; the conventional
    # dftb+ protocol is opt-in (request it with model=dftb+).
    config = _config(dftb={"type": "mrsf"})
    assert apply_dftb_model_default(config) == "dtcam-tb"
    assert config["dftb"]["model"] == "dtcam-tb"


def test_mrsf_dftb_plus_is_opt_in():
    config = _config(dftb={"type": "mrsf", "model": "dftb+"})
    assert apply_dftb_model_default(config) == "dftb+"
    assert config["dftb"]["model"] == "dftb+"


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


# --- Top-level [input] method DFTB shortcuts (expand_dftb_method_alias) ---

def _bare(method):
    return {"input": {"method": method}, "dftb": {}, "tdhf": {}, "scf": {}}


def test_mrsf_tddftb_alias_expands_to_dftb_with_roks_triplet():
    config = _bare("mrsf-tddftb")
    assert expand_dftb_method_alias(config) == "mrsf"
    assert config["input"]["method"] == "dftb"
    assert config["dftb"]["type"] == "mrsf"
    assert config["dftb"]["backend"] == "native"
    assert config["tdhf"]["type"] == "mrsf"
    assert config["scf"]["type"] == "rohf"
    assert config["scf"]["multiplicity"] == 3


def test_sf_tddftb_alias_expands_to_sf_roks_triplet():
    config = _bare("sf-tddftb")
    assert expand_dftb_method_alias(config) == "sf"
    assert config["input"]["method"] == "dftb"
    assert config["dftb"]["type"] == "sf"
    assert config["tdhf"]["type"] == "sf"
    assert config["scf"]["type"] == "rohf"
    assert config["scf"]["multiplicity"] == 3


def test_tddftb_alias_expands_to_closed_shell():
    config = _bare("tddftb")
    assert expand_dftb_method_alias(config) == "tddftb"
    assert config["input"]["method"] == "dftb"
    assert config["dftb"]["type"] == "tddftb"
    assert config["scf"]["type"] == "rhf"
    assert config["scf"]["multiplicity"] == 1


def test_non_alias_method_is_left_alone():
    config = _bare("dftb")
    assert expand_dftb_method_alias(config) == ""
    assert config["input"]["method"] == "dftb"
    # no type forced onto a plain method=dftb input
    assert "type" not in config["dftb"]


def test_alias_expansion_then_model_default_is_dtcam_tb():
    # method=mrsf-tddftb -> method=dftb + type=mrsf, then the model default
    # materializes the tuned DTCAM-TB preset.
    config = _bare("mrsf-tddftb")
    config["input"]["runtype"] = "energy"
    config["dftb"]["model"] = ""
    expand_dftb_method_alias(config)
    assert apply_dftb_model_default(config) == "dtcam-tb"
