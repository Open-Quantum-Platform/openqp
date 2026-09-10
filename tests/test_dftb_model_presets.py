"""Input-checker coverage for the named [dftb] operator model presets.

``DFTB_MODELS`` is the gate every input passes before openqp-dftb sees it: a
name missing here is rejected with "Unknown OpenQP-DFTB operator model preset"
even when the native library implements it.  These tests drive the real
``[dftb]`` section check, so they fail if a preset (or one of its accepted
spellings) is dropped from the gate, if the settings log advertises a name the
gate rejects, or if a preset stops locking the operator keys it overrides.

The native library itself is not needed: the model check runs before any
openqp-dftb call.
"""

import pytest

from oqp.utils.input_checker import CheckReport, DFTB_MODELS, _check_tb
from oqp.utils.state_labels import DFTB_KNOWN_PRESETS


def _config(dftb=None):
    config = {
        "input": {"method": "dftb", "runtype": "energy"},
        "tdhf": {"type": "mrsf", "nstate": 3},
        "properties": {"grad": []},
        "optimize": {},
        "dftb": {"backend": "native", "type": "mrsf", "model": ""},
    }
    if dftb:
        config["dftb"].update(dftb)
    return config


def _model_errors(dftb):
    report = CheckReport()
    _check_tb(_config(dftb), report, section="dftb")
    return [d for d in report.errors if d.path == "dftb.model"]


@pytest.mark.parametrize("name", ["dtcam-gap", "dtcam_gap", "dtcamgap", "DTCAM-GAP"])
def test_dtcam_gap_spellings_pass_the_model_gate(name):
    # openqp-dftb lower-cases the name and accepts exactly these three spellings.
    assert _model_errors({"model": name}) == []


def test_unknown_model_is_still_rejected():
    # Control: proves the check above is actually reached, so an accepted name
    # is a real pass rather than a skipped check.
    errors = _model_errors({"model": "dtcam-gapx"})
    assert errors
    assert "Unknown OpenQP-DFTB operator model preset" in errors[0].message


@pytest.mark.parametrize("name", DFTB_KNOWN_PRESETS)
def test_every_advertised_preset_is_accepted(name):
    # DFTB_KNOWN_PRESETS feeds the "Available presets" settings-log line; a name
    # listed there must never be one the input gate rejects.
    assert name in DFTB_MODELS
    assert _model_errors({"model": name}) == []


def test_dtcam_gap_locks_the_operator_keys_it_overrides():
    # dtcam-gap overrides the spin-pair couplings (among others); a user-tuned
    # value would be silently discarded inside openqp-dftb, so it is refused.
    errors = _model_errors({"model": "dtcam-gap", "spc_coco": 0.123456})
    assert any("fixes the operator" in d.message for d in errors)
