"""Unit tests for the perf-preset resolver (utils/perf_levels). Pure Python; no liboqp.

The resolver fills the performance input keys (sentinel 'auto') in the parsed config from
the `perf` preset; an explicit value always wins. Env-var free.
"""
import os, importlib.util
_p = os.path.join(os.path.dirname(__file__), "..", "pyoqp", "oqp", "utils", "perf_levels.py")
_spec = importlib.util.spec_from_file_location("perf_levels", _p)
P = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(P)


def _auto_cfg():
    """A fully default-expanded config: every perf key present at the 'auto' sentinel
    (mirrors the OQPConfigParser-validated dict the scripting API passes)."""
    return {
        "input": {"perf": 1},
        "scf": {"xc_c2f": "auto", "xc_phi_cache": "auto", "xc_incdft": "auto",
                "grad_cutoff": "auto"},
        "tdhf": {"resp_cutoff": "auto", "fp32": "auto", "zv_warmstart": "auto"},
    }


def test_unset_changes_nothing():
    cfg = _auto_cfg()
    report = P.resolve(cfg, P.UNSET)
    assert report == []
    assert cfg["scf"]["xc_c2f"] == "auto" and cfg["tdhf"]["resp_cutoff"] == "auto"


def test_levels_fill_config_even_from_expanded_dict():
    # Regression for the codex review: a fully default-expanded config (all keys present
    # at 'auto') must still receive the preset -- 'auto' must NOT count as an override.
    want = {
        0: dict(xc_c2f="off", grad_cutoff="1.0d-10", resp_cutoff="5e-11", zv_warmstart="off", fp32="off"),
        1: dict(xc_c2f="off", grad_cutoff="1.0d-10", resp_cutoff="1e-8",  zv_warmstart="on",  fp32="off"),
        2: dict(xc_c2f="on",  grad_cutoff="1.0d-8",  resp_cutoff="1e-8",  zv_warmstart="on",  fp32="off"),
        3: dict(xc_c2f="on",  grad_cutoff="1.0d-7",  resp_cutoff="1e-6",  zv_warmstart="on",  fp32="off"),
    }
    for lvl, w in want.items():
        cfg = _auto_cfg()
        report = P.resolve(cfg, lvl)
        assert cfg["scf"]["xc_c2f"] == w["xc_c2f"], (lvl, cfg["scf"]["xc_c2f"])
        assert cfg["scf"]["grad_cutoff"] == w["grad_cutoff"]
        assert cfg["tdhf"]["resp_cutoff"] == w["resp_cutoff"]
        assert cfg["tdhf"]["zv_warmstart"] == w["zv_warmstart"]
        assert cfg["tdhf"]["fp32"] == w["fp32"]
        assert all(s == "perf" for _, _, s in report)     # every key came from the preset


def test_explicit_value_overrides_preset():
    cfg = _auto_cfg()
    cfg["scf"]["xc_c2f"] = "off"            # explicit, at perf=2 the preset wants 'on'
    cfg["tdhf"]["resp_cutoff"] = "5e-11"    # explicit
    report = P.resolve(cfg, 2)
    assert cfg["scf"]["xc_c2f"] == "off"
    assert cfg["tdhf"]["resp_cutoff"] == "5e-11"
    assert ("scf.xc_c2f", "off", "input") in report
    assert ("tdhf.resp_cutoff", "5e-11", "input") in report
    # a still-auto key is filled by the preset
    assert cfg["scf"]["grad_cutoff"] == "1.0d-8"


def test_invalid_perf_raises():
    try:
        P.resolve(_auto_cfg(), 7)
    except ValueError:
        return
    assert False, "perf=7 should raise"


def test_validator_flags_fp32_and_incdft():
    cfg = {"scf": {"xc_incdft": "on"}, "tdhf": {"fp32": "on"}}
    w = P.validate(cfg, scf_conv=1e-6, zv_conv=1e-10)
    joined = " ".join(w)
    assert "fp32" in joined and "incdft" in joined
    assert any("stall" in x for x in w)                  # zvconv 1e-10 < 1e-9


def test_validator_pscreen_cap_loose():
    cfg = {"scf": {"pscreen": True, "pscreen_cap": "1e-6"}}
    w = P.validate(cfg)
    assert any("pscreen_cap" in x for x in w)


def test_report_format():
    cfg = _auto_cfg(); cfg["scf"]["grad_cutoff"] = "1.0d-9"
    report = P.resolve(cfg, 1)
    block = P.format_report(1, report, [])
    assert "perf = 1" in block and "grad_cutoff" in block
    assert "(input)" in block and "(preset)" in block


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn(); print("ok", fn.__name__)
    print("ALL %d PASSED" % len(fns))
