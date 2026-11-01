"""`[pt2] h0` / `contraction` left unset take the method's own value.

The schema used to default `h0=fock` and `contraction=uncontracted`.  Those
values reached `_caspt2_options` and the preflight checker indistinguishable
from user input, so `method=nevpt2` (Dyall H0) without an explicit `h0=dyall`
was refused as a contradiction -- after the reference SCF had run.  The tests
here build the configuration through the schema parser, as Molecule.get_config
does, because a bare dict never carries the schema defaults that caused it.
"""
import pytest


_H4 = "\nH 0 0 0\nH 0.03 0 0.74\nH 0 0.06 1.53\nH 0.04 0 2.31"


def _schema_config(method, pt2_extra):
    from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
    from oqp.utils.input_parser import OQPConfigParser

    parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA)
    parser.load_dict({"input": {"method": method, "system": _H4},
                      "pt2": {"reference": "casci", **pt2_extra}})
    return parser.validate()


def _errors(config, path):
    from oqp.utils.input_checker import CheckReport, _check_pt2

    report = CheckReport()
    _check_pt2(config, report)
    return [item for item in report.diagnostics
            if item.severity == "ERROR" and item.path == path]


@pytest.mark.parametrize("method,pt2_extra,expect", [
    ("nevpt2", {}, ("dyall", "none")),
    ("sc-nevpt2", {}, ("dyall", "strong")),
    ("scnevpt2", {}, ("dyall", "strong")),
    ("nevpt2", {"h0": "dyall"}, ("dyall", "none")),
    ("sc-nevpt2", {"h0": "nevpt2", "contraction": "sc"}, ("dyall", "strong")),
    # caspt2 and the QDPT family still resolve an unset h0 to Fock.
    ("caspt2", {}, ("fock", "none")),
    ("mrmp2", {}, ("fock", "none")),
    ("caspt2", {"h0": "dyall"}, ("dyall", "none")),
    ("caspt2", {"h0": "dyall", "contraction": "strong"}, ("dyall", "strong")),
])
def test_unset_h0_and_contraction_take_the_method_value(method, pt2_extra, expect):
    from oqp.library.caspt2_dyall import _caspt2_options

    config = _schema_config(method, pt2_extra)
    options = _caspt2_options(config)
    assert (options.h0, options.contraction) == expect
    assert not _errors(config, "pt2.h0")
    assert not _errors(config, "pt2.contraction")


@pytest.mark.parametrize("method,pt2_extra,path", [
    ("nevpt2", {"h0": "fock"}, "pt2.h0"),
    ("nevpt2", {"h0": "caspt2"}, "pt2.h0"),
    ("sc-nevpt2", {"h0": "fock"}, "pt2.h0"),
    ("nevpt2", {"contraction": "strong"}, "pt2.contraction"),
    ("sc-nevpt2", {"contraction": "uncontracted"}, "pt2.contraction"),
    # caspt2 names Fock, so strong contraction with h0 unset is still refused.
    ("caspt2", {"contraction": "strong"}, "pt2.contraction"),
])
def test_explicit_contradiction_is_refused_by_preflight_and_options(
        method, pt2_extra, path):
    """An explicit value that contradicts the method is an error before the
    SCF (preflight) and again in the option reader."""
    from oqp.library.caspt2_dyall import _caspt2_options

    config = _schema_config(method, pt2_extra)
    assert _errors(config, path)
    with pytest.raises(ValueError, match=path.split(".")[1]):
        _caspt2_options(config)


@pytest.mark.parametrize("method,h0,contraction", [
    ("caspt2", "fock", "uncontracted"),
    ("mrmp2", "fock", "uncontracted"),
    ("nevpt2", "dyall", "uncontracted"),
    ("sc-nevpt2", "dyall", "strong"),
])
def test_effective_config_resolves_unset_pt2_by_method(method, h0, contraction):
    """The .oqp default stripping compares effective configurations, so an
    unset h0/contraction must resolve there exactly as the runtime resolves
    it: restating the method's own value is a default, a different value is
    not."""
    from oqp.utils import oqp_input

    schema = oqp_input._load_schema_defaults()
    unset = oqp_input._effective_config(
        {"input": {"method": method}, "pt2": {}}, schema)
    restated = oqp_input._effective_config(
        {"input": {"method": method},
         "pt2": {"h0": h0, "contraction": contraction}}, schema)
    assert (unset[("pt2", "h0")], unset[("pt2", "contraction")]) == (h0, contraction)
    assert restated == unset
    other = "fock" if h0 == "dyall" else "dyall"
    changed = oqp_input._effective_config(
        {"input": {"method": method}, "pt2": {"h0": other}}, schema)
    assert changed != unset
