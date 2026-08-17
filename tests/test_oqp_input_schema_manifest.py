"""Schema-coverage gate for the stdlib-only semantic ``.oqp`` parser."""

import ast
import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).parents[1]
SCHEMA_PATH = ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
INPUT_PATH = ROOT / "pyoqp" / "oqp" / "utils" / "oqp_input.py"
LEGACY_PARSER_PATH = ROOT / "pyoqp" / "oqp" / "utils" / "input_parser.py"


def _load_oqp_input():
    spec = importlib.util.spec_from_file_location(
        "_openqp_schema_manifest_input", INPUT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _schema_keys_from_ast():
    """Read section/key names without importing native ``oqp`` dependencies."""

    tree = ast.parse(SCHEMA_PATH.read_text(encoding="utf-8"))
    schema = next(
        node.value
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Name) and target.id == "OQP_CONFIG_SCHEMA"
            for target in node.targets
        )
    )
    return {
        ast.literal_eval(section): frozenset(ast.literal_eval(key) for key in values.keys)
        for section, values in zip(schema.keys, schema.values)
    }


def _schema_defaults_from_ast():
    tree = ast.parse(SCHEMA_PATH.read_text(encoding="utf-8"))
    schema = next(
        node.value
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Name) and target.id == "OQP_CONFIG_SCHEMA"
            for target in node.targets
        )
    )
    defaults = {}
    for section_node, options in zip(schema.keys, schema.values):
        section = ast.literal_eval(section_node)
        defaults[section] = {}
        for key_node, descriptor in zip(options.keys, options.values):
            key = ast.literal_eval(key_node)
            fields = {
                ast.literal_eval(name): value
                for name, value in zip(descriptor.keys, descriptor.values)
            }
            converter = fields["type"]
            converter_name = converter.id if isinstance(converter, ast.Name) else ""
            defaults[section][key] = (
                ast.literal_eval(fields["default"]), converter_name
            )
    return defaults


def _value_text(value, converter):
    if converter == "bool":
        enabled = value if isinstance(value, bool) else str(value).lower() in {
            "1", "true", "yes", "on"
        }
        return "true" if enabled else "false"
    if value is None:
        return "null"
    if isinstance(value, (int, float)):
        return str(value)
    return json.dumps(str(value))


def _generic_input(section, key, value_text):
    if section in {"droplet", "solute_com"}:
        return (
            'mrsf(nstate=2)/bhhlyp/6-31g* geom="h2o.xyz" namd(S1) '
            "%s(%s=%s)" % (section, key, value_text)
        )
    if section == "qmmm":
        return (
            'dft/pbe0/def2-svp geom="h2o.xyz" energy '
            "qmmm(%s=%s)" % (key, value_text)
        )
    if section == "mp2":
        return (
            'mp2/cc-pvdz geom="h2o.xyz" energy '
            "mp2(%s=%s)" % (key, value_text)
        )
    if section == "cc":
        # The [cc] solver controls only mean anything on a coupled-cluster
        # route, the same way [mp2] needs an mp2 route.
        return (
            'ccsd_t/cc-pvdz geom="h2o.xyz" energy '
            "cc(%s=%s)" % (key, value_text)
        )
    if section == "dftb":
        route = 'dftb geom="h2o.xyz" energy'
    elif section == "tdhf":
        route = 'mrsf/bhhlyp/6-31g* geom="h2o.xyz" energy'
    elif section in {"cas", "casscf", "ci", "state_average"}:
        # The active space, CI solver, orbital optimizer and state-average
        # weights are all owned by the multiconfigurational route.
        route = 'casscf/sto-3g geom="h2o.xyz" energy'
    elif section == "fci":
        route = 'fci/sto-3g geom="h2o.xyz" energy'
    elif section == "pt2":
        route = 'caspt2/sto-3g geom="h2o.xyz" energy'
    else:
        route = 'dft/pbe0/def2-svp geom="h2o.xyz" energy'
    return "%s %s(%s=%s)" % (route, section, key, value_text)


def test_trivial_crossing_following_is_opt_in():
    """Standard FSSH must not relabel the active root without an explicit opt-in."""

    assert _schema_defaults_from_ast()["md"]["trivial"] == ("False", "bool")


def test_namd_scientific_safety_defaults_are_minimal_input_defaults():
    defaults = _schema_defaults_from_ast()["md"]

    assert defaults["seed"] == ("0", "int")
    assert defaults["rng_stream"] == ("1", "int")
    assert defaults["first_hop_step"] == ("1", "int")
    assert defaults["thrshe"] == ("0.1", "float")
    assert defaults["nacme_check"] == ("baeck_an", "str")
    assert defaults["nve_gate"] == ("warn", "str")


def test_every_schema_keyword_has_exactly_one_semantic_input_owner():
    oqp_input = _load_oqp_input()
    schema = _schema_keys_from_ast()

    assert set(schema) == oqp_input.SECTION_NAMES
    assert schema == oqp_input.OQP_SCHEMA_KEYS

    owner_counts = Counter()
    for section, keys in schema.items():
        owners = oqp_input.SCHEMA_KEY_OWNERS[section]
        assert set(owners) == set(keys)
        owner_counts.update(owners.values())

    # 348 generic keys since the SA-CASSCF gradients added [casscf]
    # gradient_state, zvector_tol and zvector_degeneracy_tol.
    assert owner_counts == {
        "generic": 348,
        "route_driver": 142,
        "legacy_only": 20,
        "intentional_forbidden": 1,
    }
    assert oqp_input.INTENTIONALLY_FORBIDDEN_SCHEMA_KEYS == {
        "qmmm": frozenset({"istate"})
    }


def test_every_schema_section_except_xtb_is_fully_owned():
    """Enforce the ownership gate for everything xTB does not block.

    ``test_every_schema_keyword_has_exactly_one_semantic_input_owner`` is
    quarantined because ``[xtb]`` has no concise converter yet, and a whole-test
    xfail also hides regressions in the sections that *are* wired.  This keeps
    the rest of the schema enforced; delete it when the xTB entry is unblocked.
    """

    oqp_input = _load_oqp_input()
    schema = _schema_keys_from_ast()
    unowned = set(schema) - oqp_input.SECTION_NAMES

    assert unowned == {"xtb"}, unowned
    assert not oqp_input.SECTION_NAMES - set(schema)
    for section in sorted(set(schema) & oqp_input.SECTION_NAMES):
        assert schema[section] == oqp_input.OQP_SCHEMA_KEYS[section], section
        assert set(oqp_input.SCHEMA_KEY_OWNERS[section]) == set(schema[section])


def test_all_generic_schema_keys_survive_parse_render_reparse_and_lower():
    oqp_input = _load_oqp_input()
    defaults = _schema_defaults_from_ast()
    concise_dftb_defaults = {("dftb", "backend"), ("dftb", "parameter_path")}

    checked = []
    for section, keys in oqp_input.GENERIC_SCHEMA_KEYS.items():
        for key in sorted(keys):
            default, converter = defaults[section][key]
            if section == "tdhf" and key in {"nstate_s", "nstate_t"}:
                # The legacy schema's zero sentinel is not a valid explicit
                # semantic state count; canonical .oqp requires a positive
                # integer whenever either split-manifold count is supplied.
                default = 1
            text = _generic_input(
                section, key, _value_text(default, converter)
            )
            first = oqp_input.parse_canonical_oqp(text)
            canonical = oqp_input.render_canonical_oqp(first)
            reparsed = oqp_input.parse_canonical_oqp(canonical)
            lowered = oqp_input.lower_to_legacy(reparsed)
            if (section, key) in concise_dftb_defaults:
                assert key not in lowered.get(section, {})
            else:
                assert key in lowered.get(section, {}), (
                    section, key, canonical, lowered
                )
            checked.append((section, key))

    # This includes the native multiconfigurational sections plus the DFTB,
    # coupled-cluster, D4, and SCF controls now present on main.
    assert len(checked) == 352


def test_geometric_backend_is_canonical_only_through_opt_driver_options():
    oqp_input = _load_oqp_input()

    assert oqp_input.LEGACY_ONLY_SCHEMA_KEYS == {
        "optimize": frozenset({"optimizer", "step_size", "step_tol", "mep_maxit"}),
        "geometric": frozenset({
            "coordsys", "trust", "tmax", "convergence_set", "prefix",
            "hessian", "irc_direction", "constraints_file", "enforce", "conmethod",
        }),
        "neb": frozenset({"k", "maxg", "avgg", "climb", "align", "optep"}),
    }
    spec = oqp_input.parse_canonical_oqp(
        'dft/pbe0/def2-svp opt(S0,lib=geometric,coordsys=dlc) geom="h2o.xyz"'
    )
    assert oqp_input.lower_to_legacy(spec)["optimize"]["lib"] == "geometric"
    try:
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp opt(S0) geometric(coordsys=tric) geom="h2o.xyz"'
        )
    except oqp_input.OQPInputError as exc:
        assert "traditional sectioned .inp" in str(exc)
    else:
        raise AssertionError("exact geometric(...) remains a legacy section call")


def test_native_optimization_and_meci_defaults_are_automatic():
    defaults = _schema_defaults_from_ast()

    assert defaults["optimize"]["lib"] == ("oqp", "str")
    assert defaults["oqp"]["coordsys"] == ("auto", "str")
    assert defaults["oqp"]["auto_recovery"] == ("True", "bool")
    assert defaults["optimize"]["meci_search"] == ("auto", "str")


def test_route_driver_manifest_matches_public_driver_coverage():
    oqp_input = _load_oqp_input()
    owners = oqp_input.ROUTE_DRIVER_SCHEMA_KEYS

    # gap_sigma and mecp_search are [optimize] schema keys emitted only by the
    # crossing drivers, so the manifest owns them while DRIVER_OPTIONS keeps
    # them out of the drivers that never read them.
    assert owners["optimize"] == (
        set(oqp_input._OPT_OPTIONS)
        | set(oqp_input._CROSSING_OPTIONS)
        | set(oqp_input._MECP_ONLY_OPTIONS)
        | {"lib", "istate", "jstate", "kstate", "states", "imult", "jmult"}
    )
    assert owners["neb"] == {"product", "nimage"}
    assert owners["oqp"] == set().union(*oqp_input.OQP_DRIVER_OPTIONS.values())
    assert owners["hess"] == (
        set(oqp_input.DRIVER_OPTIONS["hess"]) | {"state"}
    )
    assert owners["nac"] == (
        set(oqp_input.DRIVER_OPTIONS["nac"])
        | set(oqp_input.DRIVER_OPTIONS["nacme"])
        | {"bp", "states"}
    )
    assert owners["md"] == set(oqp_input.DRIVER_OPTIONS["namd"])
    assert owners["ekt"] == set(oqp_input.DRIVER_OPTIONS["ekt"])

    internally_lowered = {
        "input": {"charge", "basis", "functional", "method", "runtype", "system", "system2"},
        "dftb": {"type", "reference_multiplicity", "target_multiplicity"},
        "scf": {"type", "multiplicity"},
        "tdhf": {"type", "multiplicity", "nstate", "target"},
        "properties": {"grad"},
    }
    for section, keys in internally_lowered.items():
        assert owners[section] == keys


def test_legacy_tests_section_rejects_unknown_options_too():
    spec = importlib.util.spec_from_file_location(
        "_openqp_legacy_input_parser_for_manifest", LEGACY_PARSER_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    schema = {
        "input": {"system": {"type": str, "default": ""}},
        "tests": {"exception": {"type": bool, "default": "False"}},
    }
    parser = module.OQPConfigParser(schema=schema)
    parser.load_dict({"input": {"system": "H 0 0 0"}, "tests": {"exeption": "true"}})

    try:
        parser.validate()
    except ValueError as exc:
        assert str(exc) == "Unknown option: tests.exeption"
    else:
        raise AssertionError("tests.* typos must not bypass schema validation")

    valid = module.OQPConfigParser(schema=schema)
    valid.load_dict(
        {"input": {"system": "H 0 0 0"}, "tests": {"exception": "true"}}
    )
    assert valid.validate()["tests"]["exception"] is True
