"""Pure-Python tests for the markerless semantic .oqp input language."""

import configparser
import importlib.util
import sys
from pathlib import Path

import pytest


# Importing oqp normally loads liboqp.  The semantic parser intentionally has no
# native dependency, so exercise it directly to keep these tests fast and usable
# in source-only documentation/editor environments.
ROOT = Path(__file__).parents[1]
MODULE_PATH = ROOT / "pyoqp" / "oqp" / "utils" / "oqp_input.py"
SPEC = importlib.util.spec_from_file_location("_openqp_semantic_input", MODULE_PATH)
oqp_input = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = oqp_input
SPEC.loader.exec_module(oqp_input)

OQPInputError = oqp_input.OQPInputError


def _parse(text, source_dir=None):
    spec = oqp_input.parse_canonical_oqp(text)
    return spec, oqp_input.lower_to_legacy(spec, source_dir=source_dir)


def test_afqmc_route_model_is_a_ground_state_energy_method():
    spec, legacy = _parse(
        'afqmc/sto-3g geom="ch2.xyz" mult=3 energy afqmc(walkers=16,steps=20)'
    )
    assert legacy["input"]["method"] == "afqmc"
    assert legacy["scf"]["multiplicity"] == "3"
    assert legacy["afqmc"]["nwalkers"] == "16"
    assert "afqmc/sto-3g" in oqp_input.render_canonical_oqp(spec)


def test_afqmc_route_model_rejects_a_functional_and_non_energy_drivers():
    with pytest.raises(OQPInputError):
        _parse('afqmc/bhhlyp/sto-3g geom="ch2.xyz" energy')
    with pytest.raises(OQPInputError):
        _parse('afqmc/sto-3g geom="ch2.xyz" grad')


def test_afqmc_section_promotes_a_response_route_to_an_afqmc_run():
    # The route says which trial machinery to solve first; the afqmc(...)
    # section says the job is an AFQMC job. Without the promotion the method
    # would stay tdhf and [afqmc] would never be read.
    _, legacy = _parse(
        'mrsf(nstate=3)/bhhlyp/sto-3g geom="ch2.xyz" energy '
        'afqmc(trial=mrsf_state,state=2)'
    )
    assert legacy["input"]["method"] == "afqmc"
    assert legacy["tdhf"]["type"] == "mrsf"
    assert legacy["scf"]["type"] == "rohf"
    assert legacy["scf"]["multiplicity"] == "3"
    assert legacy["afqmc"]["trial"] == "mrsf_state"
    assert legacy["afqmc"]["state"] == "2"


def test_afqmc_section_promotes_a_reference_route_to_an_afqmc_run():
    _, legacy = _parse(
        'rohf/sto-3g geom="ch2.xyz" mult=3 energy afqmc(trial=mean_field)'
    )
    assert legacy["input"]["method"] == "afqmc"
    assert legacy["scf"]["type"] == "rohf"
    assert legacy["afqmc"]["trial"] == "mean_field"


def test_afqmc_section_call_and_concise_aliases_lower_to_legacy():
    spec, legacy = _parse(
        'mrsf(nstate=3)/bhhlyp/sto-3g geom="ch2.xyz" energy '
        'afqmc(trial=mrsf_state,state=2,walkers=32,steps=20,'
        'dt=0.01,chol=1e-9,accumulate_after=4)'
    )
    assert legacy["afqmc"]["trial"] == "mrsf_state"
    assert legacy["afqmc"]["state"] == "2"
    assert legacy["afqmc"]["nwalkers"] == "32"
    assert legacy["afqmc"]["nsteps"] == "20"
    assert legacy["afqmc"]["timestep"] == "0.01"
    assert legacy["afqmc"]["chol_tol"] == "1e-09"
    canonical = oqp_input.render_canonical_oqp(spec)
    assert "afqmc(trial=mrsf_state" in canonical


def test_mrsf_s1_opt_hides_reference_and_root_bookkeeping(tmp_path):
    spec, legacy = _parse(
        'mrsf(nstate=5)/bhhlyp/6-31g* geom="h2o.xyz" charge=0 '
        'opt(S1,maxit=100) scf(conv=1e-8)',
        tmp_path,
    )

    assert spec.physical_method == "MRSF-TDDFT"
    assert spec.reference_method == "ROHF"
    assert spec.reference_multiplicity == 3
    assert legacy["input"] == {
        "system": str((tmp_path / "h2o.xyz").resolve()),
        "charge": "0",
        "basis": "6-31g*",
        "functional": "bhhlyp",
        "method": "tdhf",
        "runtype": "optimize",
    }
    assert legacy["scf"]["type"] == "rohf"
    assert legacy["scf"]["multiplicity"] == "3"
    assert legacy["tdhf"]["multiplicity"] == "1"
    assert legacy["optimize"]["istate"] == "2"
    assert legacy["optimize"]["maxit"] == "100"
    assert legacy["scf"]["conv"] == "1e-08"


@pytest.mark.parametrize(
    "label,root,multiplicity",
    [("S0", "1", "1"), ("S2", "3", "1"), ("T0", "1", "3"), ("T1", "2", "3"), ("Q2", "3", "5")],
)
def test_mrsf_physical_state_mapping(label, root, multiplicity):
    _, legacy = _parse(
        'mrsf(nstate=5)/bhhlyp/6-31g* geom="h2o.xyz" grad(%s)' % label
    )
    assert legacy["properties"]["grad"] == root
    assert legacy["tdhf"]["multiplicity"] == multiplicity
    assert legacy["scf"]["multiplicity"] == "3"


def test_response_spin_is_selected_only_by_physical_state():
    _, legacy = _parse(
        'mrsf(nstate=4)/bhhlyp/6-31g* geom="h2o.xyz" energy(T0)'
    )
    assert legacy["scf"]["multiplicity"] == "3"
    assert legacy["tdhf"]["multiplicity"] == "3"
    assert "spin" not in legacy["tdhf"]

    _, conventional = _parse(
        'tddft(nstate=4)/pbe0/def2-svp geom="h2o.xyz" grad(T0)'
    )
    assert conventional["scf"]["multiplicity"] == "1"
    assert conventional["tdhf"]["multiplicity"] == "3"
    assert conventional["properties"]["grad"] == "1"

    with pytest.raises(OQPInputError, match="does not take spin"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=4,spin=triplet)/bhhlyp/6-31g* geom="h2o.xyz" grad(S1)'
        )


def test_mrsf_rejects_numeric_multiplicity_bookkeeping():
    with pytest.raises(OQPInputError, match="internal high-spin reference is automatic"):
        oqp_input.parse_canonical_oqp(
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" multiplicity=3 energy()'
        )
    with pytest.raises(OQPInputError, match="does not take multiplicity"):
        oqp_input.parse_canonical_oqp(
            'mrsf(multiplicity=1)/bhhlyp/6-31g* geom="h2o.xyz" energy()'
        )


def test_mult_is_canonical_and_long_spelling_is_only_an_input_alias():
    short = oqp_input.parse_canonical_oqp(
        'hf/6-31g* geom="radical.xyz" mult=2 energy()'
    )
    long = oqp_input.parse_canonical_oqp(
        'hf/6-31g* geom="radical.xyz" multiplicity=2 energy()'
    )
    assert short.options == long.options == {"geom": "radical.xyz", "mult": 2}
    assert "\nmult=2\n" in oqp_input.render_canonical_oqp(long)
    assert "multiplicity=" not in oqp_input.render_canonical_oqp(long)
    with pytest.raises(OQPInputError, match="not a route option"):
        oqp_input.parse_canonical_oqp(
            'mrsf(type=rpa)/bhhlyp/6-31g* geom="h2o.xyz" energy()'
        )


def test_route_options_are_small_and_section_options_are_not_silently_dropped():
    _, mp2 = _parse(
        'mp2(variant=scs-mp2,same_spin_scale=0.33)/cc-pvdz geom="h2o.xyz" energy()'
    )
    assert mp2["mp2"] == {"variant": "scs-mp2", "same_spin_scale": "0.33"}

    with pytest.raises(OQPInputError, match=r"tdhf\(nvdav=\.\.\.\)"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=5,nvdav=30)/bhhlyp/6-31g* geom="h2o.xyz" energy()'
        )
    with pytest.raises(OQPInputError, match="not a route option"):
        oqp_input.parse_canonical_oqp(
            'dft(nstate=5)/pbe0/def2-svp geom="h2o.xyz" energy()'
        )
    with pytest.raises(OQPInputError, match="route-owned"):
        oqp_input.parse_canonical_oqp(
            'mp2(variant=scs-mp2)/cc-pvdz geom="h2o.xyz" energy() '
            'mp2(variant=sos-mp2)'
        )


def test_mp2_reference_is_route_owned_and_lowers_to_scf_type():
    spec, uhf = _parse(
        'mp2(reference=uhf,variant=scs-mp2)/6-31g geom="h2o.xyz" '
        'mult=1 energy'
    )
    assert spec.reference_method == "UHF"
    assert uhf["scf"] == {"type": "uhf", "multiplicity": "1"}
    assert uhf["mp2"] == {"variant": "scs-mp2"}
    assert "reference" not in uhf["mp2"]
    assert "mp2(reference=uhf,variant=scs-mp2)/6-31g" in (
        oqp_input.render_canonical_oqp(spec)
    )

    _, rhf = _parse(
        'mp2(reference=rhf)/cc-pvdz geom="h2o.xyz" mult=1 energy'
    )
    _, rohf = _parse(
        'mp2(reference=rohf)/cc-pvdz geom="o2.xyz" mult=3 energy'
    )
    assert rhf["scf"] == {"type": "rhf", "multiplicity": "1"}
    assert rohf["scf"] == {"type": "rohf", "multiplicity": "3"}
    assert "mp2" not in rhf
    assert "mp2" not in rohf

    with pytest.raises(OQPInputError, match="must be rhf, rohf, or uhf"):
        oqp_input.parse_canonical_oqp(
            'mp2(reference=rks)/cc-pvdz geom="h2o.xyz" energy'
        )
    with pytest.raises(OQPInputError, match="reference=rhf requires mult=1"):
        oqp_input.parse_canonical_oqp(
            'mp2(reference=rhf)/cc-pvdz geom="o2.xyz" mult=3 energy'
        )
    with pytest.raises(OQPInputError, match="reference=rohf requires"):
        oqp_input.parse_canonical_oqp(
            'mp2(reference=rohf)/cc-pvdz geom="h2o.xyz" mult=1 energy'
        )


def test_umrsf_uses_triplet_uhf_reference():
    spec, legacy = _parse(
        'umrsf(nstate=3)/bhhlyp/6-31g* geom="radical.xyz" energy()'
    )
    assert spec.reference_method == "UHF"
    assert legacy["scf"] == {"type": "uhf", "multiplicity": "3"}
    assert legacy["tdhf"]["type"] == "umrsf"


def test_explicit_hf_reference_routes_are_not_silently_changed():
    _, rohf = _parse('rohf/6-31g* geom="o2.xyz" mult=3 energy()')
    _, uhf = _parse('uhf/6-31g* geom="h2.xyz" mult=1 energy()')
    assert rohf["scf"] == {"type": "rohf", "multiplicity": "3"}
    assert uhf["scf"] == {"type": "uhf", "multiplicity": "1"}

    with pytest.raises(OQPInputError, match="rhf requires mult=1"):
        oqp_input.parse_canonical_oqp('rhf/6-31g* geom="o2.xyz" mult=3 energy()')
    with pytest.raises(OQPInputError, match="rohf requires an open-shell mult"):
        oqp_input.parse_canonical_oqp('rohf/6-31g* geom="h2o.xyz" energy()')


def test_mrsf_tdhf_alias_has_a_basis_only_route():
    spec, legacy = _parse(
        'mrsf-tdhf(nstate=3)/6-31g* geom="h2o.xyz" opt(S0)'
    )
    assert spec.physical_method == "MRSF-TDHF"
    assert "functional" not in legacy["input"]
    assert legacy["tdhf"]["type"] == "mrsf"


@pytest.mark.parametrize(
    "route,td_type,scf_type",
    [
        ("sf-tdhf(nstate=3)/6-31g*", "sf", "rohf"),
        ("umrsf-tdhf(nstate=3)/6-31g*", "umrsf", "uhf"),
        ("cis(nstate=3)/6-31g*", "tda", "rhf"),
    ],
)
def test_hf_response_variants_have_basis_only_routes(route, td_type, scf_type):
    _, legacy = _parse('%s geom="h2o.xyz" energy()' % route)
    assert "functional" not in legacy["input"]
    assert legacy["tdhf"]["type"] == td_type
    assert legacy["scf"]["type"] == scf_type


@pytest.mark.parametrize(
    "route,expected",
    [
        ("rks/pbe0/def2-svp", "rhf"),
        ("uks/pbe0/def2-svp", "uhf"),
        ("roks/pbe0/def2-svp", "rohf"),
    ],
)
def test_explicit_dft_reference_routes(route, expected):
    mult = 1 if expected == "rhf" else 3
    _, legacy = _parse('%s geom="o2.xyz" mult=%d energy()' % (route, mult))
    assert legacy["scf"]["type"] == expected


def test_one_primary_driver_only_and_no_grad_opt_fallback():
    with pytest.raises(OQPInputError, match="Exactly one primary driver"):
        oqp_input.parse_canonical_oqp(
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" grad(S1) opt(S0)'
        )


def test_no_driver_defaults_to_implicit_energy_and_renderer_keeps_it_implicit():
    spec, legacy = _parse('dft/pbe0/def2-svp geom="h2o.xyz" charge=0')
    assert spec.driver.name == "energy"
    assert spec.driver.explicit is False
    assert legacy["input"]["runtype"] == "energy"
    rendered = oqp_input.render_canonical_oqp(spec)
    assert "\nenergy\n" not in rendered
    assert rendered == 'dft/pbe0/def2-svp\ngeom="h2o.xyz"\n'


def test_explicit_bare_energy_renders_as_the_same_implicit_default():
    implicit, _ = _parse('dft/pbe0/def2-svp geom="h2o.xyz"')
    explicit, _ = _parse('dft/pbe0/def2-svp geom="h2o.xyz" energy()')
    assert oqp_input.render_canonical_oqp(explicit) == (
        oqp_input.render_canonical_oqp(implicit)
    )


def test_energy_state_selector_is_not_dropped():
    spec, legacy = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" energy(T0)'
    )
    assert legacy["tdhf"]["multiplicity"] == "3"
    assert "\nenergy(T0)\n" in oqp_input.render_canonical_oqp(spec)


@pytest.mark.parametrize("driver", ["opt", "opt()"])
def test_mrsf_opt_without_state_defaults_to_s0(driver):
    _, legacy = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" %s' % driver
    )
    assert legacy["tdhf"]["multiplicity"] == "1"
    assert legacy["optimize"]["istate"] == "1"


@pytest.mark.parametrize(
    ("route", "short_driver", "explicit_driver", "section", "key"),
    [
        ("hf/6-31g", "grad", "grad(S0)", "properties", "grad"),
        ("dft/pbe0/6-31g", "grad", "grad(S0)", "properties", "grad"),
        ("hf/6-31g", "opt", "opt(S0)", "optimize", "istate"),
        ("dft/pbe0/6-31g", "opt", "opt(S0)", "optimize", "istate"),
    ],
)
def test_hf_dft_ground_driver_s0_is_the_omitted_default(
    route, short_driver, explicit_driver, section, key
):
    short = oqp_input.parse_canonical_oqp(
        '%s geom="h2o.xyz" %s' % (route, short_driver)
    )
    explicit = oqp_input.parse_canonical_oqp(
        '%s geom="h2o.xyz" %s' % (route, explicit_driver)
    )

    assert short.driver == explicit.driver
    assert oqp_input.render_canonical_oqp(short) == oqp_input.render_canonical_oqp(
        explicit
    )
    assert oqp_input.lower_to_legacy(short) == oqp_input.lower_to_legacy(explicit)
    assert oqp_input.lower_to_legacy(explicit)[section][key] == "0"


def test_concise_defaults_do_not_need_to_be_written():
    short = oqp_input.parse_canonical_oqp(
        'dft/pbe0/6-31g geom="h2o.xyz" opt'
    )
    explicit = oqp_input.parse_canonical_oqp(
        'dft/pbe0/6-31g geom="h2o.xyz" charge=0 mult=1 '
        'opt(S0,maxit=30)'
    )

    assert short.options == explicit.options
    assert short.driver == explicit.driver
    assert oqp_input.render_canonical_oqp(short) == oqp_input.render_canonical_oqp(
        explicit
    )
    assert "charge=" not in oqp_input.render_canonical_oqp(explicit)
    assert "mult=" not in oqp_input.render_canonical_oqp(explicit)
    assert "maxit" not in oqp_input.render_canonical_oqp(explicit)
    assert "maxit" not in oqp_input.lower_to_legacy(explicit)["optimize"]


def test_dftb_bundled_parameter_defaults_need_no_section_call():
    short = oqp_input.parse_canonical_oqp(
        'dftb geom="h2o.xyz" energy'
    )
    explicit_defaults = oqp_input.parse_canonical_oqp(
        'dftb geom="h2o.xyz" charge=0 energy '
        'dftb(backend=native,parameter_path="")'
    )

    assert oqp_input.render_canonical_oqp(short) == oqp_input.render_canonical_oqp(
        explicit_defaults
    )
    assert all(call.name != "dftb" for call in explicit_defaults.modifiers)
    assert "parameter_path" not in oqp_input.lower_to_legacy(short)["dftb"]


def test_single_line_and_line_oriented_inputs_are_identical():
    single = oqp_input.parse_canonical_oqp(
        'dft/pbe0/6-31g geom="h2o.xyz" charge=0 opt(S0) scf(conv=1e-8)'
    )
    multiline = oqp_input.parse_canonical_oqp(
        """\
dft/pbe0/6-31g
opt(S0)
scf(
  conv=1e-8
)
geom="h2o.xyz"
"""
    )

    assert multiline.model == single.model
    assert multiline.functional == single.functional
    assert multiline.basis == single.basis
    assert multiline.options == single.options
    assert multiline.driver == single.driver
    assert multiline.modifiers == single.modifiers
    assert oqp_input.lower_to_legacy(multiline) == oqp_input.lower_to_legacy(single)
    assert oqp_input.render_canonical_oqp(multiline) == (
        "dft/pbe0/6-31g\n"
        "opt\n"
        "scf(conv=1e-08)\n"
        'geom="h2o.xyz"\n'
    )


def test_inline_geometry_renders_one_atom_per_line_and_reparses():
    single_line = oqp_input.parse_canonical_oqp(
        'hf/sto-3g energy geom="O 0 0 0\\nH 0 0 1\\nH 0 1 0"'
    )
    multiline = oqp_input.parse_canonical_oqp(
        '''\
hf/sto-3g
energy
geom="""
O 0 0 0
H 0 0 1
H 0 1 0
"""
'''
    )

    assert multiline.options["geom"] == single_line.options["geom"]
    rendered = oqp_input.render_canonical_oqp(single_line)

    assert rendered == (
        "hf/sto-3g\n"
        'geom="""\n'
        "O 0 0 0\n"
        "H 0 0 1\n"
        "H 0 1 0\n"
        '"""\n'
    )
    reparsed = oqp_input.parse_canonical_oqp(rendered)
    assert reparsed.options["geom"] == single_line.options["geom"]
    assert oqp_input.lower_to_legacy(reparsed) == oqp_input.lower_to_legacy(
        single_line
    )


def test_harmless_whitespace_and_bare_calls_normalize_to_compact_input():
    spec, legacy = _parse(
        'mrsf(nstate = 3)/bhhlyp/6-31g* geom = "h2o.xyz" '
        'opt (S0, maxit = 100) scf (conv = 1e-8)'
    )
    assert legacy["optimize"]["istate"] == "1"
    assert legacy["optimize"]["maxit"] == "100"
    assert legacy["scf"]["conv"] == "1e-08"
    assert oqp_input.render_canonical_oqp(spec) == (
        "mrsf(nstate=3)/bhhlyp/6-31g*\n"
        "opt(S0,maxit=100)\n"
        "scf(conv=1e-08)\n"
        'geom="h2o.xyz"\n'
    )


def test_zero_argument_drivers_and_simple_modifiers_render_without_parentheses():
    spec, legacy = _parse(
        'dft/pbe0/def2-svp geom="h2o.xyz" energy() pcm nmr d4'
    )
    rendered = oqp_input.render_canonical_oqp(spec)
    assert rendered == (
        "dft/pbe0/def2-svp\n"
        "pcm\n"
        "nmr\n"
        "d4\n"
        'geom="h2o.xyz"\n'
    )
    assert oqp_input.render_canonical_oqp(
        oqp_input.parse_canonical_oqp(rendered)
    ) == rendered
    assert legacy["pcm"]["enabled"] == "True"
    assert legacy["properties"]["nmr_gauge"] == "giao"

    mrsf = oqp_input.parse_canonical_oqp(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" opt()'
    )
    assert "\nopt\n" in oqp_input.render_canonical_oqp(mrsf)


@pytest.mark.parametrize("geometry", ["h2o.xyz", '"water molecule.pdb"'])
def test_geometry_file_may_follow_the_route_positionally(geometry):
    spec, legacy = _parse(
        "mrsf(nstate=3)/bhhlyp/6-31g* %s opt" % geometry
    )
    expected = geometry.strip('"')
    assert spec.options["geom"] == expected
    assert legacy["input"]["system"] == expected
    assert oqp_input.render_canonical_oqp(spec) == (
        "mrsf(nstate=3)/bhhlyp/6-31g*\n"
        "opt\n"
        "geom=%s\n" % oqp_input._render_value(expected)
    )


def test_positional_and_named_geometry_cannot_be_combined():
    with pytest.raises(OQPInputError, match="Duplicate top-level option: geom"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp h2o.xyz geom="other.xyz" energy'
        )


def test_pcm_and_nmr_are_modifiers_not_primary_drivers():
    _, legacy = _parse(
        'dft/pbe0/def2-svp geom="h2o.xyz" grad(S0) '
        'pcm(water,model=ddpcm) nmr(gauge=giao)'
    )
    assert legacy["input"]["runtype"] == "grad"
    assert legacy["pcm"] == {
        "enabled": "True",
        "solvent": "water",
        "model": "ddpcm",
    }
    assert legacy["properties"]["scf_prop"] == "nmr"
    assert legacy["properties"]["nmr_gauge"] == "giao"
    _, default_nmr = _parse(
        'dft/pbe0/def2-svp geom="h2o.xyz" energy() nmr()'
    )
    assert default_nmr["properties"]["nmr_gauge"] == "giao"
    with pytest.raises(OQPInputError, match="PCM solvent is specified twice"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" energy() '
            'pcm(water,solvent=ethanol)'
        )


def test_conventional_tddft_s0_derivative_routes_to_ground_dft():
    _, legacy = _parse(
        'tddft(nstate=4)/pbe0/def2-svp geom="h2o.xyz" opt(S0)'
    )
    assert legacy["input"]["method"] == "hf"
    assert legacy["optimize"]["istate"] == "0"


def test_meci_and_mecp_use_physical_labels():
    _, meci = _parse(
        'mrsf(nstate=5)/bhhlyp/6-31g* geom="ethene.xyz" meci(S2,S1,maxit=4)'
    )
    assert meci["optimize"]["istate"] == "2"
    assert meci["optimize"]["jstate"] == "3"

    _, mecp = _parse(
        'mrsf(nstate=5)/bhhlyp/6-31g* geom="ethene.xyz" mecp(S1,T0)'
    )
    assert mecp["optimize"]["istate"] == "2"
    assert mecp["optimize"]["jstate"] == "1"
    assert mecp["optimize"]["imult"] == "1"
    assert mecp["optimize"]["jmult"] == "3"

    with pytest.raises(OQPInputError, match="distinct states"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=5)/bhhlyp/6-31g* geom="ethene.xyz" meci(S1,S1)'
        )


def test_sf_requires_an_explicit_root_not_a_pre_run_state_label():
    with pytest.raises(OQPInputError, match="use root=N"):
        oqp_input.parse_canonical_oqp(
            'sf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" grad(S0)'
        )
    _, legacy = _parse(
        'sf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" grad(root=1)'
    )
    assert legacy["properties"]["grad"] == "1"


def test_all_legacy_sections_can_be_keyword_calls_but_semantic_keys_are_owned():
    _, legacy = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" energy() '
        'input(perf=2) scf(maxit=80) tdhf(conv=1e-7) '
        'dftgrid(rad_npts=120) guess(type=huckel)'
    )
    assert legacy["input"]["perf"] == "2"
    assert legacy["scf"]["maxit"] == "80"
    assert legacy["tdhf"]["conv"] == "1e-07"
    assert legacy["dftgrid"]["rad_npts"] == "120"
    assert legacy["guess"]["type"] == "huckel"

    with pytest.raises(OQPInputError, match="route-owned"):
        oqp_input.parse_canonical_oqp(
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" energy() scf(multiplicity=1)'
        )
    with pytest.raises(OQPInputError, match="route-owned"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" mult=3 energy() scf(type=rohf)'
        )
    with pytest.raises(OQPInputError, match="duplicates the top-level perf"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" perf=1 energy() input(perf=2)'
        )


def test_lowered_values_are_accepted_by_configparser():
    _, legacy = _parse(
        'dft/pbe0/def2-svp geom="h2o.xyz" energy() '
        'input(perf=2) pcm(water) symmetry(enabled=true)'
    )
    assert all(isinstance(value, str) for section in legacy.values() for value in section.values())
    parser = configparser.ConfigParser()
    parser.read_dict(legacy)
    assert parser.getboolean("pcm", "enabled") is True


def test_string_enums_are_not_globally_coerced_to_bool_or_null():
    _, legacy = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" namd(S1,decoherence=off) '
        'scf(init_scf=no,init_basis=none)'
    )
    assert legacy["md"]["decoherence"] == "off"
    assert legacy["scf"]["init_scf"] == "no"
    assert legacy["scf"]["init_basis"] == "none"
    _, nac = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" nac(S0,S1,align=no)'
    )
    assert nac["nac"]["align"] == "no"

    spec = oqp_input.parse_canonical_oqp(
        'hf/sto-3g geom="h2.xyz" energy() scf(init_basis=null)'
    )
    rendered = oqp_input.render_canonical_oqp(spec)
    assert "init_basis=null" in rendered
    reparsed = oqp_input.parse_canonical_oqp(rendered)
    assert reparsed.modifiers[0].kwargs["init_basis"] is None


def test_paths_resolve_from_oqp_directory_not_process_cwd(tmp_path):
    _, legacy = _parse(
        'dft/pbe0/def2-svp geom="reactant.xyz" geom2="previous.xyz" '
        'neb(S0,product="product.xyz",images=7)',
        tmp_path,
    )
    assert legacy["input"]["system"] == str((tmp_path / "reactant.xyz").resolve())
    assert legacy["input"]["system2"] == str((tmp_path / "previous.xyz").resolve())
    assert legacy["neb"]["product"] == str((tmp_path / "product.xyz").resolve())
    assert legacy["neb"]["nimage"] == "7"


def test_compact_pdb_geometry_resolves_only_path_prefix(tmp_path):
    _, legacy = _parse(
        'dft/pbe0/def2-svp geom="protein model.pdb 9 10 17-19" energy',
        tmp_path,
    )
    assert legacy["input"]["system"] == (
        str((tmp_path / "protein model.pdb").resolve()) + " 9 10 17-19"
    )


def test_inline_geometry_is_not_mistaken_for_a_relative_filename(tmp_path):
    inline = "8 0.0 0.0 0.0 1 0.0 0.0 1.0"
    _, legacy = _parse(
        'hf/sto-3g geom="%s" energy()' % inline,
        tmp_path,
    )
    assert legacy["input"]["system"] == "8 0.0 0.0 0.0\n1 0.0 0.0 1.0"


def test_multiline_inline_geom2_keeps_its_first_atom():
    _, legacy = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* '
        'geom="8 0 0 0\\n1 0 0 1" geom2="8 0 0 0.1\\n1 0 0 1.1" '
        'nacme(S0,S1)'
    )
    assert legacy["input"]["system"].splitlines()[0].startswith("8 ")
    assert legacy["input"]["system2"].splitlines()[0].startswith("8 ")


def test_plain_natural_pcm_does_not_consume_the_next_word_as_solvent(tmp_path):
    result = oqp_input.resolve_oqp_text(
        "h2o.xyz에서 DFT/PBE0/def2-SVP PCM calculation을 수행한다",
        source_path=tmp_path / "pcm.oqp",
    )
    assert result.legacy_config["pcm"]["enabled"] == "True"
    assert "solvent" not in result.legacy_config["pcm"]


def test_natural_korean_request_writes_and_reparses_resolved_file(tmp_path):
    source = tmp_path / "water.oqp"
    source.write_text(
        "h2o.xyz를 사용한다. 전하 0. MRSF-TDDFT/BHHLYP/6-31G*로 "
        "5 singlet roots를 구하고 S1 구조 최적화를 최대 반복 100회 수행한다. "
        "SCF convergence=1e-8.",
        encoding="utf-8",
    )

    result = oqp_input.resolve_oqp_file(source)

    assert result.was_natural is True
    assert result.resolved_path == tmp_path / "water.resolved.oqp"
    assert result.resolved_path.read_text(encoding="utf-8") == result.canonical_text
    assert result.canonical_text == (
        "mrsf(nstate=5)/bhhlyp/6-31g*\n"
        "opt(S1,maxit=100)\n"
        "scf(conv=1e-08)\n"
        'geom="h2o.xyz"\n'
    )
    assert result.legacy_config["input"]["system"] == str((tmp_path / "h2o.xyz").resolve())
    assert result.legacy_config["optimize"]["istate"] == "2"


def test_natural_request_with_grad_and_opt_is_ambiguous():
    with pytest.raises(OQPInputError, match="multiple primary drivers"):
        oqp_input.compile_natural_request(
            "MRSF-TDDFT/BHHLYP/6-31G* h2o.xyz에서 S1 gradient와 구조 최적화"
        )


def test_natural_energy_correction_preserves_triplet_manifold():
    result = oqp_input.resolve_oqp_text(
        "MRSF-TDDFT/BHHLYP/6-31G* h2o.xyz에서 triplet energy calculation"
    )
    assert "energy(T0)" in result.canonical_text
    assert result.legacy_config["tdhf"]["multiplicity"] == "3"


def test_action_phrase_in_comment_does_not_turn_canonical_input_into_prose():
    result = oqp_input.resolve_oqp_text(
        'mrsf/bhhlyp/6-31g* geom="h2o.xyz" energy(T0) ! energy calculation'
    )
    assert result.was_natural is False
    assert "energy(T0)" in result.canonical_text


def test_canonical_looking_syntax_error_never_falls_back_to_natural_language():
    assert oqp_input.looks_canonical(
        'mrsf(nstate=5/bhhlyp/6-31g* geom="h2o.xyz" opt(S1)'
    )
    with pytest.raises(OQPInputError, match="Unclosed"):
        oqp_input.resolve_oqp_text(
            'mrsf(nstate=5/bhhlyp/6-31g* geom="h2o.xyz" opt(S1)'
        )


def test_explicit_nstate_must_cover_mrsf_internal_root():
    with pytest.raises(OQPInputError, match="response root 3, but nstate=2"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=2)/bhhlyp/6-31g* geom="h2o.xyz" grad(S2)'
        )


def test_branching_plane_uses_supported_nac_dispatch():
    _, legacy = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" bp(S0,S1)'
    )
    assert legacy["input"]["runtype"] == "nac"
    assert legacy["nac"]["bp"] == "True"


def test_nac_family_is_same_spin_and_nacme_requires_previous_geometry():
    with pytest.raises(OQPInputError, match="same spin manifold"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" nac(S0,T0)'
        )
    with pytest.raises(OQPInputError, match="nacme requires geom2"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" nacme(S0,S1)'
        )
    with pytest.raises(OQPInputError, match="not available for MRSF-TDDFTB"):
        oqp_input.parse_canonical_oqp(
            'mrsf-tddftb(nstate=3) geom="h2o.xyz" bp(S0,S1)'
        )
    with pytest.raises(OQPInputError, match="does not define option 'dt'"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" nac(S0,S1,dt=0.5)'
        )
    with pytest.raises(OQPInputError, match="does not define option 'nproc'"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" geom2="old.xyz" '
            'nacme(S0,S1,nproc=2)'
        )


def test_ir_and_raman_are_hessian_modifiers():
    spec, legacy = _parse(
        'dft/pbe0/def2-svp geom="h2o.xyz" hess(S0) ir() raman()'
    )
    assert [call.name for call in spec.modifiers] == ["ir", "raman"]
    assert legacy["input"]["runtype"] == "hess"
    with pytest.raises(OQPInputError, match="requires hess"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" energy() ir()'
        )


def test_qmmm_call_enables_qmmm_and_resolves_local_paths(tmp_path):
    (tmp_path / "custom.xml").write_text("<ForceField/>", encoding="utf-8")
    _, legacy = _parse(
        'dft/pbe0/def2-svp geom="qm.xyz" md(S0) '
        'qmmm(pdb_file="system.pdb",qm_atoms_xyz="qm.xyz",trajectory_file="traj.pdb",'
        'forcefield_files="amber14-all.xml,custom.xml")',
        tmp_path,
    )
    assert legacy["input"]["qmmm_flag"] == "True"
    assert legacy["qmmm"]["pdb_file"] == str((tmp_path / "system.pdb").resolve())
    assert legacy["qmmm"]["qm_atoms_xyz"] == str((tmp_path / "qm.xyz").resolve())
    assert legacy["qmmm"]["trajectory_file"] == str((tmp_path / "traj.pdb").resolve())
    assert legacy["qmmm"]["forcefield_files"] == (
        "amber14-all.xml," + str((tmp_path / "custom.xml").resolve())
    )
    assert legacy["properties"]["grad"] == "0"


def test_qmmm_pdb_file_is_inferred_from_the_last_geometry_line(tmp_path):
    text = """\
mrsf(nstate=3)/bhhlyp/6-31g*
namd(S1,soc=true,soc_basis=mch,nstep=200)
qmmm(forcefield_files="amber14-all.xml,amber14/tip3p.xml",qm_atoms="0-14")
geom="chromophore_water.pdb 0-14"
"""

    spec, legacy = _parse(text, tmp_path)

    pdb = str((tmp_path / "chromophore_water.pdb").resolve())
    assert legacy["input"]["system"] == pdb + " 0-14"
    assert legacy["qmmm"]["pdb_file"] == pdb
    assert "pdb_file" not in oqp_input.render_canonical_oqp(spec)

    redundant = oqp_input.parse_canonical_oqp(
        text.replace(
            'qmmm(forcefield_files=',
            'qmmm(pdb_file="chromophore_water.pdb",forcefield_files=',
        )
    )
    assert "pdb_file" not in oqp_input.render_canonical_oqp(redundant)

    distinct = oqp_input.parse_canonical_oqp(
        'dft/pbe0/def2-svp energy '
        'qmmm(pdb_file="environment.pdb") geom="qm.xyz"'
    )
    assert 'pdb_file="environment.pdb"' in oqp_input.render_canonical_oqp(distinct)


@pytest.mark.parametrize("spelling", ["nsteps", "n_steps"])
def test_qmmm_md_step_alias_lowers_to_active_config_key(spelling):
    _, legacy = _parse(
        'dft/pbe0/def2-svp geom="qm.xyz" md(S0) qmmm(%s=2)' % spelling
    )
    assert legacy["qmmm"]["n_steps"] == "2"
    assert "nsteps" not in legacy["qmmm"]


def test_qmmm_md_step_alias_collision_is_rejected():
    with pytest.raises(OQPInputError, match="same QM/MM MD step count"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="qm.xyz" md(S0) '
            'qmmm(nsteps=2,n_steps=3)'
        )


def test_md_requires_qmmm_and_physical_state_cannot_be_overridden():
    with pytest.raises(OQPInputError, match=r"requires qmmm\(\.\.\.\)"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" md(S0)'
        )
    with pytest.raises(OQPInputError, match="internal state selector"):
        oqp_input.parse_canonical_oqp(
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" namd(S1,active=1)'
        )
    with pytest.raises(OQPInputError, match="does not define option 'nstep'"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" md(S0,nstep=10) qmmm()'
        )
    _, energy = _parse(
        'dft/pbe0/def2-svp geom="ala.pdb 9 10" energy '
        'qmmm(pdb_file="system.pdb")'
    )
    assert energy["input"]["runtype"] == "energy"
    assert energy["input"]["qmmm_flag"] == "True"

    for unsupported in ("grad(S0)", "opt(S0)"):
        with pytest.raises(OQPInputError, match="not connected"):
            oqp_input.parse_canonical_oqp(
                'dft/pbe0/def2-svp geom="h2o.xyz" %s '
                'qmmm(pdb_file="system.pdb")' % unsupported
            )

    with pytest.raises(OQPInputError, match="not connected"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" grad(S0) '
            'input(qmmm_flag=true)'
        )


def test_soc_namd_numeric_active_surface_is_not_rewritten_as_mch_state():
    _, legacy = _parse(
        'mrsf(nstate=2)/bhhlyp/6-31g* geom="h2co.xyz" '
        'namd(active=5,soc=true,nstep=1)'
    )
    assert legacy["tdhf"]["nstate"] == "2"
    assert legacy["md"]["active"] == "5"
    assert legacy["md"]["soc"] == "True"
    assert "init_state" not in legacy["md"]


def test_d4_is_a_no_argument_modifier():
    _, legacy = _parse(
        'dft/pbe0/def2-svp geom="h2o.xyz" energy() d4()'
    )
    assert legacy["input"]["d4"] == "True"
    with pytest.raises(OQPInputError, match=r"d4\(\) does not accept"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" energy() d4(enabled=true)'
        )


def test_soc_uses_equal_nstate_or_explicit_singlet_triplet_counts():
    _, equal = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" soc'
    )
    assert equal["tdhf"]["nstate"] == "3"
    assert "nstate_s" not in equal["tdhf"]
    assert "nstate_t" not in equal["tdhf"]

    _, split = _parse(
        'mrsf/bhhlyp/6-31g* geom="h2o.xyz" soc(ns=3,nt=5)'
    )
    assert split["tdhf"]["nstate_s"] == "3"
    assert split["tdhf"]["nstate_t"] == "5"

    _, exact_split = _parse(
        'mrsf/bhhlyp/6-31g* geom="h2o.xyz" soc '
        'tdhf(nstate_s=3,nstate_t=5)'
    )
    assert exact_split["tdhf"]["nstate_s"] == "3"
    assert exact_split["tdhf"]["nstate_t"] == "5"

    with pytest.raises(OQPInputError, match="requires ns and nt together"):
        oqp_input.parse_canonical_oqp(
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" soc(ns=3)'
        )
    with pytest.raises(OQPInputError, match="Do not combine route nstate"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" soc(ns=3,nt=5)'
        )

    with pytest.raises(OQPInputError, match="specified twice"):
        oqp_input.parse_canonical_oqp(
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" soc(ns=3,nt=5) '
            'tdhf(nstate_s=3,nstate_t=5)'
        )
    with pytest.raises(OQPInputError, match="route nstate"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" soc '
            'tdhf(nstate_s=3,nstate_t=5)'
        )

    with pytest.raises(OQPInputError, match="requires nstate_s and nstate_t together"):
        oqp_input.parse_canonical_oqp(
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" soc tdhf(nstate_s=3)'
        )


@pytest.mark.parametrize("key", ["nstate_s", "nstate_t"])
@pytest.mark.parametrize("value", ["0", "-1", "2.5", "true", '"2"', "many"])
def test_exact_tdhf_soc_state_counts_require_positive_integers(key, value):
    counts = {"nstate_s": "3", "nstate_t": "5"}
    counts[key] = value
    with pytest.raises(OQPInputError, match=r"tdhf %s must be a positive integer" % key):
        oqp_input.parse_canonical_oqp(
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" soc '
            'tdhf(nstate_s=%s,nstate_t=%s)'
            % (counts["nstate_s"], counts["nstate_t"])
        )


def test_concise_soc_state_counts_reject_non_integer_numbers():
    with pytest.raises(OQPInputError, match="soc ns must be a positive integer"):
        oqp_input.parse_canonical_oqp(
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" soc(ns=2.5,nt=3)'
        )


def test_alias_duplicates_fail_instead_of_overwriting_by_order():
    duplicate_cases = [
        (
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" soc(soc_2e=0) '
            'input(soc_2e=1)',
            "two-electron treatment is specified twice",
        ),
        (
            'dft/pbe0/def2-svp geom="h2o.xyz" energy nmr '
            'properties(nmr_gauge=cgo)',
            "same NMR request",
        ),
        (
            'dft/pbe0/def2-svp geom="h2o.xyz" energy '
            'properties(scf_prop=mulliken) nmr',
            "same NMR request",
        ),
        (
            'dft/pbe0/def2-svp geom="h2o.xyz" data(scf_prop=mulliken) nmr',
            "same property",
        ),
        (
            'dft/pbe0/def2-svp geom="h2o.xyz" energy d4 input(d4=false)',
            "specified twice",
        ),
        (
            'dft/pbe0/def2-svp geom="h2o.xyz" energy qmmm=false '
            'qmmm(pdb_file="system.pdb")',
            "activation is specified more than once",
        ),
        (
            'dft/pbe0/def2-svp geom="h2o.xyz" energy '
            'input(qmmm_flag=true) qmmm(pdb_file="system.pdb")',
            "activation is specified more than once",
        ),
    ]
    for text, message in duplicate_cases:
        with pytest.raises(OQPInputError, match=message):
            oqp_input.parse_canonical_oqp(text)


def test_generic_section_typos_fail_early_with_suggestions():
    with pytest.raises(OQPInputError, match=r"scf\.maxdiiss.*Did you mean 'maxdiis'"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" energy scf(maxdiiss=9)'
        )
    with pytest.raises(OQPInputError, match=r"tests\.exeption.*Did you mean 'exception'"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" energy tests(exeption=true)'
        )


def test_qmmm_legacy_numeric_state_has_an_accurate_obsolete_diagnostic():
    with pytest.raises(OQPInputError, match="obsolete numeric selector"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" energy qmmm(istate=1)'
        )


@pytest.mark.parametrize("driver", ["grad", "opt", "hess", "data"])
def test_sf_state_aware_drivers_require_an_explicit_root(driver):
    with pytest.raises(OQPInputError, match="Specify an implementation root"):
        oqp_input.parse_canonical_oqp(
            'sf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" %s' % driver
        )

    _, legacy = _parse(
        'sf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" %s(root=1)' % driver
    )
    assert legacy["tdhf"]["type"] == "sf"


def test_irc_public_options_lower_to_owning_sections():
    _, legacy = _parse(
        'dft/pbe0/def2-svp geom="ts.xyz" '
        'irc(S0,direction=reverse,step=0.08,gtol=2e-5,maxit=20,hessian=analytical)'
    )
    assert legacy["input"]["runtype"] == "irc"
    assert legacy["optimize"]["istate"] == "0"
    assert legacy["optimize"]["maxit"] == "20"
    assert "lib" not in legacy["optimize"]
    assert legacy["oqp"]["irc_direction"] == "backward"
    assert legacy["oqp"]["irc_step"] == "0.08"
    assert legacy["oqp"]["path_gtol"] == "2e-05"
    assert legacy["hess"]["type"] == "analytical"

    with pytest.raises(OQPInputError, match="does not define option 'lib'"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="ts.xyz" irc(S0,lib=geometric,step=0.1)'
        )
    with pytest.raises(OQPInputError, match="must be numerical or analytical"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="ts.xyz" irc(S0,hessian=first)'
        )


def test_neb_workflow_options_lower_to_neb_and_oqp_sections(tmp_path):
    _, legacy = _parse(
        'dft/pbe0/def2-svp geom="reactant.xyz" '
        'neb(S0,product="product.xyz",images=7,spring=0.08,climb=true,'
        'fmax=0.002,frms=0.001,align=true,opt_ends=true,dt=0.4,'
        'output="final-path.xyz",maxit=12)',
        tmp_path,
    )
    assert legacy["neb"]["product"] == str((tmp_path / "product.xyz").resolve())
    assert legacy["neb"]["nimage"] == "7"
    assert legacy["oqp"]["spring"] == "0.08"
    assert legacy["oqp"]["climb"] == "True"
    assert legacy["oqp"]["fmax"] == "0.002"
    assert legacy["oqp"]["frms"] == "0.001"
    assert legacy["oqp"]["align"] == "True"
    assert legacy["oqp"]["opt_ends"] == "True"
    assert legacy["oqp"]["neb_dt"] == "0.4"
    assert legacy["oqp"]["neb_output"] == str((tmp_path / "final-path.xyz").resolve())
    assert legacy["optimize"]["maxit"] == "12"

    with pytest.raises(OQPInputError, match="legacy geomeTRIC NEB option"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="r.xyz" '
            'neb(S0,product="p.xyz",k=0.1)'
        )
    with pytest.raises(OQPInputError, match="climb must be true or false"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="r.xyz" '
            'neb(S0,product="p.xyz",climb=0.5)'
        )
    with pytest.raises(OQPInputError, match="climb_fmax.*greater than or equal"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="r.xyz" '
            'neb(S0,product="p.xyz",climb=true,fmax=0.01,climb_fmax=0.005)'
        )
    with pytest.raises(OQPInputError, match="climb_fmax.*greater than or equal"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="r.xyz" '
            'neb(S0,product="p.xyz",fmax=0.01) oqp(climb_fmax=0.005)'
        )


def test_geometry_drivers_are_native_and_mep_aliases_map_correctly():
    with pytest.raises(OQPInputError, match="legacy SciPy-backend option"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" meci(S0,S1,optimizer=bfgs)'
        )
    with pytest.raises(OQPInputError, match="hessian must be model, numerical, or analytical"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="ts.xyz" ts(S0,hessian=first)'
        )

    _, native = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" '
        'mep(S0,points=4,step=0.12,gtol=3e-5)'
    )
    assert native["optimize"]["maxit"] == "4"
    assert native["oqp"]["mep_step"] == "0.12"
    assert native["oqp"]["path_gtol"] == "3e-05"

    _, ts = _parse(
        'dft/pbe0/def2-svp geom="ts.xyz" '
        'ts(S0,hessian=numerical,coordsys=dlc,trust=0.1,trust_max=0.3,follow=1)'
    )
    assert ts["oqp"] == {
        "init_hessian": "numerical", "coordsys": "dlc", "trust": "0.1",
        "trust_max": "0.3", "follow": "1",
    }
    with pytest.raises(OQPInputError, match="does not define option 'lib'"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="ts.xyz" ts(S0,lib=geometric)'
        )


def test_readable_opt_supports_native_recovery_and_explicit_geometric_backend():
    native_text = """
mrsf-tddftb(nstate=3)
opt(S0,lib=oqp,maxit=100,auto_recovery=true,recovery_maxit=40,recovery_trust=0.02)
dftb(model=dtcam-tb)
geom="c60.xyz"
"""
    native_spec, native = _parse(native_text)
    assert native["optimize"]["lib"] == "oqp"
    assert native["oqp"]["auto_recovery"] == "True"
    assert native["oqp"]["recovery_maxit"] == "40"
    assert native["oqp"]["recovery_trust"] == "0.02"
    rendered_native = oqp_input.render_canonical_oqp(native_spec)
    assert "\nopt(S0,lib=oqp," in rendered_native
    assert rendered_native.rstrip().endswith('geom="c60.xyz"')

    geometric_text = """
mrsf-tddftb(nstate=3)
opt(S0,lib=geometric,maxit=200,coordsys=dlc,trust=0.02,tmax=0.05,
    convergence_set=GAU,hessian=never)
dftb(model=dtcam-tb)
geom="c60.xyz"
"""
    geometric_spec, geometric = _parse(geometric_text)
    assert geometric["optimize"]["lib"] == "geometric"
    assert geometric["geometric"] == {
        "coordsys": "dlc",
        "trust": "0.02",
        "tmax": "0.05",
        "convergence_set": "GAU",
        "hessian": "never",
    }
    reparsed_geometric = oqp_input.parse_canonical_oqp(
        oqp_input.render_canonical_oqp(geometric_spec)
    )
    assert reparsed_geometric.driver == geometric_spec.driver
    assert reparsed_geometric.modifiers == geometric_spec.modifiers
    assert reparsed_geometric.options == geometric_spec.options


def test_dftb_preset_allows_explicit_first_try_trah_mixer():
    """A DTCAM operator may select TRAH without disabling SCC or the preset."""
    spec, config = _parse("""
mrsf-tddftb(nstate=3)
dftb(model=dtcam-tb,scc_mixer=trah)
geom="alanine-dipeptide.xyz"
""")
    assert config["dftb"]["model"] == "dtcam-tb"
    assert config["dftb"]["scc_mixer"] == "trah"
    assert "scc_mixer=trah" in oqp_input.render_canonical_oqp(spec)


@pytest.mark.parametrize(
    "options,message",
    [
        ("lib=scipy", "must be oqp or geometric"),
        ("lib=geometric,trust=0.2,tmax=0.1", "0 < trust <= tmax"),
        ("lib=geometric,trust_max=0.2", "native lib=oqp option"),
        ("lib=oqp,tmax=0.1", "only with opt.*lib=geometric"),
        ("lib=oqp,auto_recovery=1", "must be true or false"),
        ("lib=oqp,recovery_maxit=0", "positive integer"),
    ],
)
def test_readable_opt_rejects_backend_option_mismatches(options, message):
    with pytest.raises(OQPInputError, match=message):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp opt(S0,%s) geom="h2o.xyz"' % options
        )


def test_baeka_is_a_variadic_meci_algorithm_with_safe_defaults():
    _, automatic_two = _parse(
        'mrsf(nstate=5)/bhhlyp/6-31g* geom="guess.xyz" meci(S0,S1)'
    )
    assert automatic_two["optimize"]["meci_search"] == "auto"

    _, automatic_with_penalty_controls = _parse(
        'mrsf-tddftb(nstate=2)\n'
        'meci(S0,S1,pen_incre=1.2,max_grad=0.01,'
        'rmsd_step=0.02,max_step=0.03)\n'
        'geom="guess.xyz"'
    )
    assert automatic_with_penalty_controls["optimize"]["meci_search"] == "auto"
    assert automatic_with_penalty_controls["optimize"]["pen_incre"] == "1.2"
    assert automatic_with_penalty_controls["optimize"]["max_grad"] == "0.01"
    assert automatic_two["optimize"]["states"] == "1,2"

    _, two = _parse(
        'mrsf(nstate=5)/bhhlyp/6-31g* geom="guess.xyz" '
        'meci(S0,S1,algorithm=baeka)'
    )
    assert two["input"]["runtype"] == "meci"
    assert two["optimize"]["meci_search"] == "baeka"
    assert two["optimize"]["states"] == "1,2"
    assert two["optimize"]["pen_sigma"] == "1.0"
    assert two["optimize"]["pen_alpha"] == "0.02"
    assert two["optimize"]["pen_delta"] == "0.025"
    assert two["optimize"]["pen_jump"] == "10,10,25,25,100,100,1000,1000,3000"
    assert two["optimize"]["energy_gap"] == "0.0001"

    spec, four = _parse(
        'mrsf(nstate=6)/bhhlyp/6-31g* geom="guess.xyz" '
        'meci(S3,S1,S0,S2,algorithm=baeka,sigma=2,alpha=0.03,'
        'delta_beta=0.1,beta_schedule="5,10",gap=2e-4)'
    )
    assert four["optimize"]["states"] == "1,2,3,4"
    assert four["optimize"]["pen_sigma"] == "2"
    assert four["optimize"]["pen_alpha"] == "0.03"
    assert four["optimize"]["pen_delta"] == "0.1"
    assert four["optimize"]["pen_jump"] == "5,10"
    assert four["optimize"]["energy_gap"] == "0.0002"
    assert "algorithm=baeka" in oqp_input.render_canonical_oqp(spec)

    _, inferred_three = _parse(
        'mrsf(nstate=5)/bhhlyp/6-31g* geom="guess.xyz" '
        'meci(S0,S1,S2)'
    )
    assert inferred_three["optimize"]["meci_search"] == "baeka"
    assert inferred_three["optimize"]["states"] == "1,2,3"


def test_baeka_rejects_ambiguous_or_nonconsecutive_requests():
    with pytest.raises(
        OQPInputError,
        match="more than two states require algorithm=auto or algorithm=baeka",
    ):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=5)/bhhlyp/6-31g* geom="g.xyz" '
            'meci(S0,S1,S2,algorithm=penalty)'
        )
    with pytest.raises(OQPInputError, match="consecutive response roots"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=5)/bhhlyp/6-31g* geom="g.xyz" '
            'meci(S0,S2,algorithm=baeka)'
        )
    with pytest.raises(OQPInputError, match="specified twice"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=5)/bhhlyp/6-31g* geom="g.xyz" '
            'meci(S0,S1,algorithm=baeka,meci_search=baeka)'
        )
    with pytest.raises(OQPInputError, match="alpha must be a positive"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=5)/bhhlyp/6-31g* geom="g.xyz" '
            'meci(S0,S1,algorithm=baeka,alpha=0)'
        )
    with pytest.raises(OQPInputError, match="gap_weight is fixed at 1.0"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=5)/bhhlyp/6-31g* geom="g.xyz" '
            'meci(S0,S1,algorithm=baeka,gap_weight=2)'
        )
    with pytest.raises(OQPInputError, match="additive delta_beta"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=5)/bhhlyp/6-31g* geom="g.xyz" '
            'meci(S0,S1,algorithm=baeka,pen_incre=1.2)'
        )
    with pytest.raises(OQPInputError, match="BaekA convergence does not use"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=5)/bhhlyp/6-31g* geom="g.xyz" '
            'meci(S0,S1,algorithm=baeka,max_step=0.01)'
        )


@pytest.mark.parametrize("value", ["2.5", "true"])
def test_count_and_multiplicity_fields_reject_lossy_integer_coercion(value):
    with pytest.raises(OQPInputError, match="nstate must be a positive integer"):
        oqp_input.parse_canonical_oqp(
            f'mrsf(nstate={value})/bhhlyp/6-31g* geom="g.xyz" energy'
        )
    with pytest.raises(OQPInputError, match="mult must be a positive integer"):
        oqp_input.parse_canonical_oqp(
            f'uhf/6-31g* geom="g.xyz" mult={value} energy'
        )
    with pytest.raises(OQPInputError, match="root must be a positive integer"):
        oqp_input.parse_canonical_oqp(
            f'sf/bhhlyp/6-31g* geom="g.xyz" grad(root={value})'
        )
    with pytest.raises(OQPInputError, match="maxit must be a positive integer"):
        oqp_input.parse_canonical_oqp(
            f'dft/pbe0/6-31g* geom="g.xyz" opt(S0,maxit={value})'
        )


def test_default_normalization_does_not_hide_lossy_numeric_types():
    with pytest.raises(OQPInputError, match="charge must be an integer"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/6-31g* geom="g.xyz" charge=0.0 energy'
        )
    with pytest.raises(OQPInputError, match="maxit must be a positive integer"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/6-31g* geom="g.xyz" opt(maxit=30.0)'
        )


def test_tci_retains_legacy_multiplicative_controls():
    _, legacy = _parse(
        'mrsf(nstate=5)/bhhlyp/6-31g* geom="guess.xyz" '
        'tci(S0,S1,S2,pen_sigma=2,pen_incre=1.2,energy_gap=0.002)'
    )
    optimize = legacy["optimize"]
    assert optimize["istate"] == "1"
    assert optimize["jstate"] == "2"
    assert optimize["kstate"] == "3"
    assert optimize["pen_sigma"] == "2"
    assert optimize["pen_incre"] == "1.2"
    assert optimize["energy_gap"] == "0.002"
    assert "states" not in optimize
    assert "meci_search" not in optimize
    assert "pen_delta" not in optimize
    assert "pen_jump" not in optimize

    for option in ("meci_search=baeka", "pen_delta=0.025", "pen_jump=10"):
        with pytest.raises(OQPInputError, match="does not define option"):
            oqp_input.parse_canonical_oqp(
                'mrsf(nstate=5)/bhhlyp/6-31g* geom="guess.xyz" '
                f'tci(S0,S1,S2,{option})'
            )


def test_native_minimum_accepts_frozen_distance_constraints():
    _, legacy = _parse(
        'dft/bhhlyp/3-21g geom="hcn.xyz" '
        'opt(S0,freeze="distance(1,2);r(1,3)",coordsys=dlc,trust=0.05)'
    )
    assert legacy["input"]["runtype"] == "optimize"
    assert legacy["oqp"]["freeze"] == "distance(1,2);r(1,3)"
    with pytest.raises(OQPInputError, match="distinct positive indices"):
        oqp_input.parse_canonical_oqp(
            'dft/bhhlyp/3-21g geom="hcn.xyz" opt(S0,freeze="distance(1,1)")'
        )
    with pytest.raises(OQPInputError, match="does not define option 'freeze'"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=3)/bhhlyp/3-21g geom="hcn.xyz" '
            'meci(S0,S1,freeze="distance(1,2)")'
        )


@pytest.mark.parametrize(
    "text, message",
    [
        ('dft/pbe0/def2-svp geom="h2o.xyz" energy oqp(init_hessian=analytical)',
         "not used by the native energy workflow"),
        ('dft/pbe0/def2-svp geom="ts.xyz" irc(S0) oqp(frms=0.1)',
         "not used by the native irc workflow"),
        ('dft/pbe0/def2-svp geom="h2o.xyz" opt(S0) oqp(trust=-1)',
         "0 < trust <= trust_max"),
        ('dft/pbe0/def2-svp geom="ts.xyz" ts(S0) oqp(follow=-1)',
         "non-negative mode index"),
        ('dft/pbe0/def2-svp geom="r.xyz" neb(S0,product="p.xyz") oqp(fmax=nan)',
         "positive number"),
        ('dft/pbe0/def2-svp geom="h2o.xyz" opt(S0,trust=nan)',
         "0 < trust <= trust_max"),
    ],
)
def test_native_exact_section_controls_cannot_be_ignored_or_invalid(text, message):
    with pytest.raises(OQPInputError, match=message):
        oqp_input.parse_canonical_oqp(text)


def test_ekt_parent_state_uses_physical_mrsf_label():
    _, legacy = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" ekt(S0,ip=true,ea=false)'
    )
    assert legacy["input"]["runtype"] == "ekt"
    assert legacy["tdhf"]["target"] == "1"
    assert legacy["ekt"] == {"ip": "True", "ea": "False"}


def test_mrsf_tddftb_target_spin_reaches_dftb_adapter():
    _, legacy = _parse(
        'mrsf-tddftb(nstate=4) geom="h2o.xyz" grad(T0) '
        'dftb(parameter_path="minimal.opdftb")'
    )
    assert legacy["tdhf"]["multiplicity"] == "3"
    assert legacy["dftb"]["target_multiplicity"] == "3"
    with pytest.raises(OQPInputError, match="not quintet"):
        oqp_input.parse_canonical_oqp(
            'mrsf-tddftb(nstate=4) geom="h2o.xyz" grad(Q1)'
        )


def test_dftb_route_types_and_open_shell_reference_reach_adapter():
    _, td = _parse(
        'tddftb(nstate=3) geom="h2o.xyz" energy(S0) '
        'dftb(parameter_path="minimal.opdftb")'
    )
    assert td["dftb"]["type"] == "tddftb"
    assert td["tdhf"]["type"] == "tda"

    _, ground = _parse(
        'dftb geom="radical.xyz" mult=2 energy() '
        'dftb(parameter_path="minimal.opdftb")'
    )
    assert ground["dftb"]["reference_multiplicity"] == "2"


@pytest.mark.parametrize(
    "route,expected_type",
    [
        ("dftb", "ground"),
        ("dftb0", "ground_noscc"),
        ("tddftb(nstate=2)", "tddftb"),
        ("tda-tddftb(nstate=2)", "tddftb"),
    ],
)
def test_default_dftb_reference_multiplicity_is_not_forwarded_to_probe(route, expected_type):
    _, legacy = _parse(
        '%s geom="h2o.xyz" energy() '
        'dftb(backend=probe,parameter_path="minimal.opdftb")' % route
    )

    assert legacy["dftb"]["type"] == expected_type
    assert "reference_multiplicity" not in legacy["dftb"]


def test_explicit_dftb_reference_multiplicity_is_preserved():
    _, legacy = _parse(
        'dftb geom="h2o.xyz" mult=1 energy() '
        'dftb(parameter_path="minimal.opdftb")'
    )
    assert legacy["dftb"]["reference_multiplicity"] == "1"


def test_dftb0_route_and_alias_lower_to_non_scc_ground_state():
    spec, legacy = _parse(
        'dftb-noscc geom="h2o.xyz" grad(S0) '
        'dftb(parameter_path="minimal.opdftb")'
    )

    assert spec.model == "dftb0"
    assert spec.physical_method == "DFTB0"
    assert oqp_input.render_canonical_oqp(spec).startswith("dftb0\n")
    assert legacy["input"]["method"] == "dftb"
    assert legacy["dftb"]["type"] == "ground_noscc"
    assert legacy["properties"]["grad"] == "0"


def test_tda_tddftb_route_is_canonical_tda_singlet_response():
    spec, legacy = _parse(
        'tddftb-tda(nstate=3) geom="h2o.xyz" grad(S1) '
        'dftb(parameter_path="minimal.opdftb")'
    )

    assert spec.model == "tda-dftb"
    assert spec.physical_method == "TD-DFTB (TDA)"
    assert oqp_input.render_canonical_oqp(spec).startswith("tda-tddftb(nstate=3)\n")
    assert legacy["tdhf"]["type"] == "tda"
    assert legacy["tdhf"]["multiplicity"] == "1"
    assert legacy["dftb"]["type"] == "tddftb"
    assert legacy["dftb"]["target_multiplicity"] == "1"
    assert legacy["properties"]["grad"] == "1"

    with pytest.raises(OQPInputError, match="supports singlet targets only"):
        oqp_input.parse_canonical_oqp(
            'tda-tddftb(nstate=3) geom="h2o.xyz" grad(T0)'
        )


def test_conventional_tddftb_triplet_is_rejected_with_spin_flip_guidance():
    with pytest.raises(OQPInputError, match="use sf-tddftb or mrsf-tddftb"):
        oqp_input.parse_canonical_oqp(
            'tddftb(nstate=3) geom="h2o.xyz" grad(T1)'
        )


@pytest.mark.parametrize(
    "filename,expected_model,expected_type",
    [
        ("H2_DFTB_ENERGY.oqp", "dftb", "ground"),
        ("H2O_TDDFTB_ENERGY.oqp", "tddftb", "tddftb"),
        ("CH2_SF-TDDFTB_ENERGY.oqp", "sf-dftb", "sf"),
        ("CH2_MRSF-TDDFTB_ENERGY.oqp", "mrsf-dftb", "mrsf"),
    ],
)
def test_minimal_dftb_family_examples_parse(filename, expected_model, expected_type):
    example_dir = ROOT / "examples" / "DFTB"
    text = (example_dir / filename).read_text(encoding="utf-8")

    spec = oqp_input.parse_canonical_oqp(text)
    legacy = oqp_input.lower_to_legacy(spec, source_dir=example_dir)

    assert spec.model == expected_model
    assert legacy["input"]["method"] == "dftb"
    assert legacy["dftb"]["type"] == expected_type


def test_sf_tddftb_route_uses_high_spin_reference_and_explicit_root():
    spec, legacy = _parse(
        'sftddftb(nstate=3) geom="h2o.xyz" grad(root=2) '
        'dftb(parameter_path="minimal.opdftb")'
    )

    assert spec.model == "sf-dftb"
    assert spec.physical_method == "SF-TDDFTB"
    assert spec.reference_method == "ROHF"
    assert spec.reference_multiplicity == 3
    assert oqp_input.render_canonical_oqp(spec).startswith("sf-tddftb(nstate=3)\n")
    assert legacy["scf"] == {"type": "rohf", "multiplicity": "3"}
    assert legacy["tdhf"]["type"] == "sf"
    assert legacy["tdhf"]["multiplicity"] == "1"
    assert legacy["dftb"]["type"] == "sf"
    assert legacy["dftb"]["target_multiplicity"] == "1"
    assert "reference_multiplicity" not in legacy["dftb"]
    assert legacy["properties"]["grad"] == "2"

    with pytest.raises(OQPInputError, match="use root=N"):
        oqp_input.parse_canonical_oqp(
            'sf-tddftb(nstate=3) geom="h2o.xyz" grad(S1)'
        )


def test_thermo_uses_supported_hessian_execution_path():
    _, legacy = _parse(
        'dft/pbe0/def2-svp geom="h2o.xyz" thermo(temperature=298.15)'
    )
    assert legacy["input"]["runtype"] == "hess"
    assert legacy["hess"]["state"] == "0"
    assert legacy["hess"]["temperature"] == "298.15"


def test_mrsf_prop_has_an_explicit_or_default_physical_state():
    _, explicit = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" prop(S1)'
    )
    _, default = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" prop'
    )
    assert explicit["properties"]["grad"] == "2"
    assert default["properties"]["grad"] == "1"


def test_primary_call_rejects_unknown_convenience_option_with_section_hint():
    with pytest.raises(OQPInputError, match="exact legacy section call"):
        oqp_input.parse_canonical_oqp(
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" opt(S0,rad_npts=100)'
        )


def test_method_driver_capability_errors_are_early_and_actionable():
    with pytest.raises(OQPInputError, match=r"MP2 currently supports energy\(\) only"):
        oqp_input.parse_canonical_oqp(
            'mp2/cc-pvdz geom="h2o.xyz" grad(S0)'
        )
    with pytest.raises(OQPInputError, match=r"UMRSF currently supports energy\(\) only"):
        oqp_input.parse_canonical_oqp(
            'umrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" opt(S0)'
        )
    with pytest.raises(OQPInputError, match="soc currently requires an MRSF route"):
        oqp_input.parse_canonical_oqp(
            'tddft(nstate=3)/pbe0/def2-svp geom="h2o.xyz" soc()'
        )


def test_unavailable_or_misspelled_derivative_types_fail_early():
    with pytest.raises(OQPInputError, match="analytical NAC is unavailable"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" '
            'nac(S0,S1,type=analytical)'
        )
    with pytest.raises(OQPInputError, match="type must be numerical or analytical"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" hess(S0,type=analyticla)'
        )


def test_common_input_mistakes_get_actionable_corrections():
    with pytest.raises(OQPInputError, match="leading '#' route marker"):
        oqp_input.parse_canonical_oqp(
            '# mrsf/bhhlyp/6-31g* geom="h2o.xyz" opt(S0)'
        )
    with pytest.raises(OQPInputError, match="leading '#' route marker"):
        oqp_input.resolve_oqp_text(
            '# mrsf/bhhlyp/6-31g* geom="h2o.xyz" opt(S0)'
        )
    with pytest.raises(OQPInputError, match="legacy .inp"):
        oqp_input.resolve_oqp_text('[input]\nmethod=hf\nsystem=h2o.xyz')
    with pytest.raises(OQPInputError, match="Did you mean 'mult'"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="h2o.xyz" multiplitcy=1 energy()'
        )
    with pytest.raises(OQPInputError, match="Raw istate"):
        oqp_input.parse_canonical_oqp(
            'mrsf/bhhlyp/6-31g* geom="h2o.xyz" istate=1 opt(S0)'
        )
