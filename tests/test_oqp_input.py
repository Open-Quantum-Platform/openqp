"""Pure-Python tests for the markerless semantic .oqp input language."""

import configparser
import importlib.util
import sys
from pathlib import Path

import pytest


# Importing oqp normally loads liboqp.  The semantic parser intentionally has no
# native dependency, so exercise it directly to keep these tests fast and usable
# in source-only documentation/editor environments.
MODULE_PATH = Path(__file__).parents[1] / "pyoqp" / "oqp" / "utils" / "oqp_input.py"
SPEC = importlib.util.spec_from_file_location("_openqp_semantic_input", MODULE_PATH)
oqp_input = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = oqp_input
SPEC.loader.exec_module(oqp_input)

OQPInputError = oqp_input.OQPInputError


def _parse(text, source_dir=None):
    spec = oqp_input.parse_canonical_oqp(text)
    return spec, oqp_input.lower_to_legacy(spec, source_dir=source_dir)


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
    assert " mult=2 " in oqp_input.render_canonical_oqp(long)
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


def test_no_driver_defaults_to_explicitly_rendered_energy():
    spec, legacy = _parse('dft/pbe0/def2-svp geom="h2o.xyz" charge=0')
    assert spec.driver.name == "energy"
    assert spec.driver.explicit is False
    assert legacy["input"]["runtype"] == "energy"
    assert oqp_input.render_canonical_oqp(spec).endswith(" energy\n")


@pytest.mark.parametrize("driver", ["opt", "opt()"])
def test_mrsf_opt_without_state_defaults_to_s0(driver):
    _, legacy = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" %s' % driver
    )
    assert legacy["tdhf"]["multiplicity"] == "1"
    assert legacy["optimize"]["istate"] == "1"


def test_harmless_whitespace_and_bare_calls_normalize_to_compact_input():
    spec, legacy = _parse(
        'mrsf(nstate = 3)/bhhlyp/6-31g* geom = "h2o.xyz" '
        'opt (S0, maxit = 100) scf (conv = 1e-8)'
    )
    assert legacy["optimize"]["istate"] == "1"
    assert legacy["optimize"]["maxit"] == "100"
    assert legacy["scf"]["conv"] == "1e-08"
    assert oqp_input.render_canonical_oqp(spec) == (
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" '
        'opt(S0,maxit=100) scf(conv=1e-08)\n'
    )


def test_zero_argument_drivers_and_simple_modifiers_render_without_parentheses():
    spec, legacy = _parse(
        'dft/pbe0/def2-svp geom="h2o.xyz" energy() pcm nmr d4'
    )
    rendered = oqp_input.render_canonical_oqp(spec)
    assert rendered == (
        'dft/pbe0/def2-svp geom="h2o.xyz" energy pcm nmr d4\n'
    )
    assert oqp_input.render_canonical_oqp(
        oqp_input.parse_canonical_oqp(rendered)
    ) == rendered
    assert legacy["pcm"]["enabled"] == "True"
    assert legacy["properties"]["nmr_gauge"] == "giao"

    mrsf = oqp_input.parse_canonical_oqp(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" opt()'
    )
    assert oqp_input.render_canonical_oqp(mrsf).endswith(" opt\n")


@pytest.mark.parametrize("geometry", ["h2o.xyz", '"water molecule.pdb"'])
def test_geometry_file_may_follow_the_route_positionally(geometry):
    spec, legacy = _parse(
        "mrsf(nstate=3)/bhhlyp/6-31g* %s opt" % geometry
    )
    expected = geometry.strip('"')
    assert spec.options["geom"] == expected
    assert legacy["input"]["system"] == expected
    assert oqp_input.render_canonical_oqp(spec) == (
        'mrsf(nstate=3)/bhhlyp/6-31g* geom=%s opt\n'
        % oqp_input._render_value(expected)
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
        'mrsf(nstate=5)/bhhlyp/6-31g* geom="h2o.xyz" charge=0 '
        'opt(S1,maxit=100) scf(conv=1e-08)\n'
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
        'irc(S0,direction=reverse,step=0.08,maxit=20,hessian=analytical,lib=oqp)'
    )
    assert legacy["input"]["runtype"] == "irc"
    assert legacy["optimize"]["istate"] == "0"
    assert legacy["optimize"]["maxit"] == "20"
    assert legacy["optimize"]["lib"] == "oqp"
    assert legacy["oqp"]["irc_direction"] == "reverse"
    assert legacy["oqp"]["irc_step"] == "0.08"
    assert legacy["hess"]["type"] == "analytical"

    _, geometric = _parse(
        'dft/pbe0/def2-svp geom="ts.xyz" '
        'irc(S0,direction=forward,hessian=first,lib=geometric)'
    )
    assert geometric["geometric"]["hessian"] == "first"
    assert "hess" not in geometric

    with pytest.raises(OQPInputError, match="irc step=.*only with lib=oqp"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="ts.xyz" irc(S0,lib=geometric,step=0.1)'
        )
    with pytest.raises(OQPInputError, match="must be numerical or analytical"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="ts.xyz" irc(S0,lib=oqp,hessian=first)'
        )


def test_neb_workflow_options_lower_to_neb_and_oqp_sections(tmp_path):
    _, legacy = _parse(
        'dft/pbe0/def2-svp geom="reactant.xyz" '
        'neb(S0,product="product.xyz",images=7,spring=0.08,climb=true,maxit=12)',
        tmp_path,
    )
    assert legacy["neb"]["product"] == str((tmp_path / "product.xyz").resolve())
    assert legacy["neb"]["nimage"] == "7"
    assert legacy["oqp"]["spring"] == "0.08"
    assert legacy["oqp"]["climb"] == "True"
    assert legacy["optimize"]["maxit"] == "12"

    _, geometric = _parse(
        'dft/pbe0/def2-svp geom="reactant.xyz" '
        'neb(S0,product="product.xyz",images=7,lib=geometric,k=0.8,maxg=0.2,climb=0.5,align=true)',
        tmp_path,
    )
    assert geometric["neb"]["k"] == "0.8"
    assert geometric["neb"]["maxg"] == "0.2"
    assert geometric["neb"]["align"] == "True"
    assert geometric["neb"]["climb"] == "0.5"
    with pytest.raises(OQPInputError, match="does not belong to lib=geometric"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="r.xyz" '
            'neb(S0,product="p.xyz",lib=geometric,spring=0.1)'
        )


def test_optimizer_options_are_backend_specific_and_mep_aliases_map_correctly():
    with pytest.raises(OQPInputError, match="optimizer is used only with lib=scipy"):
        oqp_input.parse_canonical_oqp(
            'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" meci(S0,S1,optimizer=bfgs)'
        )
    with pytest.raises(OQPInputError, match="ts hessian=.*only with lib=geometric"):
        oqp_input.parse_canonical_oqp(
            'dft/pbe0/def2-svp geom="ts.xyz" ts(S0,lib=oqp,hessian=first)'
        )

    _, scipy = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" '
        'mep(S0,lib=scipy,points=4,step=0.12)'
    )
    assert scipy["optimize"]["mep_maxit"] == "4"
    assert scipy["optimize"]["step_size"] == "0.12"

    _, native = _parse(
        'mrsf(nstate=3)/bhhlyp/6-31g* geom="h2o.xyz" '
        'mep(S0,lib=oqp,points=4,step=0.12)'
    )
    assert native["optimize"]["maxit"] == "4"
    assert native["oqp"]["mep_step"] == "0.12"


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

    _, ground = _parse(
        'dftb geom="radical.xyz" mult=2 energy() '
        'dftb(parameter_path="minimal.opdftb")'
    )
    assert ground["dftb"]["reference_multiplicity"] == "2"


@pytest.mark.parametrize(
    "route,expected_type",
    [
        ("dftb", "auto"),
        ("dftb0", "ground_noscc"),
        ("tddftb(nstate=2)", "tddftb"),
        ("tda-tddftb(nstate=2)", "tda"),
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
    assert oqp_input.render_canonical_oqp(spec).startswith("dftb0 ")
    assert legacy["input"]["method"] == "dftb"
    assert legacy["dftb"]["type"] == "ground_noscc"
    assert legacy["properties"]["grad"] == "0"


def test_tda_tddftb_route_preserves_tda_and_physical_triplet_mapping():
    spec, legacy = _parse(
        'tddftb-tda(nstate=3) geom="h2o.xyz" grad(T0) '
        'dftb(parameter_path="minimal.opdftb")'
    )

    assert spec.model == "tda-dftb"
    assert spec.physical_method == "TDA-TDDFTB"
    assert oqp_input.render_canonical_oqp(spec).startswith("tda-tddftb(nstate=3) ")
    assert legacy["tdhf"]["type"] == "tda"
    assert legacy["tdhf"]["multiplicity"] == "3"
    assert legacy["dftb"]["type"] == "tda"
    assert legacy["dftb"]["target_multiplicity"] == "3"
    assert legacy["properties"]["grad"] == "1"

    with pytest.raises(OQPInputError, match="singlet or triplet"):
        oqp_input.parse_canonical_oqp(
            'tda-tddftb(nstate=3) geom="h2o.xyz" grad(Q0)'
        )


def test_sf_tddftb_route_uses_high_spin_reference_and_explicit_root():
    spec, legacy = _parse(
        'sftddftb(nstate=3) geom="h2o.xyz" grad(root=2) '
        'dftb(parameter_path="minimal.opdftb")'
    )

    assert spec.model == "sf-dftb"
    assert spec.physical_method == "SF-TDDFTB"
    assert spec.reference_method == "ROHF"
    assert spec.reference_multiplicity == 3
    assert oqp_input.render_canonical_oqp(spec).startswith("sf-tddftb(nstate=3) ")
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
