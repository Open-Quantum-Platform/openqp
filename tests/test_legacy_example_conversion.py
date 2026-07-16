"""Regression coverage for deterministic legacy-example migration."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
CONVERTER_PATH = ROOT / "tools" / "convert_legacy_examples.py"


def _load_converter():
    name = "_test_openqp_legacy_example_converter"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, CONVERTER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


converter = _load_converter()


def _run():
    return converter.convert_all(EXAMPLES)


def _result(run, relative: str):
    source = (EXAMPLES / relative).resolve()
    return next(result for result in run.results if result.source == source)


def test_every_current_legacy_example_is_convertible_and_committed():
    run = _run()
    legacy = sorted(EXAMPLES.rglob("*.inp"))

    assert not run.blockers, [
        "%s: %s" % (item.source.relative_to(ROOT), item.reason)
        for item in run.blockers
    ]
    assert len(run.results) == len(legacy)
    assert {result.source for result in run.results} == {
        path.resolve() for path in legacy
    }
    for result in run.results:
        assert result.target == result.source.with_suffix(".oqp")
        assert result.target.is_file(), result.target.relative_to(ROOT)
        assert result.target.read_text(encoding="utf-8") == result.text
        assert "\n" not in result.text.rstrip("\n")
        reparsed = converter.oqp_input.parse_canonical_oqp(result.text)
        assert converter.oqp_input.render_canonical_oqp(reparsed) == result.text
        # Compilation is also exercised: a syntactically valid line must still
        # lower into an executable legacy request.
        converter.oqp_input.lower_to_legacy(
            reparsed, source_dir=result.target.parent
        )

    assert not list(EXAMPLES.rglob("*.resolved.oqp"))


def test_shared_xyz_catalog_is_real_cartesian_geometry_only():
    run = _run()
    assets = run.catalog.ordered_assets()

    assert len(assets) == 27
    assert all(asset.name.endswith(".xyz") for asset in assets)
    assert all("pdb" not in asset.name.lower() for asset in assets)
    for asset in assets:
        path = EXAMPLES / "geometries" / asset.name
        assert path.read_text(encoding="utf-8") == asset.text
        lines = asset.text.splitlines()
        assert int(lines[0]) == len(lines[2:])
        assert converter._parse_cartesian_geometry("\n".join(lines[2:])) is not None


def test_response_roots_become_unambiguous_physical_states():
    run = _run()

    singlet = _result(
        run, "MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_GRADIENT.inp"
    ).text
    triplet = _result(
        run, "other/h2o_rohf_mrsf-t_6-31g_bhhlyp.inp"
    ).text
    quintet = _result(
        run, "other/h2o_rohf_mrsf-q_6-31g_bhhlyp.inp"
    ).text
    spin_flip = _result(
        run, "SF-TDDFT/H2O_BHHLYP-SFTDDFT_GRADIENT.inp"
    ).text

    assert "grad(S2)" in singlet
    assert "grad(T2)" in triplet
    assert "grad(Q2)" in quintet
    assert "grad(root=3)" in spin_flip


def test_native_workflows_use_public_concise_options():
    run = _run()
    constrained = _result(
        run, "OPT/HCN_RHF-DFT_CONSTRAINED_OQP.inp"
    ).text
    baeka = _result(
        run, "OPT/C2H4_BHHLYP-MRSFTDDFT_BAEKA_OQP.inp"
    ).text
    legacy_tci = _result(
        run, "OPT/C2H4_BHHLYP-MRSFTDDFT_TCI_OQP.inp"
    ).text
    mep = _result(
        run, "OPT/C2H4_BHHLYP-MRSFTDDFT_MEP.inp"
    ).text

    assert 'opt(S0,' in constrained
    assert 'freeze="distance(1,2)"' in constrained
    assert "geometric" not in constrained.lower()
    assert "meci(S0,S1,S2," in baeka
    assert "algorithm=baeka" in baeka
    assert "delta_beta=0.025" in baeka
    assert "beta_schedule=" in baeka
    assert "tci(S0,S1,S2," in legacy_tci
    assert "pen_incre=1.2" in legacy_tci
    assert "algorithm=baeka" not in legacy_tci
    assert "delta_beta=" not in legacy_tci
    assert "beta_schedule=" not in legacy_tci
    assert "mep(S1,maxit=2,step=0.1,gtol=0.003)" in mep
    assert "scipy" not in mep.lower()


def test_nmr_qmmm_and_nacme_edge_cases_preserve_intent():
    run = _run()
    nmr = _result(run, "NMR/H2O_RHF-NMR.inp")
    compact_qmmm = _result(run, "QMMM/ala.inp")
    qmmm_without_system = _result(run, "QMMM/run.inp")
    soc_namd = _result(
        run, "QMMM/H2CO-water_BHHLYP-SOC-NAMD-QMMM.inp"
    )
    nacme = _result(
        run, "other/h2o_nacme_rohf_mrsf-s_6-31g_bhhlyp.inp"
    ).text

    assert "properties(scf_prop=nmr)" in nmr.text
    assert " nmr" not in nmr.text
    lowered_nmr = converter.oqp_input.lower_to_legacy(
        converter.oqp_input.parse_canonical_oqp(nmr.text),
        source_dir=nmr.target.parent,
    )
    assert lowered_nmr["properties"]["scf_prop"] == "nmr"
    # Absence preserves the established legacy CGO default; bare concise nmr
    # would instead request the new GIAO default.
    assert "nmr_gauge" not in lowered_nmr["properties"]

    assert 'geom="ala.pdb 9 10 17 18 19"' in compact_qmmm.text
    lowered_compact = converter.oqp_input.lower_to_legacy(
        converter.oqp_input.parse_canonical_oqp(compact_qmmm.text),
        source_dir=compact_qmmm.target.parent,
    )
    assert lowered_compact["input"]["system"] == (
        str((compact_qmmm.target.parent / "ala.pdb").resolve())
        + " 9 10 17 18 19"
    )
    assert 'geom="water_dimer.pdb"' in qmmm_without_system.text
    assert "synthesized geom from qmmm.pdb_file" in qmmm_without_system.warnings
    lowered_md = converter.oqp_input.lower_to_legacy(
        converter.oqp_input.parse_canonical_oqp(qmmm_without_system.text),
        source_dir=qmmm_without_system.target.parent,
    )
    assert "md(S0)" in qmmm_without_system.text
    assert lowered_md["input"]["runtype"] == "md"
    assert "namd(active=5" in soc_namd.text
    lowered_soc_namd = converter.oqp_input.lower_to_legacy(
        converter.oqp_input.parse_canonical_oqp(soc_namd.text),
        source_dir=soc_namd.target.parent,
    )
    assert lowered_soc_namd["md"]["active"] == "5"
    assert "init_state" not in lowered_soc_namd["md"]
    assert "geom2=" in nacme
    assert "nacme(S0,S1)" in nacme


def test_ump2_reference_survives_example_conversion():
    run = _run()
    ump2 = _result(run, "MP2/h2o_ump2_6-31g.inp")

    assert ump2.text.startswith("mp2(reference=uhf)/6-31g ")
    lowered = converter.oqp_input.lower_to_legacy(
        converter.oqp_input.parse_canonical_oqp(ump2.text),
        source_dir=ump2.target.parent,
    )
    assert lowered["scf"]["type"] == "uhf"
    assert "reference" not in lowered.get("mp2", {})
