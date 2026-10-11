"""Explicit CASPT2 direct energies pinned to the existing dense contraction."""
from types import SimpleNamespace

import numpy as np
import pytest

from test_caspt2_dyall import _run_caspt2


_H4 = "\nH 0 0 0\nH 0.03 0 0.74\nH 0 0.06 1.53\nH 0.04 0 2.31"
_CAS = {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1",
        "max_det": "5000", "max_memory": "512"}


@pytest.mark.parametrize("engine", ["direct", "python", "fortran"])
@pytest.mark.parametrize("h0", [None, "fock", "caspt2"])
def test_options_and_preflight_accept_single_state_energy(engine, h0):
    from oqp.library.caspt2_dyall import _caspt2_options
    from oqp.utils.input_checker import CheckReport, _check_pt2

    config = {"input": {"method": "caspt2", "runtype": "energy"},
              "pt2": {"engine": engine, "reference": "casci"}}
    if h0 is not None:
        config["pt2"]["h0"] = h0
    options = _caspt2_options(config)
    assert options.engine == ("direct" if engine == "python" else engine)
    assert options.h0 == "fock"
    report = CheckReport()
    _check_pt2(config, report)
    assert not [item for item in report.diagnostics if item.severity == "ERROR"]


@pytest.mark.parametrize("engine", ["direct", "python", "fortran"])
@pytest.mark.parametrize("runtype", [None, "energy", "grad", "optimize"])
def test_schema_preserves_runtype_before_direct_validation(engine, runtype):
    from oqp.library.caspt2_dyall import _caspt2_options
    from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
    from oqp.utils.input_checker import CheckReport, _check_pt2
    from oqp.utils.input_parser import OQPConfigParser

    parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA)
    raw = {"input": {"method": "caspt2", "system": _H4},
           "pt2": {"reference": "casci", "engine": engine}}
    if runtype is not None:
        raw["input"]["runtype"] = runtype
    parser.load_dict(raw)
    config = parser.validate()
    assert config["input"]["runtype"] == (runtype or "energy")
    report = CheckReport()
    _check_pt2(config, report)
    engine_errors = [item for item in report.diagnostics
                     if item.severity == "ERROR" and item.path == "pt2.engine"]
    if runtype in {"grad", "optimize"}:
        assert engine_errors
        with pytest.raises(ValueError, match="single-state CASPT2 energy"):
            _caspt2_options(config)
    else:
        assert not engine_errors
        assert _caspt2_options(config).h0 == "fock"


@pytest.mark.parametrize("engine", ["direct", "python", "fortran", "dense"])
@pytest.mark.parametrize("targets", ["", "0", "1", "0,1"])
def test_single_state_target_count_matches_runtime(engine, targets):
    from oqp.library.caspt2_dyall import _caspt2_options, _reference_roots
    from oqp.utils.input_checker import CheckReport, _check_pt2

    config = {"input": {"method": "caspt2", "runtype": "energy"},
              "ci": {"nroot": 2},
              "pt2": {"engine": engine, "reference": "casci",
                      "target_roots": targets}}
    options = _caspt2_options(config)
    report = CheckReport()
    _check_pt2(config, report)
    errors = [item for item in report.diagnostics if item.severity == "ERROR"]
    if targets == "0,1":
        assert any(item.path == "pt2.target_roots" for item in errors)
        with pytest.raises(ValueError, match="single-state"):
            _reference_roots(options)
    else:
        assert not errors
        assert _reference_roots(options) == [int(targets or "0")]


@pytest.mark.parametrize("method,variant,h0,contraction", [
    ("mrmp2", "auto", "fock", "none"),
    ("mrmp2", "mrmp2", "fock", "none"),
    ("nevpt2", "auto", "dyall", "none"),
    ("sc-nevpt2", "auto", "dyall", "strong"),
    ("scnevpt2", "auto", "dyall", "strong"),
])
@pytest.mark.parametrize("targets", ["", "0", "1", "0,1"])
def test_other_single_state_pt2_targets_match_runtime(method, variant, h0,
                                                     contraction, targets):
    from oqp.library.caspt2_dyall import _caspt2_options, _reference_roots
    from oqp.utils.input_checker import CheckReport, _check_pt2

    config = {"input": {"method": method, "runtype": "energy"},
              "ci": {"nroot": 2},
              "pt2": {"engine": "auto", "reference": "casci", "variant": variant,
                      "h0": h0, "contraction": contraction,
                      "target_roots": targets}}
    options = _caspt2_options(config)
    assert options.variant == "caspt2", "these methods share the single-state runtime"
    report = CheckReport()
    _check_pt2(config, report)
    errors = [item for item in report.diagnostics if item.severity == "ERROR"]
    if targets == "0,1":
        assert any(item.path == "pt2.target_roots" for item in errors)
        with pytest.raises(ValueError, match="single-state"):
            _reference_roots(options)
    else:
        assert not errors
        assert _reference_roots(options) == [int(targets or "0")]


@pytest.mark.parametrize("method", ["casci", "casscf"])
@pytest.mark.parametrize("variant", ["caspt2", "mrmp2"])
@pytest.mark.parametrize("targets", ["1", "0,1"])
def test_single_state_reference_report_targets_match_runtime(method, variant, targets):
    from oqp.library.caspt2_dyall import _caspt2_options, _reference_roots
    from oqp.utils.input_checker import CheckReport, _check_pt2

    # The diagnostic caller selects its PT2 method through variant, rather
    # than input.method. It must obey the same single-state contract.
    config = {"input": {"method": method, "runtype": "energy"},
              "ci": {"nroot": 2},
              "pt2": {"reference_report": True, "reference": method,
                      "variant": variant,
                      "target_roots": targets}}
    options = _caspt2_options(config)
    report = CheckReport()
    _check_pt2(config, report)
    errors = [item for item in report.diagnostics if item.severity == "ERROR"]
    if targets == "0,1":
        assert any(item.path == "pt2.target_roots" for item in errors)
        with pytest.raises(ValueError, match="single-state"):
            _reference_roots(options)
    else:
        assert not errors
        assert _reference_roots(options) == [1]


@pytest.mark.parametrize("changes", [
    {"pt2": {"h0": "dyall"}},
    {"pt2": {"contraction": "strong"}},
    {"pt2": {"contraction": "ic"}},
    {"pt2": {"ipea_shift": "0.25"}},
    {"input": {"method": "ms-caspt2"}},
    {"input": {"method": "xms-caspt2"}},
    {"input": {"runtype": "grad"}},
])
def test_unsupported_requests_fail_before_calculation(changes):
    from oqp.library.caspt2_dyall import _caspt2_options
    from oqp.utils.input_checker import CheckReport, _check_pt2

    config = {"input": {"method": "caspt2", "runtype": "energy"},
              "pt2": {"engine": "direct", "reference": "casci"}}
    for section, values in changes.items():
        config[section].update(values)
    with pytest.raises(ValueError):
        _caspt2_options(config)
    report = CheckReport()
    _check_pt2(config, report)
    assert any(item.severity == "ERROR" and item.path == "pt2.engine"
               for item in report.diagnostics)


@pytest.mark.parametrize("engine", ["direct", "fortran"])
@pytest.mark.parametrize("shift", [{}, {"level_shift": "0.1"},
                                  {"imaginary_shift": "0.1"}, {"edshft": "0.02"}])
def test_signed_couplings_and_energy_match_dense(tmp_path, monkeypatch, engine, shift):
    import oqp
    from oqp.library import caspt2_dyall as dense
    from oqp.library import qdpt2_direct as direct

    assert oqp.BACKEND_AVAILABLE, "this regression requires the installed native backend"
    if engine == "fortran":
        from oqp import lib
        assert hasattr(lib, "qdpt2_stream_kernel"), "explicit Fortran engine must be exported"
    captured = {}
    original = direct.direct_qdpt2

    def pin(*args, **kwargs):
        h1e, eri, coeffs, energies, dets, eps, density, nc, na, ne, norb, enuc, roots, options = args
        result = original(*args, **kwargs)
        fdets, index, hfull, h0, external = dense._build_operators(
            h1e, eri, eps, nc, na, ne, norb, 5000, "fock", max_memory=512)
        vec = dense._embed_reference(coeffs[:, roots[0]], dets, nc, na, norb, index)
        expected_v, _, expected_e2, _ = dense._first_order(hfull, h0, vec, external, options)
        # Pin the native kernel on the production caller's own integrals and
        # signed CI coefficients. The bent, uneven geometry has no end-exchange
        # symmetry, so determinant/index permutations change these couplings.
        core = (1 << nc) - 1
        active = (1 << na) - 1
        ka = np.array([core | ((int(d) & active) << nc) for d in dets], dtype=np.uint64)
        kb = np.array([core | ((int(d) >> na) << nc) for d in dets], dtype=np.uint64)
        c = coeffs[:, roots].copy()
        # qdpt2_stream_kernel caps its team at nsup / 4. CAS(4,4)
        # supplies 36 determinants, so the four-thread call enters that path.
        assert c.shape == (36, 1)
        for threads in (1, 4):
            stream = direct._stream_fortran(h1e, eri, eps, norb, nc, na,
                                            ka, kb, c, options.max_terms, threads,
                                            max_memory=512)
            assert stream is not None, "missing native direct contraction"
            a, b, e0, v = stream
            positions = np.array([index[int(x) | (int(y) << norb)] for x, y in zip(a, b)])
            assert set(positions).issubset(set(external))
            coupled = set(external[np.abs((hfull @ vec)[external]) > 1e-12])
            assert coupled.issubset(set(positions))
            dense_v = (hfull @ vec)[positions]
            assert np.any(dense_v < -1e-8) and np.any(dense_v > 1e-8)
            np.testing.assert_allclose(v[0], dense_v, rtol=0, atol=1e-10)
            np.testing.assert_allclose(e0, h0.diagonal[positions], rtol=0, atol=1e-10)
        # Pin the NumPy enumerator independently, including completeness of
        # the coupled external space and the sign after reference contraction.
        a, b, val, e0, src = direct._stream_all(
            h1e, eri, eps, norb, nc, na, ka, kb, 1)
        merged = {}
        for x, y, coupling, energy, source in zip(a, b, val, e0, src):
            position = index[int(x) | (int(y) << norb)]
            assert position in external
            old_v, old_e = merged.get(position, (0.0, energy))
            assert abs(old_e - energy) < 1e-10
            merged[position] = (old_v + coupling * c[int(source), 0], energy)
        coupled = set(external[np.abs((hfull @ vec)[external]) > 1e-12])
        assert coupled.issubset(merged)
        positions = np.array(list(merged), dtype=int)
        np.testing.assert_allclose([merged[i][0] for i in positions],
                                   (hfull @ vec)[positions], rtol=0, atol=1e-10)
        np.testing.assert_allclose([merged[i][1] for i in positions],
                                   h0.diagonal[positions], rtol=0, atol=1e-10)
        np.testing.assert_allclose(result["e2"][0], expected_e2, rtol=0, atol=1e-10)
        assert np.linalg.norm(expected_v) > 1e-6
        captured["e2"] = expected_e2
        return result

    monkeypatch.setattr(direct, "direct_qdpt2", pin)
    cas = {**_CAS, "active_electrons": "4", "active_orbitals": "4",
           "frozen_core": "0"}
    runner = _run_caspt2(tmp_path, "pin_" + engine, system=_H4, basis="6-31g",
                         cas=cas, pt2={"reference": "casci", "engine": engine,
                                      "frozen": "0", "nproc": "1", **shift})
    np.testing.assert_allclose(runner.mol.data["OQP::CASPT2_STATE_SPECIFIC_CORRECTIONS"],
                               [captured["e2"]], rtol=0, atol=1e-10)
    log = (tmp_path / ("pin_" + engine + ".log")).read_text()
    assert "PyOQP PT2 engine:" in log and engine in log
    assert "generated external singles/doubles" in log
    assert "E(CASPT2" in log


def test_default_auto_keeps_dense_path(tmp_path, monkeypatch):
    from oqp.library import qdpt2_direct

    def unexpected(*args, **kwargs):
        pytest.fail("CASPT2 auto must retain the dense contraction")

    monkeypatch.setattr(qdpt2_direct, "direct_qdpt2", unexpected)
    _run_caspt2(tmp_path, "auto_dense", system=_H4, basis="sto-3g",
                cas=_CAS, pt2={"reference": "casci", "engine": "auto"})


@pytest.mark.parametrize("engine", ["direct", "fortran"])
def test_repeated_molecules_frozen_core_and_excited_root(tmp_path, engine):
    # Increase and decrease the orbital space in one process; use an excited
    # reference to check that target_roots is retained by the new dispatch.
    cases = [
        ("\nH 0 0 0\nH 0 0 0.74", {**_CAS, "frozen_core": "0"}, {}),
        (_H4, _CAS, {"target_roots": "1", "nroot": "2"}),
        ("\nLi 0 0 0\nH 0 0 1.6", _CAS, {"frozen": "1"}),
        ("\nH 0 0 0\nH 0 0 0.74", {**_CAS, "frozen_core": "0"}, {}),
    ]
    for i, (system, cas, extra) in enumerate(cases):
        results = []
        for route in ("dense", engine):
            runner = _run_caspt2(tmp_path, f"repeat_{i}_{route}", system=system,
                                 basis="sto-3g", cas=cas,
                                 ci={"nroot": "2"} if i == 1 else None,
                                 pt2={"reference": "casci", "engine": route,
                                      "nproc": "1", **extra})
            results.append(np.asarray(runner.mol.energies, dtype=float))
        np.testing.assert_allclose(results[1], results[0], rtol=0, atol=1e-9)


@pytest.mark.parametrize("engine", ["direct", "fortran"])
def test_direct_guards_are_explicit_errors(tmp_path, engine):
    from oqp.library.qdpt2_direct import direct_qdpt2

    with pytest.raises(ValueError, match="63 orbitals"):
        direct_qdpt2(None, None, None, None, None, None, None, 0, 0,
                      (0, 0), 64, 0, [0], SimpleNamespace(ipea_shift=0), False)
    with pytest.raises(ValueError, match="max_terms"):
        _run_caspt2(tmp_path, "term_guard", system=_H4, basis="sto-3g",
                    cas=_CAS, pt2={"reference": "casci", "engine": engine,
                                  "max_terms": "1", "nproc": "1"})


@pytest.mark.parametrize("engine", ["direct", "fortran"])
@pytest.mark.parametrize("reference", ["casci", "casscf"])
@pytest.mark.parametrize("has_reference_record", [False, True])
def test_term_guard_preserves_orbitals_for_same_molecule_retry(
        tmp_path, monkeypatch, engine, reference, has_reference_record):
    from oqp.library import caspt2_dyall

    runner = _run_caspt2(
        tmp_path, "before_failure", system=_H4, basis="sto-3g", cas=_CAS,
        pt2={"reference": reference, "engine": "dense", "nproc": "1"})
    mol = runner.mol
    key = "OQP::CASPT2_REFERENCE_MO"
    if has_reference_record:
        # A pre-existing gradient reference belongs to the completed operation.
        previous_reference = np.asarray(mol.data[key]).copy() + 0.125
        mol.data[key] = previous_reference
    else:
        del mol.data[key]
    mol.config["pt2"].update(engine=engine, max_terms=1)

    original_setup = caspt2_dyall._caspt2_setup
    original_semicanonicalize = caspt2_dyall._semicanonicalize
    snapshots = []
    semicanonical_snapshots = []

    def record_reference(*args, **kwargs):
        # Run the real CASCI/CASSCF setup and capture its orbitals before PT2's
        # semicanonicalization; no contraction or guard is mocked.
        setup = original_setup(*args, **kwargs)
        snapshots.append(np.asarray(mol.data["OQP::VEC_MO_A"]).copy())
        return setup

    def record_semicanonical(*args, **kwargs):
        result = original_semicanonicalize(*args, **kwargs)
        shape = np.asarray(mol.data["OQP::VEC_MO_A"]).shape
        semicanonical_snapshots.append(result[0].T.reshape(shape).copy())
        return result

    monkeypatch.setattr(caspt2_dyall, "_caspt2_setup", record_reference)
    monkeypatch.setattr(caspt2_dyall, "_semicanonicalize", record_semicanonical)
    with pytest.raises(ValueError, match="max_terms"):
        caspt2_dyall.native_caspt2_energy(mol)
    np.testing.assert_array_equal(mol.data["OQP::VEC_MO_A"], snapshots[-1])
    if has_reference_record:
        np.testing.assert_array_equal(mol.data[key], previous_reference)
    else:
        with pytest.raises(AttributeError, match="CASPT2_REFERENCE_MO"):
            mol.data[key]

    # The same molecule succeeds after raising the real term-count ceiling.
    mol.config["pt2"]["max_terms"] = 30_000_000
    caspt2_dyall.native_caspt2_energy(mol)
    retry_energy = np.asarray(mol.energies).copy()
    nbf = int(mol.data.get_basis()["nbf"])
    np.testing.assert_array_equal(
        mol.data[key], snapshots[-1].reshape((nbf, nbf)).T)
    np.testing.assert_array_equal(
        mol.data["OQP::VEC_MO_A"], semicanonical_snapshots[-1])
    monkeypatch.setattr(caspt2_dyall, "_caspt2_setup", original_setup)
    monkeypatch.setattr(caspt2_dyall, "_semicanonicalize", original_semicanonicalize)
    fresh = _run_caspt2(
        tmp_path, "fresh_after_failure", system=_H4, basis="sto-3g", cas=_CAS,
        pt2={"reference": reference, "engine": engine, "nproc": "1"})
    np.testing.assert_allclose(retry_energy, fresh.mol.energies, rtol=0, atol=1e-9)


def test_word_boundary_empty_external_and_memory_budget():
    from oqp.library.caspt2_dyall import CASPT2Options
    from oqp.library.qdpt2_direct import direct_qdpt2, _stream_fortran

    # Exact supported word boundary with dimensionally valid operands. An
    # empty-electron reference has no generated external configurations.
    n = 63
    h = np.zeros((n, n))
    g = np.zeros((n, n, n, n))
    eps = np.linspace(-1.0, 1.0, n)
    coeffs = np.ones((1, 1))
    options = CASPT2Options(engine="direct", nproc=1)
    result = direct_qdpt2(h, g, coeffs, np.zeros(1), [0], eps,
                         np.zeros((n, n)), 0, 0, (0, 0), n, 0,
                         [0], options, False, max_memory=512)
    assert result["n_external"] == 0 and result["e2"][0] == 0
    # Reject the working-memory requirement before allocating streamed terms.
    word = np.array([15], dtype=np.uint64)
    with pytest.raises(ValueError, match="max_memory"):
        _stream_fortran(h, g, eps, n, 0, 4, word, word,
                         coeffs, 30_000_000, 1, max_memory=1)


def test_populated_word_boundary_has_analytic_signed_couplings():
    from oqp.library.caspt2_dyall import CASPT2Options
    from oqp.library.qdpt2_direct import direct_qdpt2, _stream_all, _stream_fortran

    n, nc, na = 63, 61, 2
    eps = np.linspace(-3.0, 3.0, n)
    h = np.diag(eps)
    t0, t1 = 0.17, -0.23
    h[61, 0] = h[0, 61] = t0
    h[61, 1] = h[1, 61] = t1
    g = np.zeros((n, n, n, n))
    # Both spins occupy orbital 62 and every core orbital, leaving 61 empty.
    word = ((1 << nc) - 1) | (1 << 62)
    keys = np.array([word], dtype=np.uint64)
    coeffs = np.ones((1, 1))
    eref = 2 * (eps[:nc].sum() + eps[62])
    expected = {}
    for hole, signed in ((0, t0), (1, -t1)):
        excited = (word ^ (1 << hole)) | (1 << 61)
        for pair in ((excited, word), (word, excited)):
            expected[pair] = (signed, eref + eps[61] - eps[hole])

    a, b, v, e0, src = _stream_all(h, g, eps, n, nc, na, keys, keys, 1)
    assert np.all(src == 0)
    native = _stream_fortran(h, g, eps, n, nc, na, keys, keys,
                             coeffs, 30_000_000, 1, max_memory=512)
    assert native is not None, "the 63-orbital boundary requires the native kernel"
    an, bn, en, vn = native
    for aa, bb, vv, ee in ((a, b, v, e0), (an, bn, vn[0], en)):
        actual = {(int(x), int(y)): (coupling, energy)
                  for x, y, coupling, energy in zip(aa, bb, vv, ee)
                  if abs(coupling) > 1e-14}
        assert actual.keys() == expected.keys()
        assert all(x & (1 << 62) and y & (1 << 62) for x, y in actual)
        for pair, values in expected.items():
            np.testing.assert_allclose(actual[pair], values, rtol=0, atol=1e-10)

    e2 = -2 * (t0**2 / (eps[61] - eps[0]) + t1**2 / (eps[61] - eps[1]))
    for engine in ("direct", "fortran"):
        result = direct_qdpt2(h, g, coeffs, np.array([eref]), [0b1010], eps,
                             np.zeros((n, n)), nc, na, (1, 1), n, 0,
                             [0], CASPT2Options(engine=engine, nproc=1), False,
                             max_memory=512)
        np.testing.assert_allclose(result["e2"], [e2], rtol=0, atol=1e-12)


@pytest.mark.parametrize("engine", ["direct", "fortran"])
def test_casscf_reference_matches_dense(tmp_path, engine):
    results = []
    for route in ("dense", engine):
        runner = _run_caspt2(tmp_path, "casscf_" + route, system=_H4,
                             basis="sto-3g", cas=_CAS,
                             pt2={"reference": "casscf", "engine": route,
                                  "nproc": "1"})
        results.append(np.asarray(runner.mol.energies, dtype=float))
    np.testing.assert_allclose(results[1], results[0], rtol=0, atol=1e-9)


def test_explicit_fortran_refuses_missing_kernel(tmp_path, monkeypatch):
    from oqp.library import qdpt2_direct

    # The bridge returns None when the export is absent; the caller must fail
    # rather than silently substituting the NumPy engine.
    monkeypatch.setattr(qdpt2_direct, "_stream_fortran", lambda *a, **k: None)
    with pytest.raises(RuntimeError, match="lacks qdpt2_stream_kernel"):
        _run_caspt2(tmp_path, "missing_native", system=_H4, basis="sto-3g",
                    cas=_CAS, pt2={"reference": "casci", "engine": "fortran",
                                  "nproc": "1"})


def test_caspt2_numpy_process_pool_matches_dense(tmp_path, monkeypatch):
    from oqp.library import qdpt2_direct

    original = qdpt2_direct._stream_all
    selected = []

    def record(*args):
        # The production reference must reach the actual process-pool branch,
        # rather than the small-support serial fallback.
        selected.append((args[6].size, args[8]))
        return original(*args)

    monkeypatch.setattr(qdpt2_direct, "_stream_all", record)
    cas = {**_CAS, "active_electrons": "4", "active_orbitals": "4",
           "frozen_core": "0"}
    energies = []
    for route, workers in (("dense", "1"), ("direct", "1"), ("direct", "2")):
        runner = _run_caspt2(tmp_path, f"pool_{route}_{workers}", system=_H4,
                             basis="6-31g", cas=cas,
                             pt2={"reference": "casci", "engine": route,
                                  "nproc": workers, "frozen": "0"})
        energies.append(np.asarray(runner.mol.energies, dtype=float))
    assert any(support >= 8 and workers == 2 for support, workers in selected)
    np.testing.assert_allclose(energies[1], energies[0], rtol=0, atol=1e-9)
    np.testing.assert_allclose(energies[2], energies[0], rtol=0, atol=1e-9)
