"""Explicit CASPT2 direct energies pinned to the existing dense contraction."""
from types import SimpleNamespace

import numpy as np
import pytest

from test_caspt2_dyall import _run_caspt2


_H4 = "\nH 0 0 0\nH 0.03 0 0.74\nH 0 0.06 1.53\nH 0.04 0 2.31"
_CAS = {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1",
        "max_det": "5000", "max_memory": "512"}


@pytest.mark.parametrize("engine", ["direct", "python", "fortran"])
def test_options_and_preflight_accept_single_state_energy(engine):
    from oqp.library.caspt2_dyall import _caspt2_options
    from oqp.utils.input_checker import CheckReport, _check_pt2

    config = {"input": {"method": "caspt2", "runtype": "energy"},
              "pt2": {"engine": engine, "reference": "casci"}}
    options = _caspt2_options(config)
    assert options.engine == ("direct" if engine == "python" else engine)
    report = CheckReport()
    _check_pt2(config, report)
    assert not [item for item in report.diagnostics if item.severity == "ERROR"]


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
        from oqp.library.casci import _lib_backend
        lib, _ = _lib_backend()
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
        for threads in (1, 4):
            stream = direct._stream_fortran(h1e, eri, eps, norb, nc, na,
                                            ka, kb, c, options.max_terms, threads,
                                            max_memory=512)
            assert stream is not None, "missing native direct contraction"
            a, b, e0, v = stream
            positions = np.array([index[int(x) | (int(y) << norb)] for x, y in zip(a, b)])
            dense_v = (hfull @ vec)[positions]
            assert np.any(dense_v < -1e-8) and np.any(dense_v > 1e-8)
            np.testing.assert_allclose(v[0], dense_v, rtol=0, atol=1e-10)
            np.testing.assert_allclose(e0, h0.diagonal[positions], rtol=0, atol=1e-10)
        np.testing.assert_allclose(result["e2"][0], expected_e2, rtol=0, atol=1e-10)
        assert np.linalg.norm(expected_v) > 1e-6
        captured["e2"] = expected_e2
        return result

    monkeypatch.setattr(direct, "direct_qdpt2", pin)
    runner = _run_caspt2(tmp_path, "pin_" + engine, system=_H4, basis="sto-3g",
                         cas=_CAS, pt2={"reference": "casci", "engine": engine,
                                       "nproc": "1", **shift})
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


def test_direct_guards_are_explicit_errors(tmp_path):
    from oqp.library.qdpt2_direct import direct_qdpt2

    with pytest.raises(ValueError, match="63 orbitals"):
        direct_qdpt2(None, None, None, None, None, None, None, 0, 0,
                      (0, 0), 64, 0, [0], SimpleNamespace(ipea_shift=0), False)
    with pytest.raises(ValueError, match="max_terms"):
        _run_caspt2(tmp_path, "term_guard", system=_H4, basis="sto-3g",
                    cas=_CAS, pt2={"reference": "casci", "engine": "direct",
                                  "max_terms": "1", "nproc": "1"})
