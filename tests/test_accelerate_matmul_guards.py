"""Finite matrix products must ignore stale Accelerate flags and fail closed."""

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def _load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def guarded_modules():
    saved_oqp_modules = {
        name: module
        for name, module in sys.modules.items()
        if name == "oqp" or name.startswith("oqp.")
    }
    metadata_tests = _load_module(
        "openqp_matmul_guard_metadata_tests", "tests/test_symmetry_metadata.py"
    )
    modules = (
        _load_module(
            "openqp_matmul_guard_symmetry", "pyoqp/oqp/library/symmetry.py"
        ),
        _load_module(
            "openqp_matmul_guard_transition_density",
            "pyoqp/oqp/analysis/transition_density.py",
        ),
        metadata_tests.load_molecule_module(),
    )
    yield modules
    for name in list(sys.modules):
        if name == "oqp" or name.startswith("oqp."):
            del sys.modules[name]
    sys.modules.update(saved_oqp_modules)


@pytest.mark.parametrize("warning_text", [
    "divide by zero encountered in matmul",
    "overflow encountered in matmul",
    "invalid value encountered in matmul",
])
def test_stale_accelerate_warning_is_suppressed(
    guarded_modules, monkeypatch, warning_text
):
    for module in guarded_modules:
        original_matmul = module.np.matmul

        def stale_fenv_matmul(left, right, *, _matmul=original_matmul):
            warnings.warn(warning_text, RuntimeWarning, stacklevel=2)
            return _matmul(left, right)

        with monkeypatch.context() as patch:
            patch.setattr(module.np, "matmul", stale_fenv_matmul)
            with warnings.catch_warnings():
                warnings.simplefilter("error", RuntimeWarning)
                result = module._finite_matmul(
                    np.eye(2), np.eye(2), "test matrix product"
                )
        np.testing.assert_array_equal(result, np.eye(2))


@pytest.mark.parametrize("bad_value", [np.nan, np.inf, -np.inf])
def test_nonfinite_matrix_product_fails_closed(guarded_modules, bad_value):
    left = np.eye(2)
    left[0, 0] = bad_value
    for module in guarded_modules:
        with pytest.raises(FloatingPointError, match="non-finite"):
            module._finite_matmul(left, np.eye(2), "test matrix product")


@pytest.mark.parametrize("bad_value", [np.nan, np.inf, -np.inf])
def test_nonfinite_backend_result_fails_closed(
    guarded_modules, monkeypatch, bad_value
):
    for module in guarded_modules:
        with monkeypatch.context() as patch:
            patch.setattr(
                module.np,
                "matmul",
                lambda _left, _right: np.full((2, 2), bad_value),
            )
            with pytest.raises(FloatingPointError, match="non-finite"):
                module._finite_matmul(
                    np.eye(2), np.eye(2), "test matrix product"
                )


def test_unrelated_runtime_warning_is_not_suppressed(guarded_modules, monkeypatch):
    for module in guarded_modules:
        original_matmul = module.np.matmul

        def unrelated_warning(left, right, *, _matmul=original_matmul):
            warnings.warn("unrelated numerical warning", RuntimeWarning, stacklevel=2)
            return _matmul(left, right)

        with monkeypatch.context() as patch:
            patch.setattr(module.np, "matmul", unrelated_warning)
            with pytest.raises(RuntimeWarning, match="unrelated numerical warning"):
                with warnings.catch_warnings():
                    warnings.simplefilter("error", RuntimeWarning)
                    module._finite_matmul(
                        np.eye(2), np.eye(2), "test matrix product"
                    )
