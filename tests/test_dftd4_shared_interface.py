"""Dependency-free checks for the replaceable, charge-aware DFT-D4 seam."""

from __future__ import annotations

import array
import ast
from pathlib import Path
from types import SimpleNamespace

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SINGLE_POINT = ROOT / "pyoqp" / "oqp" / "library" / "single_point.py"


class _FakeFFI:
    def new(self, cdecl, initializer=None):
        if cdecl == "int[]":
            return array.array("i", initializer)
        if cdecl == "double[]":
            values = [0.0] * initializer if isinstance(initializer, int) else initializer
            return array.array("d", values)
        if cdecl == "double*":
            return array.array("d", [0.0])
        if cdecl == "int*":
            return array.array("i", [0])
        raise AssertionError(f"unexpected cdecl: {cdecl}")

    @staticmethod
    def buffer(value, size):
        assert size <= value.buffer_info()[1] * value.itemsize
        return memoryview(value)


class _V2Library:
    def __init__(self):
        self.calls = []

    def oqp_dftd4_disp_v2(self, *args):
        self.calls.append(args)
        energy, gradient, status = args[-3:]
        energy[0] = -0.25
        for index in range(len(gradient)):
            gradient[index] = float(index + 1)
        status[0] = 0


class _LegacyLibrary:
    def __init__(self):
        self.calls = []

    def oqp_dftd4_disp(self, *args):
        self.calls.append(args)
        energy, gradient, status = args[-3:]
        energy[0] = -0.125
        status[0] = 0


def _load_dftd4_functions(library):
    """Execute only the two small D4 helpers, avoiding a native oqp import."""
    tree = ast.parse(SINGLE_POINT.read_text(encoding="utf-8"), filename=str(SINGLE_POINT))
    selected = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "_D4_DAMPING_KEYS"
            for target in node.targets
        ):
            selected.append(node)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in {
            "_dftd4_damping_values", "dftd4_native_disp"
        }:
            selected.append(node)
    namespace = {
        "np": np,
        "oqp": SimpleNamespace(ffi=_FakeFFI(), lib=library),
    }
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(SINGLE_POINT), "exec"), namespace)
    return namespace


def test_v2_receives_charge_and_explicit_damping():
    library = _V2Library()
    functions = _load_dftd4_functions(library)
    energy, gradient = functions["dftd4_native_disp"](
        [8, 1], np.zeros((2, 3)), "pbe", True,
        total_charge=-1.0,
        damping_params={
            "s6": 1.0, "s8": 2.0, "s9": 3.0,
            "a1": 4.0, "a2": 5.0, "alp": 16.0,
        },
    )

    assert energy == -0.25
    np.testing.assert_array_equal(gradient, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    call = library.calls[0]
    assert call[3] == -1.0
    assert call[6] == 1
    assert list(call[7]) == [1.0, 2.0, 3.0, 4.0, 5.0, 16.0]


def test_legacy_fallback_never_discards_charge_or_parameters():
    library = _LegacyLibrary()
    functions = _load_dftd4_functions(library)
    neutral_energy, _ = functions["dftd4_native_disp"](
        [1], np.zeros((1, 3)), "pbe", False
    )
    assert neutral_energy == -0.125
    assert len(library.calls) == 1

    try:
        functions["dftd4_native_disp"](
            [1], np.zeros((1, 3)), "pbe", False, total_charge=1.0
        )
    except RuntimeError as exc:
        assert "charge-aware DFT-D4" in str(exc)
    else:
        raise AssertionError("legacy DFT-D4 silently discarded molecular charge")


def test_shared_stack_packaging_contract_is_declared():
    external = (ROOT / "external" / "CMakeLists.txt").read_text(encoding="utf-8")
    source = (ROOT / "source" / "CMakeLists.txt").read_text(encoding="utf-8")
    header = (ROOT / "include" / "oqp.h").read_text(encoding="utf-8")
    manifest = (
        ROOT / "external" / "dftd4-corresponding-source" / "README.md"
    ).read_text(encoding="utf-8")

    assert '-DBUILD_SHARED_LIBS=ON' in external
    assert 'd4stack-${_OQP_DFTD4_LINKAGE}1' in external
    assert '-omp${ENABLE_OPENMP}' in external
    assert 'libdftd4.a' not in source
    assert 'libmulticharge.a' not in source
    assert 'libmctc-lib.a' not in source
    assert 'INSTALL_RPATH_USE_LINK_PATH FALSE' in source
    assert 'oqp_dftd4_disp_v2' in header
    assert 'mctc-lib-0.4.2/' in manifest
    assert 'multicharge-0.3.0/' in manifest
    assert 'dftd4-3.7.0/' in manifest


def test_high_level_callers_forward_actual_input_charge():
    source = SINGLE_POINT.read_text(encoding="utf-8")
    assert "total_charge = float(mol.config.get('input', {}).get('charge', 0))" in source
    assert "total_charge = float(self.mol.config.get('input', {}).get('charge', 0))" in source
    assert "total_charge=total_charge, damping_params=self.d4_param" in source


if __name__ == "__main__":
    test_v2_receives_charge_and_explicit_damping()
    test_legacy_fallback_never_discards_charge_or_parameters()
    test_shared_stack_packaging_contract_is_declared()
    test_high_level_callers_forward_actual_input_charge()
    print("DFT-D4 shared-interface unit tests passed")
