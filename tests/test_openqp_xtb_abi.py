"""C ABI v3 argument-order checks for the OpenQP-xTB adapter (no library).

Monkeypatches ctypes.CDLL with a recorder so OpenQPXTBAdapter._run_native
"calls" a fake openqp_xtb_state_gradient, then asserts the recorded argument
list against the pinned C ABI v3 of the openqp-xtb design spec (section 8a):

* the argument order is the inherited v2 openqp_dftb_state_gradient order;
* the int64 lc_gamma flag became lc_gamma_kind {yukawa:0, erf:1, ok:2} at the
  SAME position (right after lc_ground_state, right before zvector);
* exactly five new by-value scalars sit immediately AFTER the zvector flag
  and BEFORE scc_tolerance: int64 model (1 = gfn1), int64 dispersion, int64
  halogen_bond, int64 third_order, double spin_scale;
* the DFTB call (same recorder) is unchanged: the xtb argument list is
  exactly 5 entries longer.
"""

import ctypes
import importlib
from pathlib import Path
import sys
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
PYOQP = ROOT / "pyoqp"
if PYOQP.is_dir() and str(PYOQP) not in sys.path:
    sys.path.insert(0, str(PYOQP))


def _import_or_skip(module_name):
    try:
        return importlib.import_module(module_name)
    except RuntimeError as exc:
        text = str(exc)
        if (
            "Cannot locate OpenQP runtime files" in text
            or "OPENQP_ROOT does not contain" in text
            or "cannot load MPI library" in text
        ):
            for name in list(sys.modules):
                if name == "oqp" or name.startswith("oqp."):
                    sys.modules.pop(name, None)
            raise unittest.SkipTest(f"OpenQP native runtime is not available: {exc}") from exc
        raise


class _RecorderFunc:
    """Callable standing in for a shared-library symbol; records every call."""

    def __init__(self):
        self.calls = []
        self.restype = "unset"

    def __call__(self, *args):
        self.calls.append(args)


class _RecorderLib:
    def __init__(self, prefix):
        for suffix in ("state_gradient", "states_overlap", "soc_matrix"):
            setattr(self, f"{prefix}_{suffix}", _RecorderFunc())


class _FakeMol:
    """Water-sized molecule stand-in; enough surface for _run_native."""

    def __init__(self, method, section):
        self.config = {
            "input": {"method": method, "runtype": "energy", "charge": 0},
            "tdhf": {"type": "mrsf", "nstate": 2},
            "properties": {"grad": [0]},
            "scf": {"type": "rohf"},
            method: dict(section),
        }
        self.data = {"natom": 3}

    def get_atoms(self):
        return [8, 1, 1]

    def get_system(self):
        return [0.0, 0.0, 0.0, 0.0, 1.43, 1.11, 0.0, -1.43, 1.11]


# Distinctive values so positions are provable from the recorded call.
_SCC_TOLERANCE = 3.25e-9
_SPIN_SCALE = 1.25

_COMMON = {
    "backend": "native",
    "type": "mrsf",
    # existing file so the explicit library_path resolves; CDLL is mocked.
    "library_path": str(Path(__file__).resolve()),
    "parameter_path": "params_placeholder",
    "scc_tolerance": _SCC_TOLERANCE,
    "zvector": True,
}

_XTB_SECTION = dict(
    _COMMON,
    lc_gamma="ok",
    model="gfn1",
    dispersion=False,
    halogen_bond=True,
    third_order=False,
    spin_scale=_SPIN_SCALE,
)

_DFTB_SECTION = dict(_COMMON, lc_gamma="erf")

# v2 argument layout (openqp_dftb_state_gradient) up to the model splice:
#  0 natom  1 atoms*  2 coords*  3 param  4 len  5 method  6 len
#  7 state  8 nstate  9 need_grad  10 charge  11 scc_mixer  12 scc_history
# 13 max_scc_iterations  14 response_max_iterations  15 response_max_subspace
# 16 response_solver  17 reference_multiplicity  18 target_multiplicity
# 19 spin_complete  20 lc_ground_state  21 lc_gamma_kind  22 zvector
# [v3 model block spliced HERE for xtb]  then scc_tolerance, ...
_IDX_LC_GAMMA = 21
_IDX_ZVECTOR = 22
_IDX_MODEL = 23           # xtb only
_IDX_SCC_TOL_DFTB = 23
_IDX_SCC_TOL_XTB = 28


def _record_state_gradient(adapter_cls, mol):
    lib = _RecorderLib(adapter_cls.SYMBOL_PREFIX)
    adapter = adapter_cls(mol)
    with mock.patch.object(ctypes, "CDLL", lambda path: lib):
        adapter._run_native("mrsf", 1, need_grad=False)
    func = getattr(lib, f"{adapter_cls.SYMBOL_PREFIX}_state_gradient")
    # _native_library must have pinned restype on the looked-up symbol.
    assert func.restype is None
    (args,) = func.calls
    return args


class OpenQPXTBAbiV3Tests(unittest.TestCase):
    def setUp(self):
        self.openqp_dftb = _import_or_skip("oqp.library.openqp_dftb")
        self.openqp_xtb = _import_or_skip("oqp.library.openqp_xtb")

    def test_model_args_order_and_types(self):
        adapter = self.openqp_xtb.OpenQPXTBAdapter(_FakeMol("xtb", _XTB_SECTION))
        args = adapter._model_args()
        self.assertEqual(5, len(args))
        model, dispersion, halogen_bond, third_order, spin_scale = args
        for scalar in (model, dispersion, halogen_bond, third_order):
            self.assertIsInstance(scalar, ctypes.c_int64)
        self.assertIsInstance(spin_scale, ctypes.c_double)
        self.assertEqual(1, model.value)          # gfn1
        self.assertEqual(0, dispersion.value)     # False in _XTB_SECTION
        self.assertEqual(1, halogen_bond.value)   # True
        self.assertEqual(0, third_order.value)    # False
        self.assertAlmostEqual(_SPIN_SCALE, spin_scale.value)

        # dftb contributes nothing (its C ABI is unchanged).
        dftb_adapter = self.openqp_dftb.OpenQPDFTBAdapter(_FakeMol("dftb", _DFTB_SECTION))
        self.assertEqual([], dftb_adapter._model_args())

    def test_model_key_only_accepts_gfn1(self):
        adapter = self.openqp_xtb.OpenQPXTBAdapter(
            _FakeMol("xtb", dict(_XTB_SECTION, model="gfn2")))
        with self.assertRaisesRegex(ValueError, "model must be gfn1"):
            adapter._model_args()

    def test_state_gradient_argument_positions(self):
        xtb_args = _record_state_gradient(
            self.openqp_xtb.OpenQPXTBAdapter, _FakeMol("xtb", _XTB_SECTION))
        dftb_args = _record_state_gradient(
            self.openqp_dftb.OpenQPDFTBAdapter, _FakeMol("dftb", _DFTB_SECTION))

        # exactly five extra by-value scalars for xtb, nothing else moved.
        self.assertEqual(len(dftb_args) + 5, len(xtb_args))

        # shared prefix: same ctypes types up to and including zvector.
        for idx in range(_IDX_ZVECTOR + 1):
            self.assertIs(type(dftb_args[idx]), type(xtb_args[idx]), f"arg {idx}")

        # lc_gamma_kind at the v2 position: dftb erf -> 1, xtb ok -> 2.
        self.assertIsInstance(xtb_args[_IDX_LC_GAMMA], ctypes.c_int64)
        self.assertEqual(1, dftb_args[_IDX_LC_GAMMA].value)
        self.assertEqual(2, xtb_args[_IDX_LC_GAMMA].value)

        # zvector flag right before the splice point.
        self.assertIsInstance(xtb_args[_IDX_ZVECTOR], ctypes.c_int64)
        self.assertEqual(1, xtb_args[_IDX_ZVECTOR].value)

        # dftb: scc_tolerance directly follows zvector (v2 layout unchanged).
        self.assertIsInstance(dftb_args[_IDX_SCC_TOL_DFTB], ctypes.c_double)
        self.assertAlmostEqual(_SCC_TOLERANCE, dftb_args[_IDX_SCC_TOL_DFTB].value)

        # xtb: the five model scalars sit AFTER zvector and BEFORE scc_tolerance.
        model, dispersion, halogen_bond, third_order, spin_scale = \
            xtb_args[_IDX_MODEL:_IDX_MODEL + 5]
        for scalar in (model, dispersion, halogen_bond, third_order):
            self.assertIsInstance(scalar, ctypes.c_int64)
        self.assertIsInstance(spin_scale, ctypes.c_double)
        self.assertEqual(1, model.value)
        self.assertEqual(0, dispersion.value)
        self.assertEqual(1, halogen_bond.value)
        self.assertEqual(0, third_order.value)
        self.assertAlmostEqual(_SPIN_SCALE, spin_scale.value)

        self.assertIsInstance(xtb_args[_IDX_SCC_TOL_XTB], ctypes.c_double)
        self.assertAlmostEqual(_SCC_TOLERANCE, xtb_args[_IDX_SCC_TOL_XTB].value)

        # the tail (scc_tolerance onward) is identical between the backends.
        self.assertEqual(len(dftb_args) - _IDX_SCC_TOL_DFTB,
                         len(xtb_args) - _IDX_SCC_TOL_XTB)
        for off in range(len(dftb_args) - _IDX_SCC_TOL_DFTB):
            self.assertIs(type(dftb_args[_IDX_SCC_TOL_DFTB + off]),
                          type(xtb_args[_IDX_SCC_TOL_XTB + off]), f"tail {off}")

    def test_states_overlap_and_soc_keep_v2_symbols(self):
        # v2 signatures are kept: only the symbol prefix differs.
        lib = _RecorderLib("openqp_xtb")
        mol = _FakeMol("xtb", _XTB_SECTION)
        adapter = self.openqp_xtb.OpenQPXTBAdapter(mol)
        import numpy as np
        nbf, noca, nocb, nstate = 3, 2, 1, 1
        mo = np.eye(nbf).ravel()
        x = np.zeros(noca * (nbf - nocb) * nstate)
        xyz = np.asarray(mol.get_system(), dtype=float)
        with mock.patch.object(ctypes, "CDLL", lambda path: lib):
            adapter.states_overlap(xyz, xyz, mo, mo, x, x,
                                   noca=noca, nocb=nocb, multiplicity=1)
            adapter.soc_matrix(mo, x, x, noca=noca, nocb=nocb)
        self.assertEqual(1, len(lib.openqp_xtb_states_overlap.calls))
        self.assertEqual(1, len(lib.openqp_xtb_soc_matrix.calls))
        # v2 states_overlap signature: 21 arguments (see openqp_dftb adapter).
        self.assertEqual(21, len(lib.openqp_xtb_states_overlap.calls[0]))
        self.assertEqual(18, len(lib.openqp_xtb_soc_matrix.calls[0]))


if __name__ == "__main__":
    unittest.main()
