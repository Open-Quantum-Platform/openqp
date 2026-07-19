"""Regression tests for versioned openqp-dftb C ABI adapter layouts."""

from __future__ import annotations

import ctypes
import importlib.util
from pathlib import Path
import sys
import types
import unittest
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "pyoqp" / "oqp" / "library" / "openqp_dftb.py"


def _load_adapter_module():
    """Load the adapter with a dependency-light oqp package stub."""
    oqp = types.ModuleType("oqp")
    oqp.__path__ = []
    oqp.oqp_root = ""

    periodic_table = types.ModuleType("oqp.periodic_table")
    periodic_table.ELEMENTS_NAME = ["X", "H"]

    utils = types.ModuleType("oqp.utils")
    utils.__path__ = []
    constants = types.ModuleType("oqp.utils.constants")
    constants.ANGSTROM_TO_BOHR = 0.529177210903
    file_utils = types.ModuleType("oqp.utils.file_utils")
    file_utils.dump_log = lambda *args, **kwargs: None
    state_labels = types.ModuleType("oqp.utils.state_labels")
    state_labels.DFTB_CAP_STRUCTURED_TRACE = 1
    state_labels.DFTB_CAP_STATE_SPECTRUM = 2
    state_labels.DFTB_CAP_SCC_FINAL_TRUST = 4
    state_labels.canonical_dftb_type = lambda value: str(value).lower()
    state_labels.resolved_dftb_type = (
        lambda config: str(config.get("dftb", {}).get("type", "ground")).lower()
    )
    state_labels.dftb_method_name = lambda config: "DFTB"
    state_labels.public_state_label = (
        lambda config, root, **kwargs: "state %d" % int(root)
    )

    trace_spec = importlib.util.spec_from_file_location(
        "oqp.utils.dftb_trace",
        ROOT / "pyoqp" / "oqp" / "utils" / "dftb_trace.py",
    )
    assert trace_spec is not None and trace_spec.loader is not None
    dftb_trace = importlib.util.module_from_spec(trace_spec)
    trace_spec.loader.exec_module(dftb_trace)
    utils.dftb_trace = dftb_trace

    stubs = {
        "oqp": oqp,
        "oqp.periodic_table": periodic_table,
        "oqp.utils": utils,
        "oqp.utils.constants": constants,
        "oqp.utils.file_utils": file_utils,
        "oqp.utils.state_labels": state_labels,
        "oqp.utils.dftb_trace": dftb_trace,
    }
    saved = {name: sys.modules.get(name) for name in stubs}
    sys.modules.update(stubs)
    try:
        spec = importlib.util.spec_from_file_location(
            "_openqp_dftb_abi1_test_module", MODULE_PATH
        )
        assert spec is not None
        assert spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        return module
    finally:
        for name, previous in saved.items():
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous


ADAPTER_MODULE = _load_adapter_module()


class _FakeMolecule:
    def __init__(self, *, runtype="grad", dftb_updates=None):
        dftb = {
            "backend": "native",
            "type": "mrsf",
            "parameter_path": "/tmp/minimal.opdftb",
            "scc_mixer": "auto",
            "scc_history": 9,
            "max_scc_iterations": 321,
            "response_max_iterations": 41,
            "response_max_subspace": 73,
            "response_solver": "davidson",
            "reference_multiplicity": 3,
            "target_multiplicity": 1,
            "spin_complete": True,
            "lc_ground_state": False,
            "lc_gamma": "yukawa",
            "zvector": True,
            "scc_tolerance": 2.0e-9,
            "scc_mixing": 0.27,
            "scc_max_step": 0.44,
            "response_tolerance": 3.0e-7,
            "spc": 0.5,
            "spc_coco": 0.51,
            "spc_ovov": 0.52,
            "spc_coov": 0.53,
            "omega": 0.31,
            "cam_alpha": 0.02,
            "cam_beta": 0.98,
            "mrsf_shift_oo": 0.01,
            "mrsf_shift_co": 0.02,
            "mrsf_shift_ov": 0.03,
            "mrsf_shift_cv": 0.04,
        }
        if dftb_updates:
            dftb.update(dftb_updates)
        self.config = {
            "input": {
                "charge": -1,
                "runtype": runtype,
                "qmmm_flag": False,
            },
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "nstate": 2},
            "properties": {"grad": [2]},
            "optimize": {},
            "dftb": dftb,
        }
        self.data = {"natom": 2}

    @staticmethod
    def get_atoms():
        return np.array([1, 1], dtype=np.int64)

    @staticmethod
    def get_system():
        return np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.4], dtype=np.float64)


class _RecordingABI1Function:
    def __init__(self):
        self.calls = []
        self.restype = object()

    def __call__(self, *args):
        self.calls.append(args)
        if len(args) != 63:
            raise AssertionError(f"ABI v1 requires exactly 63 arguments, got {len(args)}")

        ctypes.cast(args[42], ctypes.POINTER(ctypes.c_double))[0] = -1.25
        ctypes.cast(args[43], ctypes.POINTER(ctypes.c_double))[0] = -1.05
        ctypes.cast(args[44], ctypes.POINTER(ctypes.c_double))[0] = 0.20
        ctypes.cast(args[45], ctypes.POINTER(ctypes.c_double))[0] = 0.02
        for index, value in enumerate((1.0, 2.0, 3.0, 4.0, 5.0, 6.0)):
            args[46][index] = value
        ctypes.cast(args[47], ctypes.POINTER(ctypes.c_int64))[0] = 2
        args[48][0], args[48][1] = -1.05, -0.95
        args[49][0], args[49][1] = 0.02, 0.03
        args[50][0], args[50][1] = 0.1, -0.1
        ctypes.cast(args[51], ctypes.POINTER(ctypes.c_int64))[0] = 2
        ctypes.cast(args[52], ctypes.POINTER(ctypes.c_int64))[0] = 1
        ctypes.cast(args[53], ctypes.POINTER(ctypes.c_int64))[0] = 1
        ctypes.cast(args[54], ctypes.POINTER(ctypes.c_int64))[0] = 1
        args[56][0], args[56][1] = -0.5, 0.2
        for index, value in enumerate((1.0, 0.0, 0.0, 1.0)):
            args[57][index] = value
        args[59][0], args[59][1] = 0.7, 0.8
        args[60].value = b""
        ctypes.cast(args[62], ctypes.POINTER(ctypes.c_int64))[0] = 0


class _FakeLibrary:
    def __init__(self):
        self.openqp_dftb_state_gradient = _RecordingABI1Function()
        self.openqp_dftb_states_overlap = object()
        self.openqp_dftb_soc_matrix = object()


class _IntegerProbe:
    def __init__(self, value):
        self.value = value
        self.restype = object()

    def __call__(self):
        return self.value


class _RecordingCurrentFunction:
    def __init__(self):
        self.calls = []
        self.restype = object()

    def __call__(self, *args):
        self.calls.append(args)
        if len(args) != 74:
            raise AssertionError(
                f"ABI v2/v3 requires exactly 74 arguments, got {len(args)}"
            )

        ctypes.cast(args[53], ctypes.POINTER(ctypes.c_double))[0] = -1.25
        ctypes.cast(args[54], ctypes.POINTER(ctypes.c_double))[0] = -1.05
        ctypes.cast(args[55], ctypes.POINTER(ctypes.c_double))[0] = 0.20
        ctypes.cast(args[56], ctypes.POINTER(ctypes.c_double))[0] = 0.02
        for index, value in enumerate((1.0, 2.0, 3.0, 4.0, 5.0, 6.0)):
            args[57][index] = value
        ctypes.cast(args[58], ctypes.POINTER(ctypes.c_int64))[0] = 2
        args[59][0], args[59][1] = -1.05, -0.95
        args[60][0], args[60][1] = 0.02, 0.03
        args[61][0], args[61][1] = 0.1, -0.1
        ctypes.cast(args[62], ctypes.POINTER(ctypes.c_int64))[0] = 2
        ctypes.cast(args[63], ctypes.POINTER(ctypes.c_int64))[0] = 1
        ctypes.cast(args[64], ctypes.POINTER(ctypes.c_int64))[0] = 1
        ctypes.cast(args[65], ctypes.POINTER(ctypes.c_int64))[0] = 1
        args[67][0], args[67][1] = -0.5, 0.2
        for index, value in enumerate((1.0, 0.0, 0.0, 1.0)):
            args[68][index] = value
        args[70][0], args[70][1] = 0.7, 0.8
        args[71].value = b""
        ctypes.cast(args[73], ctypes.POINTER(ctypes.c_int64))[0] = 0


class _FakeCurrentLibrary:
    def __init__(self, *, capabilities=7):
        self.openqp_dftb_state_gradient = _RecordingCurrentFunction()
        self.openqp_dftb_capi_abi_version = _IntegerProbe(3)
        if capabilities is not None:
            self.openqp_dftb_capi_capabilities = _IntegerProbe(capabilities)


class OpenQPDFTBABITests(unittest.TestCase):
    def _adapter_with_abi1(self, **mol_kwargs):
        mol = _FakeMolecule(**mol_kwargs)
        adapter = ADAPTER_MODULE.OpenQPDFTBAdapter(mol)
        library = _FakeLibrary()
        mol._openqp_dftb_cache["__native_library__"] = library
        mol._openqp_dftb_cache["__native_abi_version__"] = 1
        return adapter, library

    def _adapter_with_current(self, *, capabilities=7, **mol_kwargs):
        mol = _FakeMolecule(**mol_kwargs)
        adapter = ADAPTER_MODULE.OpenQPDFTBAdapter(mol)
        library = _FakeCurrentLibrary(capabilities=capabilities)
        mol._openqp_dftb_cache["__native_library__"] = library
        mol._openqp_dftb_cache["__native_abi_version__"] = 3
        mol._openqp_dftb_cache["__native_capabilities__"] = (
            0 if capabilities is None else capabilities
        )
        # These ABI layout tests load the adapter with a deliberately minimal
        # package stub, not the input-checker module that materializes the
        # production model default.
        mol._dftb_model_default_resolved = True
        return adapter, library

    def test_unversioned_library_is_classified_as_abi1(self):
        adapter, library = self._adapter_with_abi1()
        adapter.mol._openqp_dftb_cache.clear()
        adapter._native_library_path = lambda: Path("/tmp/libopenqp_dftb_c.fake")

        with mock.patch.object(ADAPTER_MODULE.ctypes, "CDLL", return_value=library):
            loaded = adapter._native_library()

        self.assertIs(loaded, library)
        self.assertEqual(adapter.mol._openqp_dftb_cache["__native_abi_version__"], 1)
        self.assertEqual(adapter.mol._openqp_dftb_cache["__native_capabilities__"], 0)
        self.assertIsNone(library.openqp_dftb_state_gradient.restype)
        self.assertEqual(len(library.openqp_dftb_state_gradient.argtypes), 63)

    def test_optional_capability_symbol_is_probed_and_cached(self):
        adapter, library = self._adapter_with_current()
        adapter.mol._openqp_dftb_cache.clear()
        adapter._native_library_path = lambda: Path("/tmp/libopenqp_dftb_c.fake")

        with mock.patch.object(ADAPTER_MODULE.ctypes, "CDLL", return_value=library):
            adapter._native_library()

        self.assertEqual(
            adapter.mol._openqp_dftb_cache["__native_capabilities__"], 7
        )
        self.assertIs(
            library.openqp_dftb_capi_capabilities.restype, ctypes.c_int64
        )
        self.assertEqual(len(library.openqp_dftb_state_gradient.argtypes), 74)

    def test_current_ctypes_signature_keeps_qmmm_fields_before_outputs(self):
        signature = ADAPTER_MODULE._state_gradient_argtypes(3)
        self.assertEqual(len(signature), 74)
        self.assertIs(signature[51], ctypes.c_int64)
        self.assertEqual(signature[52], ctypes.POINTER(ctypes.c_double))
        self.assertEqual(signature[53], ctypes.POINTER(ctypes.c_double))

    def test_abi4_signature_and_displaced_reference_anchor(self):
        signature = ADAPTER_MODULE._state_gradient_argtypes(4)
        self.assertEqual(len(signature), 78)
        # ABI-v4 inserts nbf/geometry/MOs/threshold immediately before the
        # historical external-potential pair and output tail.
        self.assertIs(signature[51], ctypes.c_int64)
        self.assertIs(signature[52], ctypes.POINTER(ctypes.c_double))
        self.assertIs(signature[53], ctypes.POINTER(ctypes.c_double))
        self.assertIs(signature[54], ctypes.c_double)
        self.assertIs(signature[55], ctypes.c_int64)

        adapter, _ = self._adapter_with_current()
        adapter.mol._openqp_dftb_cache["__native_abi_version__"] = 4
        anchor = {
            "signature": adapter._reference_anchor_signature("mrsf"),
            "nbf": 2,
            "coords": adapter.mol.get_system().copy(),
            "mos": np.eye(2),
        }
        adapter.mol._openqp_dftb_reference_anchor = anchor

        # Re-solving an unchanged geometry deliberately starts a strict SCC
        # solve; anchoring is only for the next displaced optimization point.
        nbf, coords, mos = adapter._native_reference_inputs(
            "mrsf", 4, need_grad=True
        )
        self.assertEqual(nbf, 0)
        self.assertEqual(coords.size, 1)
        self.assertEqual(mos.size, 1)

        adapter.mol.get_system = lambda: np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.5], dtype=np.float64
        )
        nbf, coords, mos = adapter._native_reference_inputs(
            "mrsf", 4, need_grad=True
        )
        self.assertEqual(nbf, 2)
        np.testing.assert_array_equal(coords, anchor["coords"])
        np.testing.assert_array_equal(mos, np.eye(2).ravel(order="F"))

    def test_missing_capability_symbol_means_zero_even_for_abi3(self):
        adapter, library = self._adapter_with_current(capabilities=None)
        adapter.mol._openqp_dftb_cache.clear()
        adapter._native_library_path = lambda: Path("/tmp/libopenqp_dftb_c.fake")

        with mock.patch.object(ADAPTER_MODULE.ctypes, "CDLL", return_value=library):
            adapter._native_library()

        self.assertEqual(
            adapter.mol._openqp_dftb_cache["__native_capabilities__"], 0
        )

    def test_cache_generation_preserves_native_capabilities(self):
        adapter, library = self._adapter_with_current()
        sentinel = object()
        adapter._cache_key = lambda *args, **kwargs: (
            ("H", "H"), b"new-geometry", (b"embedding",)
        )
        adapter._run_native = lambda *args, **kwargs: sentinel

        result = adapter._run_state("mrsf", 1, need_grad=False)

        self.assertIs(result, sentinel)
        self.assertIs(
            adapter.mol._openqp_dftb_cache["__native_library__"], library
        )
        self.assertEqual(
            adapter.mol._openqp_dftb_cache["__native_abi_version__"], 3
        )
        self.assertEqual(
            adapter.mol._openqp_dftb_cache["__native_capabilities__"], 7
        )

    def test_scc_failure_escalates_only_for_current_geometry(self):
        adapter, _ = self._adapter_with_current()
        result = object()
        calls = []

        def run_native(*args, **kwargs):
            mixer = adapter.dftb["scc_mixer"]
            calls.append(mixer)
            if mixer != "trust":
                raise RuntimeError(
                    "openqp-dftb restricted open-shell SCC cycle did not converge")
            return result

        adapter._run_native = run_native

        first = adapter._run_state("mrsf", 1, need_grad=False)
        self.assertIs(first, result)
        self.assertEqual(calls, ["auto", "broyden", "trust"])
        self.assertEqual(adapter.dftb["scc_mixer"], "auto")

        # Change geometry/cache generation.  Even though the previous geometry
        # needed TRAH, this one must start from the configured primary again.
        adapter.mol.get_system = lambda: np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.5], dtype=np.float64)
        calls.clear()
        second = adapter._run_state("mrsf", 1, need_grad=False)
        self.assertIs(second, result)
        self.assertEqual(calls, ["auto", "broyden", "trust"])
        self.assertEqual(adapter.dftb["scc_mixer"], "auto")

    def test_non_scc_native_failure_is_not_retried(self):
        adapter, _ = self._adapter_with_current()
        calls = []

        def run_native(*args, **kwargs):
            calls.append(adapter.dftb["scc_mixer"])
            raise RuntimeError("openqp-dftb parameter set is incomplete")

        adapter._run_native = run_native

        with self.assertRaisesRegex(RuntimeError, "parameter set"):
            adapter._run_state("mrsf", 1, need_grad=False)

        self.assertEqual(calls, ["auto"])
        self.assertEqual(adapter.dftb["scc_mixer"], "auto")

    def test_dtcam_auto_recovery_skips_duplicate_broyden_pass(self):
        adapter, _ = self._adapter_with_current(
            dftb_updates={"model": "dtcam-tb"}
        )
        self.assertEqual(
            adapter._scc_recovery_ladder("auto"), ("auto", "trust")
        )

    def test_explicit_trah_recovers_with_broyden_at_displaced_geometry(self):
        adapter, _ = self._adapter_with_current(
            dftb_updates={"scc_mixer": "trah"}
        )
        result = object()
        calls = []

        def run_native(*args, **kwargs):
            calls.append(adapter.dftb["scc_mixer"])
            if adapter.dftb["scc_mixer"] == "trah":
                raise RuntimeError(
                    "openqp-dftb restricted open-shell SCC cycle did not converge"
                )
            return result

        adapter._run_native = run_native
        self.assertIs(adapter._run_state("mrsf", 1, need_grad=False), result)
        self.assertEqual(calls, ["trah", "broyden"])
        self.assertEqual(adapter.dftb["scc_mixer"], "trah")

    def test_displaced_reference_starts_explicit_trah_on_broyden(self):
        adapter, _ = self._adapter_with_current(
            dftb_updates={"scc_mixer": "trah"}
        )
        result = object()
        calls = []
        adapter.mol._openqp_dftb_cache["__native_abi_version__"] = 4
        adapter._native_reference_inputs = lambda *args, **kwargs: (
            2, np.zeros(6), np.eye(2).ravel(order="F")
        )

        def run_native(*args, **kwargs):
            calls.append(adapter.dftb["scc_mixer"])
            return result

        adapter._run_native = run_native
        self.assertIs(adapter._run_state("mrsf", 1, need_grad=True), result)
        self.assertEqual(calls, ["broyden"])
        self.assertEqual(adapter.dftb["scc_mixer"], "trah")

    def test_large_explicit_trah_starts_on_charge_space_broyden(self):
        adapter, _ = self._adapter_with_current(
            dftb_updates={"scc_mixer": "trah"}
        )
        adapter.natom = 60
        result = object()
        calls = []

        def run_native(*args, **kwargs):
            calls.append(adapter.dftb["scc_mixer"])
            return result

        adapter._run_native = run_native
        self.assertIs(adapter._run_state("mrsf", 1, need_grad=False), result)
        self.assertEqual(calls, ["broyden"])
        self.assertEqual(adapter.dftb["scc_mixer"], "trah")

    def test_precontract_unversioned_layout_is_rejected_safely(self):
        adapter, library = self._adapter_with_abi1()
        adapter.mol._openqp_dftb_cache.clear()
        adapter._native_library_path = lambda: Path("/tmp/libopenqp_dftb_c.fake")
        del library.openqp_dftb_states_overlap
        del library.openqp_dftb_soc_matrix

        with mock.patch.object(ADAPTER_MODULE.ctypes, "CDLL", return_value=library):
            with self.assertRaisesRegex(
                RuntimeError, r"lacks the stable ABI-v1 marker.*cannot be identified safely"
            ):
                adapter._native_library()

    def test_abi1_call_uses_historical_63_argument_layout(self):
        adapter, library = self._adapter_with_abi1(
            dftb_updates={
                "c_mrsf": 0.71,
                "response_global_hybrid": True,
                "onsite_exchange_scale": 0.19,
            }
        )

        result = adapter._run_native("mrsf", 2, need_grad=True)

        self.assertEqual(len(library.openqp_dftb_state_gradient.calls), 1)
        args = library.openqp_dftb_state_gradient.calls[0]
        self.assertEqual(len(args), 63)

        # Stable ABI-v1 layout documented when ABI v2 was introduced in
        # openqp-dftb commit 26ddd5a (its exact predecessor is 4beaa864).
        self.assertIsInstance(args[0], ctypes.c_int64)
        self.assertEqual(args[0].value, 2)
        self.assertEqual([args[1][i] for i in range(2)], [1, 1])
        self.assertEqual([args[2][i] for i in range(6)], [0.0, 0.0, 0.0, 0.0, 0.0, 1.4])
        self.assertEqual(args[3], b"/tmp/minimal.opdftb")
        self.assertIsInstance(args[4], ctypes.c_int32)
        self.assertEqual(args[4].value, len(args[3]))
        self.assertEqual(args[5], b"mrsf")
        self.assertEqual(args[7].value, 2)    # root
        self.assertEqual(args[8].value, 2)    # requested roots
        self.assertEqual(args[9].value, 1)    # need gradient
        self.assertEqual(args[10].value, -1)  # molecular charge
        self.assertEqual(args[11].value, 3)   # auto SCC mixer
        self.assertEqual(args[12].value, 9)
        self.assertEqual(args[13].value, 321)
        self.assertEqual(args[14].value, 41)
        self.assertEqual(args[15].value, 73)
        self.assertEqual(args[16].value, 1)   # Davidson
        self.assertEqual(args[17].value, 3)   # high-spin reference
        self.assertEqual(args[18].value, 1)   # target multiplicity
        self.assertAlmostEqual(args[23].value, 2.0e-9)
        self.assertAlmostEqual(args[36].value, 0.04)
        self.assertAlmostEqual(args[37].value, 0.71)
        self.assertEqual(args[38].value, 1)
        self.assertAlmostEqual(args[39].value, 0.19)
        self.assertEqual(args[40].value, 0)
        self.assertIsInstance(args[61], ctypes.c_int32)
        self.assertEqual(args[61].value, 1024)

        self.assertEqual(result.reference_energy, -1.25)
        self.assertEqual(result.state_energy, -1.05)
        self.assertEqual(result.excitation_energy, 0.20)
        self.assertEqual(result.spin_square, 0.02)
        np.testing.assert_array_equal(
            result.gradient_bohr,
            np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        )
        np.testing.assert_array_equal(result.all_state_energies, [-1.05, -0.95])
        np.testing.assert_array_equal(result.all_spin_squares, [0.02, 0.03])
        np.testing.assert_array_equal(result.relaxed_charges, [0.1, -0.1])
        np.testing.assert_array_equal(result.mo_energies, [-0.5, 0.2])
        np.testing.assert_array_equal(result.mo_coefficients, np.eye(2))
        np.testing.assert_array_equal(result.response_vectors, [[0.7, 0.8]])

    def test_abi1_skips_the_production_model_default(self):
        # A C ABI v1 library cannot carry model presets, so the production
        # default must be skipped (NOT materialized then rejected): legacy
        # no-model inputs keep their explicit-keys meaning.
        adapter, library = self._adapter_with_abi1()
        adapter.mol.config["input"]["method"] = "dftb"

        adapter._resolve_model_default()

        self.assertEqual(adapter.dftb.get("model", ""), "")
        self.assertEqual(adapter._preset_bytes, b"")
        self.assertTrue(adapter.mol._dftb_model_default_resolved)
        # And the ABI-v1 call still proceeds without a model limitation.
        adapter._run_native("mrsf", 1, need_grad=False)
        self.assertEqual(len(library.openqp_dftb_state_gradient.calls), 1)

    def test_abi1_rejects_post_v1_operator_knobs_before_call(self):
        adapter, library = self._adapter_with_abi1(
            dftb_updates={"c_mrsf_oo": 0.7, "response_omega": 0.21}
        )

        with self.assertRaisesRegex(
            RuntimeError, r"C ABI v1.*c_mrsf_oo=0\.7.*response_omega=0\.21"
        ):
            adapter._run_native("mrsf", 1, need_grad=False)

        self.assertEqual(library.openqp_dftb_state_gradient.calls, [])

    def test_abi1_forwards_qmmm_embedding_and_relaxed_charges(self):
        adapter, library = self._adapter_with_abi1()
        adapter.mol.dftb_external_potential = np.array([0.1, -0.1])

        result = adapter._run_native("ground", 0, need_grad=True)

        args = library.openqp_dftb_state_gradient.calls[0]
        self.assertEqual(args[40].value, 2)
        self.assertEqual([args[41][i] for i in range(2)], [0.1, -0.1])
        np.testing.assert_array_equal(result.relaxed_charges, [0.1, -0.1])

    def test_current_call_uses_exact_74_argument_layout_and_slots(self):
        adapter, library = self._adapter_with_current(
            dftb_updates={
                "c_mrsf": 0.71,
                "response_global_hybrid": True,
                "onsite_exchange_scale": 0.19,
                "w_scale": 0.91,
                "response_w_scale": 0.81,
                "response_omega": 0.21,
                "response_cam_alpha": 0.11,
                "response_cam_beta": 0.61,
                "c_mrsf_oo": 0.41,
                "onsite_ss": 0.01,
                "onsite_sp": 0.02,
                "onsite_pp": 0.03,
                "model": "dtcam-tb",
            }
        )

        result = adapter._run_native("mrsf", 2, need_grad=True)

        calls = library.openqp_dftb_state_gradient.calls
        self.assertEqual(len(calls), 1)
        args = calls[0]
        self.assertEqual(len(args), 74)
        self.assertAlmostEqual(args[37].value, 0.71)
        self.assertEqual(args[38].value, 1)
        self.assertAlmostEqual(args[39].value, 0.19)
        self.assertAlmostEqual(args[40].value, 0.91)
        self.assertAlmostEqual(args[41].value, 0.81)
        self.assertAlmostEqual(args[42].value, 0.21)
        self.assertAlmostEqual(args[43].value, 0.11)
        self.assertAlmostEqual(args[44].value, 0.61)
        self.assertAlmostEqual(args[45].value, 0.41)
        self.assertAlmostEqual(args[46].value, 0.01)
        self.assertAlmostEqual(args[47].value, 0.02)
        self.assertAlmostEqual(args[48].value, 0.03)
        self.assertEqual(args[49], b"dtcam-tb")
        self.assertEqual(args[50].value, len(args[49]))
        self.assertEqual(args[51].value, 0)
        self.assertIsInstance(args[72], ctypes.c_int32)
        self.assertEqual(args[72].value, 1024)
        self.assertEqual(result.reference_energy, -1.25)
        np.testing.assert_array_equal(result.all_state_energies, [-1.05, -0.95])

    def test_current_layout_keeps_qmmm_potential_before_output_pointers(self):
        adapter, library = self._adapter_with_current()
        adapter.mol.dftb_external_potential = np.array([0.1, -0.1])

        adapter._run_native("ground", 0, need_grad=True)

        args = library.openqp_dftb_state_gradient.calls[0]
        # ABI-v2/v3 puts these between the preset fields (49--50) and the
        # reference-energy output pointer (53). Keeping the position explicit
        # prevents a native crash from a silently shifted ctypes call.
        self.assertEqual(args[51].value, 2)
        self.assertEqual([args[52][i] for i in range(2)], [0.1, -0.1])
        self.assertAlmostEqual(
            ctypes.cast(args[53], ctypes.POINTER(ctypes.c_double))[0], -1.25
        )


if __name__ == "__main__":
    unittest.main()
