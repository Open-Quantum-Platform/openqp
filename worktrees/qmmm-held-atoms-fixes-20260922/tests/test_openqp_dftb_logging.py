"""Pure-Python regression coverage for OpenQP-DFTB settings logging."""

import importlib.util
import os
import sys
import types
from pathlib import Path


ROOT = Path(__file__).parents[1]


def _load_adapter(monkeypatch, calls):
    oqp = types.ModuleType("oqp")
    oqp.__path__ = []
    oqp.oqp_root = ""
    monkeypatch.setitem(sys.modules, "oqp", oqp)

    utils = types.ModuleType("oqp.utils")
    utils.__path__ = []
    monkeypatch.setitem(sys.modules, "oqp.utils", utils)

    periodic = types.ModuleType("oqp.periodic_table")
    periodic.ELEMENTS_NAME = {1: "H"}
    monkeypatch.setitem(sys.modules, "oqp.periodic_table", periodic)

    constants = types.ModuleType("oqp.utils.constants")
    constants.ANGSTROM_TO_BOHR = 0.529177210903
    monkeypatch.setitem(sys.modules, "oqp.utils.constants", constants)

    file_utils = types.ModuleType("oqp.utils.file_utils")
    file_utils.dump_log = lambda *args, **kwargs: calls.append((args, kwargs))
    monkeypatch.setitem(sys.modules, "oqp.utils.file_utils", file_utils)

    labels_path = ROOT / "pyoqp" / "oqp" / "utils" / "state_labels.py"
    labels_spec = importlib.util.spec_from_file_location(
        "oqp.utils.state_labels", labels_path
    )
    labels = importlib.util.module_from_spec(labels_spec)
    monkeypatch.setitem(sys.modules, labels_spec.name, labels)
    labels_spec.loader.exec_module(labels)

    trace_path = ROOT / "pyoqp" / "oqp" / "utils" / "dftb_trace.py"
    trace_spec = importlib.util.spec_from_file_location(
        "oqp.utils.dftb_trace", trace_path
    )
    trace = importlib.util.module_from_spec(trace_spec)
    monkeypatch.setitem(sys.modules, trace_spec.name, trace)
    trace_spec.loader.exec_module(trace)
    utils.dftb_trace = trace

    adapter_path = ROOT / "pyoqp" / "oqp" / "library" / "openqp_dftb.py"
    adapter_spec = importlib.util.spec_from_file_location(
        "_openqp_dftb_logging_under_test", adapter_path
    )
    module = importlib.util.module_from_spec(adapter_spec)
    monkeypatch.setitem(sys.modules, adapter_spec.name, module)
    adapter_spec.loader.exec_module(module)
    return module


def _adapter_without_native_runtime(module):
    adapter = object.__new__(module.OpenQPDFTBAdapter)
    adapter.mol = types.SimpleNamespace()
    adapter.config = {
        "input": {"method": "dftb", "runtype": "energy"},
        "tdhf": {"type": "mrsf", "nstate": 3, "multiplicity": 1},
        "dftb": {"type": "mrsf", "backend": "native"},
    }
    adapter.dftb = adapter.config["dftb"]
    return adapter


def _load_input_checker(monkeypatch):
    oqp = types.ModuleType("oqp")
    oqp.__path__ = []
    utils = types.ModuleType("oqp.utils")
    utils.__path__ = []
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        size = 1
        use_mpi = False

    mpi_utils.MPIManager = MPIManager
    monkeypatch.setitem(sys.modules, "oqp", oqp)
    monkeypatch.setitem(sys.modules, "oqp.utils", utils)
    monkeypatch.setitem(sys.modules, "oqp.utils.mpi_utils", mpi_utils)

    path = ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py"
    spec = importlib.util.spec_from_file_location(
        "_openqp_dftb_input_checker_under_test", path
    )
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, spec.name, module)
    spec.loader.exec_module(module)
    return module


def test_adapter_logs_one_settings_block_for_all_states(monkeypatch):
    calls = []
    module = _load_adapter(monkeypatch, calls)
    adapter = _adapter_without_native_runtime(module)

    for _state in (1, 2, 3):
        adapter._log_settings_once(
            "mrsf",
            backend="native",
            parameter_path="/parameters/OB2W0PT3",
            library_path="/wheel/libopenqp_dftb_c.dylib",
            abi_version=3,
            capabilities=7,
        )

    assert len(calls) == 1
    assert calls[0][1]["section"] == "dftb"
    assert calls[0][1]["info"]["parameter_path"] == "/parameters/OB2W0PT3"
    assert calls[0][1]["info"]["abi_version"] == 3
    assert calls[0][1]["info"]["capabilities"] == 7


def test_excited_state_summary_uses_energy_category(monkeypatch):
    calls = []
    module = _load_adapter(monkeypatch, calls)
    adapter = _adapter_without_native_runtime(module)
    summary = {"states": "representative"}
    monkeypatch.setattr(adapter, "_excited_state_summary", lambda _base: summary)
    monkeypatch.setattr(adapter, "_resolved_method", lambda: "mrsf")
    monkeypatch.setattr(
        adapter, "_cache_key", lambda *_args, **_kwargs: ("summary",))
    monkeypatch.setattr(
        module.dftb_trace, "format_excited_state_summary",
        lambda value: "state summary" if value is summary else "unexpected",
    )
    monkeypatch.setattr(module, "dftb_method_name", lambda _config: "MRSF-TDDFTB")

    adapter._dump_excited_state_summary(types.SimpleNamespace(state=0))

    assert len(calls) == 1
    assert calls[0][1]["section"] == "dftb_state_summary"
    assert calls[0][1]["info"]["text"] == "state summary"


def test_native_diagnostics_capture_fd_and_restore_environment(monkeypatch):
    calls = []
    module = _load_adapter(monkeypatch, calls)
    monkeypatch.setenv("OPENQP_DFTB_SCC_TRACE", "previous")

    def native_writer():
        os.write(1, b"openqp_dftb_davidson_end iterations 4 elapsed_seconds 0.25\n")

    text = module._call_with_native_diagnostics(
        native_writer, print_level=2, state_spectrum=True
    )

    assert "davidson_end" in text
    assert os.environ["OPENQP_DFTB_SCC_TRACE"] == "previous"
    assert "OPENQP_DFTB_STATE_SPECTRUM" not in os.environ


def test_unadvertised_native_diagnostics_are_not_requested(monkeypatch):
    calls = []
    module = _load_adapter(monkeypatch, calls)
    seen = {}

    def inspect_environment():
        for name in module._NATIVE_DIAGNOSTIC_ENV:
            seen[name] = os.environ[name]

    module._call_with_native_diagnostics(
        inspect_environment,
        print_level=2,
        state_spectrum=False,
        structured_trace=False,
    )

    assert seen == {
        "OPENQP_DFTB_SCC_TRACE": "0",
        "OPENQP_DFTB_SOLVER_TRACE": "0",
        "OPENQP_ZVEC_TRACE": "0",
        "OPENQP_DFTB_STATE_SPECTRUM": "0",
    }


def test_bundled_parameter_resolution_is_silent(monkeypatch):
    calls = []
    module = _load_adapter(monkeypatch, calls)
    adapter = _adapter_without_native_runtime(module)
    monkeypatch.setattr(module, "_bundled_parameter_path", lambda: "/bundle/OB2W0PT3")

    assert adapter._parameter_path() == "/bundle/OB2W0PT3"
    assert adapter._parameter_path() == "/bundle/OB2W0PT3"
    assert calls == []


def test_legacy_conventional_tddftb_triplet_is_rejected(monkeypatch):
    input_checker = _load_input_checker(monkeypatch)
    for dftb_type, runtype, grad in (
        ("tddftb", "energy", []),
        ("auto", "grad", [1]),
    ):
        config = {
            "input": {"method": "dftb", "runtype": runtype},
            "tdhf": {"type": "tda", "nstate": 3, "multiplicity": 3},
            "dftb": {
                "type": dftb_type,
                "backend": "native",
                "parameter_path": "/tmp/minimal.opdftb",
                "target_multiplicity": 3,
            },
            "properties": {"grad": grad, "scf_prop": [], "nac": []},
            "optimize": {"istate": 0, "jstate": 0},
            "pcm": {"enabled": False},
        }

        report = input_checker.CheckReport()
        input_checker._check_dftb(config, report)

        assert not report.ok
        assert "singlet TDA response only" in report.to_text()
        assert "SF-TDDFTB or MRSF-TDDFTB" in report.to_text()


def test_conventional_tddftb_open_shell_reference_is_rejected(monkeypatch):
    input_checker = _load_input_checker(monkeypatch)
    config = {
        "input": {"method": "dftb", "runtype": "energy"},
        "tdhf": {"type": "tda", "nstate": 1, "multiplicity": 1},
        "dftb": {
            "type": "tddftb",
            "backend": "native",
            "parameter_path": "/tmp/minimal.opdftb",
            "reference_multiplicity": 3,
            "target_multiplicity": 1,
        },
        "properties": {"grad": [], "scf_prop": [], "nac": []},
        "optimize": {"istate": 0, "jstate": 0},
        "pcm": {"enabled": False},
    }

    report = input_checker.CheckReport()
    input_checker._check_dftb(config, report)

    assert not report.ok
    assert "closed-shell singlet reference" in report.to_text()


def test_explicit_dftb_response_type_must_match_tdhf_type(monkeypatch):
    input_checker = _load_input_checker(monkeypatch)
    for dftb_type, td_type, expected in (
        ("tddftb", "mrsf", "rpa/tda"),
        ("sf", "tda", "sf"),
        ("mrsf", "tda", "mrsf"),
    ):
        config = {
            "input": {"method": "dftb", "runtype": "energy"},
            "tdhf": {"type": td_type, "nstate": 1, "multiplicity": 1},
            "dftb": {
                "type": dftb_type,
                "backend": "native",
                "parameter_path": "/tmp/minimal.opdftb",
                "target_multiplicity": 1,
            },
            "properties": {"grad": [], "scf_prop": [], "nac": []},
            "optimize": {"istate": 0, "jstate": 0},
            "pcm": {"enabled": False},
        }

        report = input_checker.CheckReport()
        input_checker._check_dftb(config, report)

        assert not report.ok
        assert "response type conflicts" in report.to_text()
        assert expected in report.to_text()


def test_sf_reference_multiplicity_one_normalizes_to_triplet(monkeypatch):
    calls = []
    module = _load_adapter(monkeypatch, calls)
    adapter = _adapter_without_native_runtime(module)
    adapter.dftb["reference_multiplicity"] = 1

    assert adapter._reference_multiplicity("sf") == 3
    assert adapter._reference_multiplicity("mrsf") == 3
