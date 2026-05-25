import importlib.util
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def install_oqpdata_stubs():
    numpy = types.ModuleType("numpy")
    for name in (
        "void",
        "int8",
        "int16",
        "int32",
        "int64",
        "uint8",
        "uint16",
        "uint32",
        "uint64",
        "float32",
        "float64",
        "complex64",
        "complex128",
    ):
        setattr(numpy, name, type(name, (), {}))
    setattr(numpy, "dtype", lambda value: value)
    sys.modules["numpy"] = numpy

    oqp = types.ModuleType("oqp")
    setattr(oqp, "ffi", object())
    setattr(oqp, "lib", object())
    periodic_table = types.ModuleType("oqp.periodic_table")
    setattr(periodic_table, "MASSES", {})
    setattr(periodic_table, "SYMBOL_MAP", {})
    constants = types.ModuleType("oqp.utils.constants")
    setattr(constants, "ANGSTROM_TO_BOHR", 1.8897261254576558)
    sys.modules["oqp"] = oqp
    sys.modules["oqp.periodic_table"] = periodic_table
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    sys.modules["oqp.utils.constants"] = constants


def install_input_checker_stubs():
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")
    setattr(mpi_utils, "MPIManager", lambda: types.SimpleNamespace(use_mpi=False))
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils


class GpuMetcConfigTests(unittest.TestCase):
    def test_gpu_section_is_in_default_config_schema(self):
        install_oqpdata_stubs()
        oqpdata = load_module("oqpdata_under_test", "pyoqp/oqp/molecule/oqpdata.py")

        gpu_schema = oqpdata.OQP_CONFIG_SCHEMA["gpu"]

        self.assertEqual(gpu_schema["enabled"]["default"], "False")
        self.assertEqual(gpu_schema["backend"]["default"], "cuda")
        self.assertEqual(gpu_schema["target"]["default"], "metc")
        self.assertEqual(gpu_schema["device"]["default"], "0")
        self.assertEqual(gpu_schema["precision"]["default"], "float64")
        self.assertEqual(gpu_schema["fallback"]["default"], "cpu")

    def test_input_checker_rejects_unknown_gpu_backend(self):
        install_input_checker_stubs()
        checker = load_module("input_checker_under_test", "pyoqp/oqp/utils/input_checker.py")

        config = {
            "input": {
                "system": "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7",
                "method": "tdhf",
                "runtype": "energy",
                "basis": "sto-3g",
            },
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "nstate": 3},
            "gpu": {"enabled": True, "backend": "metal", "target": "metc", "precision": "float64", "fallback": "cpu"},
        }

        report = checker.check_input_values(config, raise_error=False, emit=False)

        self.assertFalse(report.ok)
        self.assertTrue(any(item.path == "gpu.backend" for item in report.errors))

    def test_input_checker_accepts_cuda_metc_gpu_target_for_mrsf(self):
        install_input_checker_stubs()
        checker = load_module("input_checker_under_test_accept", "pyoqp/oqp/utils/input_checker.py")

        config = {
            "input": {
                "system": "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7",
                "method": "tdhf",
                "runtype": "energy",
                "basis": "sto-3g",
            },
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "nstate": 3},
            "gpu": {"enabled": True, "backend": "cuda", "target": "metc", "precision": "float64", "fallback": "cpu"},
        }

        report = checker.check_input_values(config, raise_error=False, emit=False)

        self.assertTrue(report.ok, report.to_text())
    def test_runtime_gpu_config_helper_marks_mrsf_metc_as_supported(self):
        gpu = load_module("gpu_helper_under_test", "pyoqp/oqp/utils/gpu.py")

        runtime = gpu.GpuConfig.from_config({
            "input": {"method": "tdhf"},
            "tdhf": {"type": "mrsf"},
            "gpu": {"enabled": True, "backend": "cuda", "target": "metc", "device": 1},
        })

        self.assertTrue(runtime.enabled)
        self.assertEqual(runtime.device, 1)
        self.assertTrue(runtime.supports_metc({"input": {"method": "tdhf"}, "tdhf": {"type": "mrsf"}}))
    def test_cmake_exposes_optional_cuda_build_switch(self):
        root_cmake = (ROOT / "CMakeLists.txt").read_text()
        source_cmake = (ROOT / "source/CMakeLists.txt").read_text()

        self.assertIn("option(ENABLE_CUDA", root_cmake)
        self.assertIn("option(CUDA_ALLOW_UNSUPPORTED_COMPILER", root_cmake)
        self.assertIn("-allow-unsupported-compiler", root_cmake)
        self.assertIn("find_package(CUDAToolkit REQUIRED)", root_cmake)
        self.assertIn("CUDA::cudart", source_cmake)
        self.assertIn("CUDA::cublas", source_cmake)
        self.assertIn("gpu_metc_cuda.cu", source_cmake)

    def test_fortran_gpu_backend_stub_is_present(self):
        source = (ROOT / "source/gpu_backend.F90").read_text()

        self.assertIn("module gpu_backend", source)
        self.assertIn("gpu_backend_available", source)
        self.assertIn("gpu_backend_describe", source)
        self.assertIn("integer(2), intent(in) :: ids", source)
        self.assertIn("gpu_backend_metc_contract", source)
        self.assertIn("oqp_gpu_metc_contract", source)

    def test_cuda_metc_kernel_exports_real_contract_abi(self):
        source = (ROOT / "source/gpu_metc_cuda.cu").read_text()

        self.assertIn("extern \"C\" int oqp_gpu_metc_contract", source)
        self.assertIn("__global__", source)
        self.assertIn("atomicAdd", source)
        self.assertIn("mrsf_metc_kernel", source)
        self.assertIn("umrsf_metc_kernel", source)

    def test_mrsf_update_dispatches_to_gpu_backend_before_cpu_fallback(self):
        source = (ROOT / "source/tdhf_mrsf_lib.F90").read_text()

        self.assertIn("use gpu_backend, only: gpu_backend_metc_enabled", source)
        self.assertIn("call gpu_backend_metc_contract", source)
        self.assertIn("is_umrsf=.false.", source)
        self.assertIn("is_umrsf=.true.", source)


if __name__ == "__main__":
    unittest.main()
