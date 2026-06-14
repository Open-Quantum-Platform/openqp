import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT / "pyoqp" / "oqp" / "runtime.py"


def load_runtime_module():
    spec = importlib.util.spec_from_file_location("openqp_runtime_test", RUNTIME)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_runtime_root(root, suffix):
    (root / "include").mkdir(parents=True)
    (root / "lib").mkdir()
    (root / "share" / "basis_sets").mkdir(parents=True)
    (root / "include" / "oqp.h").write_text("")
    (root / "lib" / f"liboqp.{suffix}").write_text("")


class RuntimeRootResolutionTests(unittest.TestCase):
    def test_python_runtime_prefers_package_local_root_over_openqp_root(self):
        runtime = load_runtime_module()
        suffix = runtime.library_suffix()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            package_root = tmp_path / "site-packages" / "oqp"
            env_root = tmp_path / "other-openqp"
            make_runtime_root(package_root, suffix)
            make_runtime_root(env_root, suffix)

            with patch.dict(os.environ, {"OPENQP_ROOT": str(env_root)}, clear=False):
                root, actual_suffix = runtime.resolve_oqp_root(package_root=package_root)

                self.assertEqual(root, str(package_root.resolve()))
                self.assertEqual(actual_suffix, suffix)
                self.assertEqual(os.environ["OPENQP_ROOT"], str(env_root))

    def test_python_runtime_detects_source_tree_without_openqp_root(self):
        runtime = load_runtime_module()
        suffix = runtime.library_suffix()
        env = os.environ.copy()
        env.pop("OPENQP_ROOT", None)
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            source_root = tmp_path / "openqp"
            package_root = source_root / "pyoqp" / "oqp"
            package_root.mkdir(parents=True)
            make_runtime_root(source_root, suffix)

            with patch.dict(os.environ, env, clear=True):
                root, actual_suffix = runtime.resolve_oqp_root(package_root=package_root)

                self.assertEqual(root, str(source_root.resolve()))
                self.assertEqual(actual_suffix, suffix)

    def test_python_runtime_keeps_openqp_root_as_fallback(self):
        runtime = load_runtime_module()
        suffix = runtime.library_suffix()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            package_root = tmp_path / "isolated" / "oqp"
            env_root = tmp_path / "configured-openqp"
            package_root.mkdir(parents=True)
            make_runtime_root(env_root, suffix)

            with patch.dict(os.environ, {"OPENQP_ROOT": str(env_root)}, clear=False):
                root, actual_suffix = runtime.resolve_oqp_root(package_root=package_root)

                self.assertEqual(root, str(env_root.resolve()))
                self.assertEqual(actual_suffix, suffix)

    def test_non_rtld_import_resolves_package_root_without_mutating_env(self):
        runtime = load_runtime_module()
        suffix = runtime.library_suffix()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            site_root = tmp_path / "site"
            package_root = site_root / "oqp"
            (package_root / "utils").mkdir(parents=True)
            (package_root / "__init__.py").write_text(
                (ROOT / "pyoqp" / "oqp" / "__init__.py").read_text()
            )
            (package_root / "runtime.py").write_text(RUNTIME.read_text())
            (package_root / "utils" / "__init__.py").write_text("")
            (package_root / "utils" / "mpi_utils.py").write_text(
                "class MPIManager:\n"
                "    def __init__(self):\n"
                "        pass\n"
            )
            (site_root / "_oqp.py").write_text(
                "ffi = object()\n"
                "class Lib:\n"
                "    def __dir__(self):\n"
                "        return []\n"
                "lib = Lib()\n"
            )
            make_runtime_root(package_root, suffix)

            env_root = tmp_path / "other-openqp"
            env = os.environ.copy()
            env.update({
                "OPENQP_ROOT": str(env_root),
                "OQP_RTLD": "0",
                "PYTHONPATH": str(site_root),
            })
            result = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    "import os, oqp; "
                    "print(oqp.oqp_root); "
                    "print(os.environ.get('OPENQP_ROOT'))",
                ],
                check=False,
                capture_output=True,
                env=env,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            lines = result.stdout.strip().splitlines()
            self.assertEqual(lines[-2], str(package_root.resolve()))
            self.assertEqual(lines[-1], str(env_root))

    def test_cmake_forces_ilp64_when_stale_cache_sets_lp64(self):
        source = (ROOT / "CMakeLists.txt").read_text()

        self.assertIn("Forcing LINALG_LIB_INT64=ON", source)
        self.assertIn("set(LINALG_LIB_INT64 ON CACHE BOOL", source)
        self.assertIn("FORCE)", source)
