"""Keep unchanged dependencies in the established external cache namespace."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ExternalCacheCompatibilityTests(unittest.TestCase):
    def test_historical_version_key_still_selects_the_existing_artifacts(self):
        cmake = shutil.which("cmake")
        if cmake is None:
            self.skipTest("CMake is required to evaluate the external cache key")
        source = (ROOT / "external/CMakeLists.txt").read_text()
        # Evaluate the real key declarations without configuring, fetching, or
        # populating any shared cache. DFT-D4 is a separately versioned subkey.
        prefix = source.split("# Give only the DFT-D4 stack", 1)[0]
        prefix = prefix.replace("include(ExternalProject)", "")
        expected = ("ext2-libint2.7.1.1-am4-nlopt2.9.1-libxc7.0.0-tag1.0.0-"
                    "ecp1.0.7-lapack3.10.0-otr2.0.0-fmt0.3.7-ddx0.8.0")
        with tempfile.TemporaryDirectory() as directory:
            script = Path(directory) / "cache-key.cmake"
            script.write_text(prefix + '\nmessage("CACHE_KEY=${_OQP_EXTERNALS_VERSION_KEY}")\n')
            result = subprocess.run([cmake, "-P", str(script)],
                                    capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("CACHE_KEY=" + expected, result.stderr)
