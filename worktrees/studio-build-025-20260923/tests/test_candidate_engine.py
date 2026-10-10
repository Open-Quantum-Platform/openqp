import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

spec = importlib.util.spec_from_file_location(
    "ci_package", Path(__file__).parents[1] / "tools/release/ci_package.py")
candidate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(candidate)


class QualifiedEngineReuseTests(unittest.TestCase):
    def run_candidate(self, expected):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root / "engine.zip"
            with zipfile.ZipFile(archive, "w") as stream:
                stream.writestr("engine-identity.json", '{"commit":"engine-sha"}')
                stream.writestr("openqp.exe", "qualified binary")
            digest = hashlib.sha256(archive.read_bytes()).hexdigest()
            (root / "release-candidate.json").write_text(json.dumps({"engine": {"commit": "engine-sha"}}))
            env = {"STUDIO_RELEASE_CANDIDATE": "1", "STUDIO_PACKAGE_PLATFORM": "windows-x64",
                   "STUDIO_ENGINE_ARCHIVE": str(archive), "STUDIO_ENGINE_SHA256": digest if expected else "0" * 64}
            with patch.object(candidate, "ROOT", root), patch.dict(os.environ, env, clear=True), \
                 patch.object(candidate.platform, "system", return_value="Windows"), \
                 patch.object(candidate.venv.EnvBuilder, "create"), \
                 patch.object(candidate.subprocess, "run") as run:
                if expected:
                    candidate.main()
                    engine = root / ".cache/qualified-engine/openqp"
                    self.assertEqual((engine / "openqp.exe").read_text(), "qualified binary")
                    self.assertIn(str(engine), run.call_args.args[0])
                    self.assertTrue(run.call_args.kwargs["check"])
                else:
                    with self.assertRaisesRegex(SystemExit, "checksum mismatch"):
                        candidate.main()
                    run.assert_not_called()
                    self.assertFalse((root / ".cache/qualified-engine").exists())

    def test_verified_local_windows_archive_is_requalified_by_packager(self):
        self.run_candidate(True)

    def test_changed_archive_is_rejected_before_extracting_or_building(self):
        self.run_candidate(False)
