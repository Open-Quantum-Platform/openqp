"""Exercise resolution and the public-main boundary using real temporary Git repos."""
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "checkout_engine", Path(__file__).parents[1] / "tools/checkout_engine.py"
)
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)


class EngineSourceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "source"
        self.source.mkdir()
        self.git("init", "-q", "-b", "main")
        self.git("config", "user.name", "Test")
        self.git("config", "user.email", "test@example.invalid")
        self.git("commit", "-q", "--allow-empty", "-m", "public main")
        self.public = self.git("rev-parse", "HEAD")
        self.git("checkout", "-q", "-b", "unmerged")
        self.git("commit", "-q", "--allow-empty", "-m", "not on public main")
        self.git("checkout", "-q", "main")
        self.route = patch.object(engine, "SSH_GATEWAY", str(self.source))
        self.route.start()
        self.addCleanup(self.route.stop)

    def git(self, *args):
        return subprocess.check_output(
            ["git", "-C", str(self.source), *args], text=True
        ).strip()

    def test_main_resolves_to_exact_detached_commit(self):
        target = self.root / "engine"
        manifest = engine.checkout(target, "main", ssh=True)
        self.assertEqual(manifest["commit"], self.public)
        self.assertEqual(manifest["repository"], engine.GATEWAY)
        result = subprocess.run(["git", "-C", str(target), "symbolic-ref", "-q", "HEAD"])
        self.assertNotEqual(result.returncode, 0)

    def test_pinned_commit_does_not_move_when_main_advances(self):
        self.git("commit", "-q", "--allow-empty", "-m", "new main")
        manifest = engine.checkout(self.root / "engine", self.public, ssh=True)
        self.assertEqual(manifest["commit"], self.public)

    def test_rejects_ref_not_in_public_main_history(self):
        with self.assertRaises(subprocess.CalledProcessError):
            engine.checkout(self.root / "engine", "unmerged", ssh=True)

    def test_rejects_revision_expressions_and_options(self):
        for ref in ["--upload-pack=x", "main~1", "main..secret", "main\nsecret"]:
            with self.subTest(ref=ref), self.assertRaises(ValueError):
                engine.checkout(self.root / "engine", ref, ssh=True)

    def test_never_overwrites_existing_checkout(self):
        with self.assertRaises(ValueError):
            engine.checkout(self.source, "main", ssh=True)

    def test_rejects_internal_source_configuration(self):
        config = {"repository": engine.GATEWAY.replace("/openqp.git", "/internal/openqp.git"),
                  "tracking_branch": "main"}
        with patch.object(Path, "read_text", return_value=json.dumps(config)):
            with self.assertRaises(ValueError):
                engine.configuration()

    def test_ci_auth_is_ephemeral_and_scoped_to_gitlab(self):
        with patch.dict("os.environ", {"CI_JOB_TOKEN": "test-only"}, clear=True):
            env = engine.git_env()
        self.assertEqual(env["GIT_CONFIG_COUNT"], "1")
        self.assertEqual(env["GIT_CONFIG_KEY_0"],
                         "http.https://qchemlab.knu.ac.kr/.extraHeader")
        self.assertNotIn("test-only", env["GIT_CONFIG_VALUE_0"])


if __name__ == "__main__":
    unittest.main()
