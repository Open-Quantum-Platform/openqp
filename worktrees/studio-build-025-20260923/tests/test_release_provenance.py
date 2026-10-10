"""A packaging-only revision must never relabel changed application code."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools/release"))
from provenance import source_commits


class ProvenanceTests(unittest.TestCase):
    def test_packaging_change_preserves_source_but_app_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            def git(*args):
                return subprocess.check_output(["git", *args], cwd=root, text=True).strip()

            git("init", "-q")
            git("config", "user.name", "Test")
            git("config", "user.email", "test@example.invalid")
            (root / "backend").mkdir()
            app = root / "backend/app.py"
            app.write_text("version = '0.2.4'\n")
            git("add", ".")
            git("commit", "-qm", "application")
            application = git("rev-parse", "HEAD")
            (root / "tools/release").mkdir(parents=True)
            (root / "tools/release/package.py").write_text("bundles = 'nsis'\n")
            git("add", ".")
            git("commit", "-qm", "packaging")
            self.assertEqual(source_commits(root, application),
                             (application, git("rev-parse", "HEAD")))
            app.write_text("version = '0.2.5'\n")
            with self.assertRaisesRegex(ValueError, "commit and review"):
                source_commits(root, application)
            git("add", ".")
            git("commit", "-qm", "different application")
            with self.assertRaisesRegex(ValueError, "application files changed"):
                source_commits(root, application)


if __name__ == "__main__":
    unittest.main()
