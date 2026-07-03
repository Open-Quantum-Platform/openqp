import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def find_version(pattern, text, source):
    match = re.search(pattern, text, re.MULTILINE)
    if not match:
        raise AssertionError(f"Could not find version in {source}")
    return match.group(1)


class ReleaseMetadataTests(unittest.TestCase):
    def test_release_versions_stay_in_sync(self):
        pyproject = (ROOT / "pyproject.toml").read_text()
        cmake = (ROOT / "CMakeLists.txt").read_text()
        setup_py = (ROOT / "pyoqp" / "setup.py").read_text()

        version = find_version(r'^version\s*=\s*"([^"]+)"', pyproject, "pyproject.toml")
        self.assertEqual(
            find_version(r"VERSION\s+([0-9]+\.[0-9]+\.[0-9]+)", cmake, "CMakeLists.txt"),
            version,
        )
        self.assertEqual(
            find_version(r'version="([^"]+)"', setup_py, "pyoqp/setup.py"),
            version,
        )

    def test_docker_image_tag_comes_from_pyproject_version(self):
        workflow = (ROOT / ".github" / "workflows" / "docker-build.yml").read_text()

        self.assertIn("docker_tag=v{version}", workflow)
        self.assertIn("tags: openqp/openqp:${{ steps.version.outputs.docker_tag }}", workflow)
        self.assertIn("docker push openqp/openqp:${{ steps.version.outputs.docker_tag }}", workflow)
