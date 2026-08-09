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

    def test_docker_image_is_versioned_but_not_automatically_pushed(self):
        workflow = (ROOT / ".github" / "workflows" / "docker-build.yml").read_text()
        dockerfile = (ROOT / "Dockerfile").read_text()
        pyproject = (ROOT / "pyproject.toml").read_text()
        version = find_version(r'^version\s*=\s*"([^"]+)"', pyproject, "pyproject.toml")

        self.assertEqual(dockerfile.count(f"ARG OPENQP_VERSION={version}"), 2)
        self.assertIn('expected_tag = f"v{version}"', workflow)
        self.assertIn("tags: openqp/openqp:${{ steps.source.outputs.docker_tag }}", workflow)
        self.assertIn("OPENQP_VERSION=${{ steps.source.outputs.version }}", workflow)
        self.assertIn("OPENQP_REVISION=${{ steps.source.outputs.revision }}", workflow)
        self.assertIn("outputs: type=oci", workflow)
        self.assertIn("sbom: true", workflow)
        self.assertIn("provenance: mode=max", workflow)
        self.assertIn("push: false", workflow)
        self.assertNotIn("docker push ", workflow)
        self.assertNotIn("push: true", workflow)
        self.assertNotIn("load: true", workflow)

    def test_pypi_publication_requires_manual_protected_job(self):
        workflow = (ROOT / ".github" / "workflows" / "build_wheels.yml").read_text()

        self.assertIn("inputs.action == 'publish_pypi'", workflow)
        self.assertIn("name: pypi", workflow)
        self.assertIn("LicenseRef-OpenQP-Research-1.0", workflow)
        self.assertNotIn("publish_existing_artifacts", workflow)

    def test_pypi_sealed_artifact_survives_failed_job_reruns(self):
        workflow = (ROOT / ".github" / "workflows" / "build_wheels.yml").read_text()

        artifact_name = "verified-openqp-pypi-${{"
        artifact_lines = [line for line in workflow.splitlines() if artifact_name in line]
        self.assertEqual(len(artifact_lines), 2)
        self.assertTrue(all("github.run_id" in line for line in artifact_lines))
        self.assertTrue(all("github.run_attempt" not in line for line in artifact_lines))
