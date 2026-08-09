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
            find_version(
                r'^OPENQP_VERSION\s*=\s*"([^"]+)"', setup_py, "pyoqp/setup.py"
            ),
            version,
        )

    def test_legacy_pyoqp_distribution_is_retired(self):
        setup_py = (ROOT / "pyoqp" / "setup.py").read_text()
        self.assertIn("pyoqp/setup.py is retired", setup_py)
        self.assertNotIn("setup(", setup_py)

    def test_docker_image_is_versioned_but_not_automatically_pushed(self):
        workflow = (ROOT / ".github" / "workflows" / "docker-build.yml").read_text()

        self.assertIn("docker_tag=v{version}", workflow)
        self.assertIn("tags: openqp/openqp:${{ steps.version.outputs.docker_tag }}", workflow)
        self.assertIn("OPENQP_VERSION=${{ steps.version.outputs.version }}", workflow)
        self.assertIn("OPENQP_REVISION=${{ github.sha }}", workflow)
        self.assertIn("load: true", workflow)
        self.assertNotIn("docker push ", workflow)
        self.assertNotIn("push: true", workflow)

    def test_pypi_publication_requires_manual_protected_job(self):
        workflow = (ROOT / ".github" / "workflows" / "build_wheels.yml").read_text()

        self.assertIn("inputs.action == 'publish_pypi'", workflow)
        self.assertIn("name: pypi", workflow)
        self.assertIn("LicenseRef-OpenQP-Research-1.0", workflow)
        self.assertNotIn("publish_existing_artifacts", workflow)

    def test_packaged_policy_files_trigger_distribution_ci(self):
        workflow = (ROOT / ".github" / "workflows" / "build_wheels.yml").read_text()
        self.assertEqual(workflow.count('- "CONTRIBUTING.md"'), 2)
        self.assertEqual(workflow.count('- "SUSTAINABILITY.md"'), 2)

    def test_publication_requires_the_complete_configured_wheel_matrix(self):
        workflow = (ROOT / ".github" / "workflows" / "build_wheels.yml").read_text()
        self.assertEqual(workflow.count("Incomplete configured wheel matrix"), 2)
        for tag in ("cp39", "cp310", "cp311", "cp312", "cp313", "cp314"):
            self.assertGreaterEqual(workflow.count(f'"{tag}"'), 2)
        for platform in (
            "linux-x86_64", "linux-aarch64", "macos-x86_64", "macos-arm64"
        ):
            self.assertGreaterEqual(workflow.count(f'"{platform}"'), 2)

    def test_release_is_revalidated_after_environment_approval(self):
        workflow = (ROOT / ".github" / "workflows" / "build_wheels.yml").read_text()
        revalidate = workflow.index(
            "Revalidate release tag and GitHub Release after approval"
        )
        upload = workflow.index("Attach approved artifacts to GitHub Release")
        publish = workflow.index("Publish to PyPI")
        self.assertLess(revalidate, upload)
        self.assertLess(revalidate, publish)
        self.assertIn("manual_pypi_preflight.outputs.release_sha", workflow)
        self.assertIn("Release tag moved after preflight approval", workflow)

    def test_pypi_sealed_artifact_survives_failed_job_reruns(self):
        workflow = (ROOT / ".github" / "workflows" / "build_wheels.yml").read_text()

        artifact_name = "verified-openqp-pypi-${{"
        artifact_lines = [line for line in workflow.splitlines() if artifact_name in line]
        self.assertEqual(len(artifact_lines), 2)
        self.assertTrue(all("github.run_id" in line for line in artifact_lines))
        self.assertTrue(all("github.run_attempt" not in line for line in artifact_lines))
