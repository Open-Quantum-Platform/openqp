import hashlib
import importlib.util
import io
import json
import os
import tarfile
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]


def load_script(name):
    path = ROOT / ".github" / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


COLLECTOR = load_script("collect_docker_runtime")
WHEELHOUSE = load_script("record_wheelhouse")
OCI_VERIFY = load_script("verify_oci_artifact")


class DockerDistributionTests(unittest.TestCase):
    def test_base_images_are_digest_locked_and_runtime_is_minimal(self):
        lock = json.loads((ROOT / "docker/base-images.lock.json").read_text())
        dockerfile = (ROOT / "Dockerfile").read_text()
        from_lines = [
            line for line in dockerfile.splitlines() if line.startswith("FROM ")
        ]
        self.assertEqual(len(from_lines), 2)
        for stage in ("builder", "runtime"):
            image = lock["images"][stage]
            self.assertRegex(image["digest"], r"^sha256:[0-9a-f]{64}$")
            self.assertEqual(image["platform"], "linux/amd64")
            self.assertIn(
                f"FROM {image['reference']}@{image['digest']} AS {stage}",
                dockerfile,
            )
        runtime = dockerfile.split(" AS runtime", 1)[1]
        self.assertNotIn("apt-get", runtime)
        self.assertIn("--no-cache-dir --no-index", runtime)
        self.assertIn("source=/opt/openqp-wheelhouse", runtime)
        self.assertIn("USER 65532:65532", runtime)
        self.assertIn("wheelhouse-manifest.json", runtime)
        self.assertIn("runtime-library-manifest.json", runtime)
        self.assertIn("-DUSE_LIBINT=OFF", dockerfile)
        self.assertIn("-DENABLE_OPENMP=ON", dockerfile)
        self.assertIn("--base-lock=docker/base-images.lock.json", dockerfile)
        self.assertIn(
            "OQP_WHEEL_SMOKE_EXTERNAL_RUNTIME_PATH=/opt/openblas/lib",
            dockerfile,
        )
        wheel_smoke = (
            ROOT / ".github" / "scripts" / "wheel_smoke_test.py"
        ).read_text()
        self.assertIn("OQP_WHEEL_SMOKE_EXTERNAL_RUNTIME_PATH", wheel_smoke)
        self.assertIn("external runtime directory contains a DFT-D4 library", wheel_smoke)
        self.assertEqual(
            lock["images"]["builder"]["openblas_abi"],
            {
                "dynamic_arch": True,
                "interface64": True,
                "openmp": True,
                "symbol_suffix": "none",
            },
        )

    def test_workflow_is_ephemeral_build_only_with_attestations(self):
        workflow = (ROOT / ".github/workflows/docker-build.yml").read_text()
        self.assertIn("outputs: type=oci", workflow)
        self.assertIn("sbom: true", workflow)
        self.assertIn("provenance: mode=max", workflow)
        self.assertIn("push: false", workflow)
        for forbidden in (
            "push: true",
            "load: true",
            "docker/login-action",
            "actions/upload-artifact",
            "cache-to:",
            "actions/cache",
            "packages: write",
        ):
            self.assertNotIn(forbidden, workflow)
        self.assertIn("verify_oci_artifact.py", workflow)
        self.assertIn("docker/base-images.lock.json", workflow)
        self.assertIn('DOCKER_BUILD_RECORD_UPLOAD: "false"', workflow)
        for action_sha in (
            "actions/checkout@11d5960a326750d5838078e36cf38b85af677262 # v4",
            "docker/setup-buildx-action@8d2750c68a42422c14e847fe6c8ac0403b4cbd6f # v3",
            "docker/build-push-action@10e90e3645eae34f1e60eeb005ba3a3d33f178e8 # v6",
        ):
            self.assertIn(action_sha, workflow)

    def test_final_runtime_gate_is_fail_closed(self):
        smoke = (
            ROOT / ".github/scripts/container_runtime_smoke.py"
        ).read_text()
        for required in (
            "container-smoke-hidden",
            "loaded alternate libraries",
            "NLopt dependency leaked",
            "NLopt symbols/strings leaked",
            "still describes NLopt as shipped",
            "still describes static DFT-D4",
            "static archives leaked",
            "build tool leaked",
            "corresponding source missing",
            "DFT-D4 BUILD-INFO schema/runtime names are invalid",
            "DFT-D4 patch record is invalid",
            'publication_gate", {}).get("ready") is not False',
        ):
            self.assertIn(required, smoke)

    def test_ldd_parser_rejects_unresolved_and_parses_absolute_paths(self):
        output = """
linux-vdso.so.1 (0x0000)
libgfortran.so.5 => /lib/x86_64-linux-gnu/libgfortran.so.5 (0x1234)
/lib64/ld-linux-x86-64.so.2 (0x5678)
"""
        self.assertEqual(
            COLLECTOR.parse_ldd(output),
            [
                (
                    "libgfortran.so.5",
                    Path("/lib/x86_64-linux-gnu/libgfortran.so.5"),
                ),
                (
                    "ld-linux-x86-64.so.2",
                    Path("/lib64/ld-linux-x86-64.so.2"),
                ),
            ],
        )
        with self.assertRaisesRegex(RuntimeError, "unresolved"):
            COLLECTOR.parse_ldd("libgomp.so.1 => not found\n")

    def test_ldd_uses_only_explicit_package_runtime_paths(self):
        completed = mock.Mock(returncode=0, stdout="", stderr="")
        with mock.patch.object(COLLECTOR.subprocess, "run", return_value=completed) as run:
            COLLECTOR.ldd(
                Path("/tmp/owner.so"),
                (Path("/package/oqp/lib"), Path("/opt/openblas/lib")),
            )

        environment = run.call_args.kwargs["env"]
        self.assertEqual(
            environment["LD_LIBRARY_PATH"],
            os.pathsep.join(("/package/oqp/lib", "/opt/openblas/lib")),
        )

    @staticmethod
    def _write_wheel(path, name, version, tag="py3-none-any"):
        dist = name.replace("-", "_")
        with zipfile.ZipFile(path, "w") as wheel:
            wheel.writestr(
                f"{dist}-{version}.dist-info/METADATA",
                f"Metadata-Version: 2.2\nName: {name}\nVersion: {version}\n",
            )
            wheel.writestr(
                f"{dist}-{version}.dist-info/WHEEL",
                f"Wheel-Version: 1.0\nTag: {tag}\n",
            )

    def test_wheelhouse_manifest_records_exact_hashes_and_stays_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            wheelhouse = Path(temporary)
            self._write_wheel(wheelhouse / "OpenQP-1.3.0-py3-none-any.whl", "OpenQP", "1.3.0")
            self._write_wheel(wheelhouse / "numpy-2.0-py3-none-any.whl", "numpy", "2.0")
            manifest = WHEELHOUSE.record(wheelhouse, "1.3.0")
        self.assertIs(manifest["publication_gate"]["ready"], False)
        self.assertEqual(len(manifest["wheels"]), 2)
        for wheel in manifest["wheels"]:
            self.assertRegex(wheel["sha256"], r"^[0-9a-f]{64}$")

    @staticmethod
    def _json_blob(blobs, value, media_type="application/vnd.oci.image.manifest.v1+json"):
        payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
        digest = hashlib.sha256(payload).hexdigest()
        blobs[f"blobs/sha256/{digest}"] = payload
        return {
            "mediaType": media_type,
            "digest": f"sha256:{digest}",
            "size": len(payload),
        }

    @staticmethod
    def _write_oci(path, index, blobs):
        layout = b'{"imageLayoutVersion":"1.0.0"}'
        with tarfile.open(path, "w") as archive:
            for name, payload in {
                "index.json": index,
                "oci-layout": layout,
                **blobs,
            }.items():
                info = tarfile.TarInfo(name)
                info.size = len(payload)
                archive.addfile(info, io.BytesIO(payload))

    def test_oci_verifier_requires_labels_sbom_and_provenance(self):
        blobs = {}
        config = self._json_blob(
            blobs,
            {
                "architecture": "amd64",
                "os": "linux",
                "config": {
                    "Labels": {
                        "org.opencontainers.image.source": "https://github.com/Open-Quantum-Platform/openqp",
                        "org.opencontainers.image.version": "1.3.0",
                        "org.opencontainers.image.revision": "a" * 40,
                        "org.opencontainers.image.licenses": "LicenseRef-OpenQP-Research-1.0",
                    }
                },
            },
            "application/vnd.oci.image.config.v1+json",
        )
        image = self._json_blob(
            blobs,
            {"schemaVersion": 2, "config": config, "layers": []},
        )
        image["platform"] = {"os": "linux", "architecture": "amd64"}
        image_digest = image["digest"]
        statements = []
        for predicate_type in (
            "https://spdx.dev/Document",
            "https://slsa.dev/provenance/v1",
        ):
            statements.append(
                self._json_blob(
                    blobs,
                    {
                        "_type": "https://in-toto.io/Statement/v1",
                        "subject": [
                            {
                                "name": "openqp/openqp",
                                "digest": {"sha256": image_digest.split(":", 1)[1]},
                            }
                        ],
                        "predicateType": predicate_type,
                        "predicate": {},
                    },
                    "application/vnd.in-toto+json",
                )
            )
        attestation = self._json_blob(
            blobs,
            {
                "schemaVersion": 2,
                "subject": {"digest": image_digest},
                "layers": statements,
            },
        )
        attestation["platform"] = {"os": "unknown", "architecture": "unknown"}
        attestation["annotations"] = {
            "vnd.docker.reference.type": "attestation-manifest",
            "vnd.docker.reference.digest": image_digest,
        }
        index = json.dumps(
            {"schemaVersion": 2, "manifests": [image, attestation]},
            sort_keys=True,
        ).encode()
        with tempfile.TemporaryDirectory() as temporary:
            artifact = Path(temporary) / "candidate.oci.tar"
            self._write_oci(artifact, index, blobs)
            summary = OCI_VERIFY.verify(artifact, "1.3.0", "a" * 40)

            tampered_blobs = dict(blobs)
            config_path = OCI_VERIFY._blob_path(config["digest"])
            tampered_blobs[config_path] = tampered_blobs[config_path].replace(
                b'"amd64"', b'"arm64"'
            )
            tampered = Path(temporary) / "tampered.oci.tar"
            self._write_oci(tampered, index, tampered_blobs)
            with self.assertRaisesRegex(ValueError, "descriptor digest mismatch"):
                OCI_VERIFY.verify(tampered, "1.3.0", "a" * 40)

            unrelated_spdx = self._json_blob(
                blobs,
                {
                    "_type": "https://in-toto.io/Statement/v1",
                    "subject": [
                        {
                            "name": "unrelated/image",
                            "digest": {"sha256": "b" * 64},
                        }
                    ],
                    "predicateType": "https://spdx.dev/Document",
                    "predicate": {},
                },
                "application/vnd.in-toto+json",
            )
            unrelated_attestation = self._json_blob(
                blobs,
                {
                    "schemaVersion": 2,
                    "subject": {"digest": image_digest},
                    "layers": [unrelated_spdx, statements[1]],
                },
            )
            unrelated_attestation["platform"] = {
                "os": "unknown", "architecture": "unknown"
            }
            unrelated_attestation["annotations"] = dict(attestation["annotations"])
            unrelated_index = json.dumps(
                {"schemaVersion": 2, "manifests": [image, unrelated_attestation]},
                sort_keys=True,
            ).encode()
            unrelated = Path(temporary) / "unrelated-sbom.oci.tar"
            self._write_oci(unrelated, unrelated_index, blobs)
            with self.assertRaisesRegex(ValueError, "candidate image has no SPDX SBOM"):
                OCI_VERIFY.verify(unrelated, "1.3.0", "a" * 40)
        self.assertEqual(summary["image_digest"], image_digest)
        self.assertEqual(summary["platform"], "linux/amd64")


if __name__ == "__main__":
    unittest.main()
