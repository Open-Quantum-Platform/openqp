import hashlib
import importlib.util
import io
import json
import os
import struct
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
ELF_STACK = load_script("normalize_elf_stack")
CONTAINER_SMOKE = load_script("container_runtime_smoke")


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
        self.assertIn("normalize_elf_stack.py", runtime)
        self.assertIn("--openqp-package", runtime)
        self.assertNotIn("--update-record", runtime)
        stack_gate = (
            ROOT / ".github/scripts/normalize_elf_stack.py"
        ).read_text()
        self.assertNotIn("write_bytes", stack_gate)
        self.assertNotIn("-ftrampoline-impl=heap", (ROOT / "CMakeLists.txt").read_text())
        runtime_notice = (ROOT / "docker/THIRD_PARTY_RUNTIME.md").read_text()
        self.assertIn("rejects any installed OpenQP ELF library", runtime_notice)
        self.assertIn("fail-closed", runtime_notice)
        self.assertIn("wheel `RECORD` metadata", runtime_notice)
        self.assertNotIn("executable-stack\nnormalization", runtime_notice)
        self.assertIn("-DUSE_LIBINT=OFF", dockerfile)
        self.assertIn("-DENABLE_OPENMP=ON", dockerfile)
        self.assertIn("--base-lock=docker/base-images.lock.json", dockerfile)
        self.assertIn("-DOQP_SOURCE_REVISION=${OPENQP_REVISION}", dockerfile)
        self.assertIn('OPENQP_EXPECTED_REVISION="${OPENQP_REVISION}"', dockerfile)
        self.assertEqual(dockerfile.count("ARG OPENQP_REVISION\n"), 2)
        self.assertNotIn("OPENQP_REVISION=unknown", dockerfile)
        self.assertIn("/opt/openqp-build-wheelhouse", dockerfile)
        self.assertIn("--build-wheelhouse=/opt/openqp-build-wheelhouse", dockerfile)
        self.assertIn("--no-index", dockerfile.split("COPY . /opt/openqp", 1)[0])
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

    def test_elf_stack_gate_rejects_execute_permission_without_mutation(self):
        for elf_class, header_size, phentsize, flags_offset in (
            (1, 52, 32, 24),
            (2, 64, 56, 4),
        ):
            data = bytearray(header_size + phentsize)
            data[:6] = b"\x7fELF" + bytes((elf_class, 1))
            if elf_class == 1:
                struct.pack_into("<I", data, 28, header_size)
                struct.pack_into("<H", data, 42, phentsize)
                struct.pack_into("<H", data, 44, 1)
            else:
                struct.pack_into("<Q", data, 32, header_size)
                struct.pack_into("<H", data, 54, phentsize)
                struct.pack_into("<H", data, 56, 1)
            struct.pack_into("<I", data, header_size, ELF_STACK.PT_GNU_STACK)
            struct.pack_into("<I", data, header_size + flags_offset, 0x7)

            with self.subTest(elf_class=elf_class), tempfile.TemporaryDirectory() as temporary:
                library = Path(temporary) / "liboqp.so"
                library.write_bytes(data)
                original = library.read_bytes()
                with self.assertRaisesRegex(RuntimeError, "executable stack"):
                    ELF_STACK.require_non_executable_stack(library)
                self.assertEqual(library.read_bytes(), original)

                struct.pack_into("<I", data, header_size + flags_offset, 0x6)
                library.write_bytes(data)
                ELF_STACK.require_non_executable_stack(library)
                self.assertEqual(ELF_STACK.gnu_stack_flags(library.read_bytes()), 0x6)

    def test_mrsf_container_smoke_requires_convergence_and_reference_gradient(self):
        reference = [0.0, 0.0, -0.18, 0.08, -0.08, 0.09, -0.08, 0.08, 0.09]
        gradient = """Z-Vector converged
PyOQP dispersion corrected gradients
PyOQP S2
O  0.0  0.0 -0.18
H  0.08 -0.08 0.09
H -0.08  0.08 0.09
"""
        CONTAINER_SMOKE.require_mrsf_minres_result(gradient, reference)
        with self.assertRaisesRegex(AssertionError, "did not converge"):
            CONTAINER_SMOKE.require_mrsf_minres_result(
                gradient.replace("Z-Vector converged", "Z-Vector not converged"),
                reference,
            )
        with self.assertRaisesRegex(AssertionError, "does not match"):
            CONTAINER_SMOKE.require_mrsf_minres_result(
                gradient,
                [*reference[:-1], 1.0],
            )

    def test_container_smoke_reads_direct_elf_needed_entries(self):
        data = bytearray(0x240)
        data[:6] = b"\x7fELF\x02\x01"
        struct.pack_into("<Q", data, 32, 64)
        struct.pack_into("<H", data, 54, 56)
        struct.pack_into("<H", data, 56, 2)

        # One PT_LOAD maps the whole file at 0x400000; PT_DYNAMIC is at 0x100.
        struct.pack_into("<IIQQQQQQ", data, 64, 1, 4, 0, 0x400000, 0,
                         len(data), len(data), 0x1000)
        struct.pack_into("<IIQQQQQQ", data, 120, 2, 4, 0x100, 0x400100, 0,
                         0x60, 0x60, 8)

        strings = b"\0libdftd4.so.3\0libmctc-lib.so.0\0"
        data[0x180:0x180 + len(strings)] = strings
        dynamic_entries = (
            (5, 0x400180),
            (10, len(strings)),
            (1, 1),
            (1, strings.index(b"libmctc-lib.so.0")),
            (0, 0),
        )
        for index, entry in enumerate(dynamic_entries):
            struct.pack_into("<QQ", data, 0x100 + index * 16, *entry)

        with tempfile.TemporaryDirectory() as temporary:
            library = Path(temporary) / "liboqp.so"
            library.write_bytes(data)
            needed = CONTAINER_SMOKE.elf_needed(library)
        self.assertEqual(needed, ["libdftd4.so.3", "libmctc-lib.so.0"])
        CONTAINER_SMOKE.assert_stack_edges(
            Path("liboqp.so"), needed,
            {"libdftd4.so.3": Path("d4"), "libmctc-lib.so.0": Path("mctc")},
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
            build_wheelhouse = wheelhouse / "build"
            build_wheelhouse.mkdir()
            self._write_wheel(wheelhouse / "OpenQP-1.3.0-py3-none-any.whl", "OpenQP", "1.3.0")
            self._write_wheel(wheelhouse / "numpy-2.0-py3-none-any.whl", "numpy", "2.0")
            for name in (
                "scikit-build-core", "cmake", "ninja", "numpy", "cffi", "setuptools"
            ):
                filename = name.replace("-", "_") + "-1.0-py3-none-any.whl"
                self._write_wheel(build_wheelhouse / filename, name, "1.0")
            manifest = WHEELHOUSE.record(wheelhouse, "1.3.0", build_wheelhouse)
        self.assertIs(manifest["publication_gate"]["ready"], False)
        self.assertEqual(manifest["schema_version"], 2)
        self.assertEqual(len(manifest["wheels"]), 2)
        self.assertEqual(len(manifest["build_wheels"]), 6)
        for wheel in [*manifest["wheels"], *manifest["build_wheels"]]:
            self.assertRegex(wheel["sha256"], r"^[0-9a-f]{64}$")

    def test_wheel_metadata_ignores_vendored_dist_info(self):
        with tempfile.TemporaryDirectory() as temporary:
            wheel_path = Path(temporary) / "setuptools-84.0.0-py3-none-any.whl"
            self._write_wheel(wheel_path, "setuptools", "84.0.0")
            with zipfile.ZipFile(wheel_path, "a") as wheel:
                wheel.writestr(
                    "setuptools/_vendor/packaging-26.0.dist-info/METADATA",
                    "Metadata-Version: 2.2\nName: packaging\nVersion: 26.0\n",
                )
                wheel.writestr(
                    "setuptools/_vendor/packaging-26.0.dist-info/WHEEL",
                    "Wheel-Version: 1.0\nTag: py3-none-any\n",
                )

            self.assertEqual(
                WHEELHOUSE.wheel_metadata(wheel_path),
                ("setuptools", "84.0.0", ["py3-none-any"]),
            )

    def test_attestation_predicates_reject_incomplete_v02_and_blank_creator(self):
        subject = [{"name": "image", "digest": {"sha256": "a" * 64}}]
        v02 = {
            "_type": "https://in-toto.io/Statement/v0.1",
            "subject": subject,
            "predicateType": "https://slsa.dev/provenance/v0.2",
            "predicate": {
                "builder": {"id": "https://example.test/builder"},
                "buildType": "https://example.test/build",
                "invocation": {},
                "metadata": {},
                "materials": [],
            },
        }
        OCI_VERIFY._statement_subjects(v02)
        for field in ("invocation", "metadata", "materials"):
            malformed = json.loads(json.dumps(v02))
            malformed["predicate"].pop(field)
            with self.subTest(field=field), self.assertRaisesRegex(
                ValueError, field
            ):
                OCI_VERIFY._statement_subjects(malformed)

        spdx = {
            "_type": "https://in-toto.io/Statement/v1",
            "subject": subject,
            "predicateType": "https://spdx.dev/Document",
            "predicate": {
                "spdxVersion": "SPDX-2.3",
                "SPDXID": "SPDXRef-DOCUMENT",
                "dataLicense": "CC0-1.0",
                "name": "SBOM",
                "documentNamespace": "https://example.test/spdx",
                "creationInfo": {
                    "created": "2026-08-11T00:00:00Z",
                    "creators": ["Tool:  "],
                },
            },
        }
        with self.assertRaisesRegex(ValueError, "creationInfo"):
            OCI_VERIFY._statement_subjects(spdx)

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
    def _write_oci(path, index, blobs, layout_version="1.0.0"):
        entries = {"index.json": index, **blobs}
        if layout_version is not None:
            entries["oci-layout"] = json.dumps(
                {"imageLayoutVersion": layout_version}, separators=(",", ":")
            ).encode()
        with tarfile.open(path, "w") as archive:
            for name, payload in entries.items():
                info = tarfile.TarInfo(name)
                info.size = len(payload)
                archive.addfile(info, io.BytesIO(payload))

    def test_oci_verifier_requires_labels_sbom_and_provenance(self):
        blobs = {}
        layer_buffer = io.BytesIO()
        with tarfile.open(fileobj=layer_buffer, mode="w") as layer_archive:
            payload = b"openqp-rootfs-layer"
            member = tarfile.TarInfo("usr/share/openqp/layer-marker")
            member.size = len(payload)
            layer_archive.addfile(member, io.BytesIO(payload))
        layer_payload = layer_buffer.getvalue()
        layer_digest = hashlib.sha256(layer_payload).hexdigest()
        layer = {
            "mediaType": "application/vnd.oci.image.layer.v1.tar",
            "digest": f"sha256:{layer_digest}",
            "size": len(layer_payload),
        }
        blobs[f"blobs/sha256/{layer_digest}"] = layer_payload
        config = self._json_blob(
            blobs,
            {
                "architecture": "amd64",
                "os": "linux",
                "rootfs": {
                    "type": "layers",
                    "diff_ids": [f"sha256:{layer_digest}"],
                },
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
            {"schemaVersion": 2, "config": config, "layers": [layer]},
        )
        image["platform"] = {"os": "linux", "architecture": "amd64"}
        image_digest = image["digest"]
        image_subject = {
            field: image[field] for field in ("mediaType", "digest", "size")
        }
        statements = []
        for predicate_type in (
            "https://spdx.dev/Document",
            "https://slsa.dev/provenance/v1",
        ):
            if predicate_type == "https://spdx.dev/Document":
                predicate = {
                    "spdxVersion": "SPDX-2.3",
                    "SPDXID": "SPDXRef-DOCUMENT",
                    "dataLicense": "CC0-1.0",
                    "name": "OpenQP container SBOM",
                    "documentNamespace": "https://openqp.dev/spdx/container",
                    "creationInfo": {
                        "created": "2026-08-11T00:00:00Z",
                        "creators": ["Tool: BuildKit"],
                    },
                }
            else:
                predicate = {
                    "buildDefinition": {
                        "buildType": "https://mobyproject.org/buildkit@v1",
                        "externalParameters": {
                            "args": {
                                "build-arg:OPENQP_REVISION": "a" * 40,
                            }
                        },
                    },
                    "runDetails": {
                        "builder": {"id": "https://mobyproject.org/buildkit"}
                    },
                }
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
                        "predicate": predicate,
                    },
                    "application/vnd.in-toto+json",
                )
            )
        attestation_config = self._json_blob(
            blobs,
            {},
            "application/vnd.oci.image.config.v1+json",
        )
        attestation = self._json_blob(
            blobs,
            {
                "schemaVersion": 2,
                "config": attestation_config,
                "subject": dict(image_subject),
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

            nonfinite_index = (
                b'{"schemaVersion":2,"invalid":NaN,"manifests":'
                + json.dumps([image, attestation], sort_keys=True).encode()
                + b"}"
            )
            nonfinite_root = Path(temporary) / "nonfinite-index.oci.tar"
            self._write_oci(nonfinite_root, nonfinite_index, blobs)
            with self.assertRaisesRegex(
                ValueError, "non-standard JSON numeric constant"
            ):
                OCI_VERIFY.verify(nonfinite_root, "1.3.0", "a" * 40)

            config_payload = blobs[OCI_VERIFY._blob_path(config["digest"])]
            nonfinite_config_payload = config_payload[:-1] + b',"invalid":Infinity}'
            nonfinite_config_digest = hashlib.sha256(
                nonfinite_config_payload
            ).hexdigest()
            nonfinite_config = {
                **config,
                "digest": f"sha256:{nonfinite_config_digest}",
                "size": len(nonfinite_config_payload),
            }
            blobs[f"blobs/sha256/{nonfinite_config_digest}"] = (
                nonfinite_config_payload
            )
            nonfinite_image = self._json_blob(
                blobs,
                {
                    "schemaVersion": 2,
                    "config": nonfinite_config,
                    "layers": [layer],
                },
            )
            nonfinite_image["platform"] = image["platform"]
            nonfinite_descriptor = Path(temporary) / "nonfinite-config.oci.tar"
            self._write_oci(
                nonfinite_descriptor,
                json.dumps(
                    {
                        "schemaVersion": 2,
                        "manifests": [nonfinite_image, attestation],
                    },
                    sort_keys=True,
                ).encode(),
                blobs,
            )
            with self.assertRaisesRegex(
                ValueError, "non-standard JSON numeric constant"
            ):
                OCI_VERIFY.verify(nonfinite_descriptor, "1.3.0", "a" * 40)

            duplicate_payloads = (
                ("index.json", index),
                (
                    OCI_VERIFY._blob_path(image["digest"]),
                    blobs[OCI_VERIFY._blob_path(image["digest"])],
                ),
            )
            for duplicate_name, duplicate_payload in duplicate_payloads:
                duplicate = Path(temporary) / (
                    "duplicate-" + duplicate_name.replace("/", "-") + ".oci.tar"
                )
                self._write_oci(duplicate, index, blobs)
                with tarfile.open(duplicate, "a") as archive:
                    info = tarfile.TarInfo(duplicate_name)
                    info.size = len(duplicate_payload)
                    archive.addfile(info, io.BytesIO(duplicate_payload))
                with self.assertRaisesRegex(ValueError, "duplicate member path"):
                    OCI_VERIFY.verify(duplicate, "1.3.0", "a" * 40)

            bad_payload = dict(statements[0])
            bad_payload["mediaType"] = "text/plain"
            bad_payload_attestation = self._json_blob(
                blobs,
                {
                    "schemaVersion": 2,
                    "config": attestation_config,
                    "subject": dict(image_subject),
                    "layers": [bad_payload, statements[1]],
                },
            )
            bad_payload_attestation["platform"] = {
                "os": "unknown",
                "architecture": "unknown",
            }
            bad_payload_attestation["annotations"] = dict(
                attestation["annotations"]
            )
            bad_payload_index = json.dumps(
                {
                    "schemaVersion": 2,
                    "manifests": [image, bad_payload_attestation],
                },
                sort_keys=True,
            ).encode()
            bad_payload_artifact = Path(temporary) / "bad-payload-media.oci.tar"
            self._write_oci(bad_payload_artifact, bad_payload_index, blobs)
            with self.assertRaisesRegex(
                ValueError, "unsupported attestation payload media type"
            ):
                OCI_VERIFY.verify(
                    bad_payload_artifact, "1.3.0", "a" * 40
                )

            wrong_reference_attestation = dict(attestation)
            wrong_reference_attestation["annotations"] = {
                **attestation["annotations"],
                "vnd.docker.reference.digest": "sha256:" + "b" * 64,
            }
            wrong_reference_index = json.dumps(
                {
                    "schemaVersion": 2,
                    "manifests": [image, wrong_reference_attestation],
                },
                sort_keys=True,
            ).encode()
            wrong_reference_artifact = (
                Path(temporary) / "wrong-attestation-reference.oci.tar"
            )
            self._write_oci(wrong_reference_artifact, wrong_reference_index, blobs)
            with self.assertRaisesRegex(
                ValueError, "reference digest does not match candidate image"
            ):
                OCI_VERIFY.verify(
                    wrong_reference_artifact, "1.3.0", "a" * 40
                )

            invalid_subjects = (
                (
                    "missing-media",
                    {
                        key: value
                        for key, value in image_subject.items()
                        if key != "mediaType"
                    },
                ),
                ("wrong-size", {**image_subject, "size": image_subject["size"] + 1}),
                ("wrong-media", {**image_subject, "mediaType": "text/plain"}),
            )
            for label, invalid_subject in invalid_subjects:
                invalid_subject_attestation = self._json_blob(
                    blobs,
                    {
                        "schemaVersion": 2,
                        "config": attestation_config,
                        "subject": invalid_subject,
                        "layers": statements,
                    },
                )
                invalid_subject_attestation["platform"] = {
                    "os": "unknown",
                    "architecture": "unknown",
                }
                invalid_subject_attestation["annotations"] = dict(
                    attestation["annotations"]
                )
                invalid_subject_index = json.dumps(
                    {
                        "schemaVersion": 2,
                        "manifests": [image, invalid_subject_attestation],
                    },
                    sort_keys=True,
                ).encode()
                invalid_subject_artifact = (
                    Path(temporary) / f"invalid-subject-{label}.oci.tar"
                )
                self._write_oci(
                    invalid_subject_artifact, invalid_subject_index, blobs
                )
                with self.assertRaisesRegex(
                    ValueError, "subject descriptor does not match candidate image"
                ):
                    OCI_VERIFY.verify(
                        invalid_subject_artifact, "1.3.0", "a" * 40
                    )

            for layout_version, expected_error in (
                (None, "missing oci-layout"),
                ("2.0.0", "unsupported OCI image layout version"),
            ):
                invalid_layout = Path(temporary) / f"layout-{layout_version}.oci.tar"
                self._write_oci(invalid_layout, index, blobs, layout_version)
                with self.assertRaisesRegex(ValueError, expected_error):
                    OCI_VERIFY.verify(invalid_layout, "1.3.0", "a" * 40)

            missing_index_schema = Path(temporary) / "missing-index-schema.oci.tar"
            self._write_oci(
                missing_index_schema,
                json.dumps({"manifests": [image, attestation]}).encode(),
                blobs,
            )
            with self.assertRaisesRegex(
                ValueError, "OCI index schemaVersion must be 2"
            ):
                OCI_VERIFY.verify(missing_index_schema, "1.3.0", "a" * 40)

            no_schema_image = self._json_blob(
                blobs,
                {"config": config, "layers": [layer]},
            )
            no_schema_image["platform"] = {
                "os": "linux",
                "architecture": "amd64",
            }
            missing_manifest_schema = (
                Path(temporary) / "missing-manifest-schema.oci.tar"
            )
            self._write_oci(
                missing_manifest_schema,
                json.dumps(
                    {
                        "schemaVersion": 2,
                        "manifests": [no_schema_image, attestation],
                    }
                ).encode(),
                blobs,
            )
            with self.assertRaisesRegex(
                ValueError, "OCI manifest schemaVersion must be 2"
            ):
                OCI_VERIFY.verify(missing_manifest_schema, "1.3.0", "a" * 40)

            unsupported_image = dict(image)
            unsupported_image["mediaType"] = "text/plain"
            unsupported_media_index = json.dumps(
                {
                    "schemaVersion": 2,
                    "manifests": [unsupported_image, attestation],
                },
                sort_keys=True,
            ).encode()
            unsupported_media = Path(temporary) / "unsupported-media.oci.tar"
            self._write_oci(unsupported_media, unsupported_media_index, blobs)
            with self.assertRaisesRegex(
                ValueError, "unsupported OCI manifest descriptor media type"
            ):
                OCI_VERIFY.verify(unsupported_media, "1.3.0", "a" * 40)

            artifact_typed_image = dict(image)
            artifact_typed_image["mediaType"] = (
                "application/vnd.oci.artifact.manifest.v1+json"
            )
            artifact_typed_index = json.dumps(
                {
                    "schemaVersion": 2,
                    "manifests": [artifact_typed_image, attestation],
                },
                sort_keys=True,
            ).encode()
            artifact_typed_candidate = (
                Path(temporary) / "artifact-typed-image.oci.tar"
            )
            self._write_oci(
                artifact_typed_candidate, artifact_typed_index, blobs
            )
            with self.assertRaisesRegex(
                ValueError, "unrecognized OCI index manifest structure"
            ):
                OCI_VERIFY.verify(
                    artifact_typed_candidate, "1.3.0", "a" * 40
                )

            with self.assertRaisesRegex(ValueError, "not a valid tar archive"):
                OCI_VERIFY._layer_diff_id(
                    b"not-a-tar", "application/vnd.oci.image.layer.v1.tar"
                )

            bad_rootfs_config = self._json_blob(
                blobs,
                {
                    "architecture": "amd64",
                    "os": "linux",
                    "rootfs": {"type": "layers", "diff_ids": ["sha256:" + "b" * 64]},
                    "config": json.loads(
                        blobs[OCI_VERIFY._blob_path(config["digest"])]
                    )["config"],
                },
                "application/vnd.oci.image.config.v1+json",
            )
            bad_rootfs_image = self._json_blob(
                blobs,
                {"schemaVersion": 2, "config": bad_rootfs_config, "layers": [layer]},
            )
            bad_rootfs_image["platform"] = image["platform"]
            bad_rootfs_index = json.dumps(
                {
                    "schemaVersion": 2,
                    "manifests": [bad_rootfs_image, attestation],
                },
                sort_keys=True,
            ).encode()
            bad_rootfs = Path(temporary) / "bad-rootfs.oci.tar"
            self._write_oci(bad_rootfs, bad_rootfs_index, blobs)
            with self.assertRaisesRegex(ValueError, "rootfs diff_ids"):
                OCI_VERIFY.verify(bad_rootfs, "1.3.0", "a" * 40)

            image_without_platform = dict(image)
            image_without_platform.pop("platform")
            index_without_platform = json.dumps(
                {"schemaVersion": 2,
                 "manifests": [image_without_platform, attestation]},
                sort_keys=True,
            ).encode()
            no_platform = Path(temporary) / "candidate-no-platform.oci.tar"
            self._write_oci(no_platform, index_without_platform, blobs)
            OCI_VERIFY.verify(no_platform, "1.3.0", "a" * 40)

            unmarked_attestation_manifest = self._json_blob(
                blobs,
                {"schemaVersion": 2, "subject": dict(image_subject),
                 "blobs": statements},
                "application/vnd.oci.artifact.manifest.v1+json",
            )
            unmarked_index = json.dumps(
                {"schemaVersion": 2,
                 "manifests": [
                     image_without_platform, unmarked_attestation_manifest
                 ]},
                sort_keys=True,
            ).encode()
            unmarked = Path(temporary) / "candidate-unmarked-attestation.oci.tar"
            self._write_oci(unmarked, unmarked_index, blobs)
            OCI_VERIFY.verify(unmarked, "1.3.0", "a" * 40)

            for label, invalid_manifest, expected_error in (
                (
                    "missing",
                    {
                        "schemaVersion": 2,
                        "subject": dict(image_subject),
                        "layers": statements,
                    },
                    "has no config descriptor",
                ),
                (
                    "malformed",
                    {
                        "schemaVersion": 2,
                        "config": None,
                        "subject": dict(image_subject),
                        "layers": statements,
                    },
                    "config descriptor must be an object",
                ),
                (
                    "wrong-media",
                    {
                        "schemaVersion": 2,
                        "config": {
                            **attestation_config,
                            "mediaType": "text/plain",
                        },
                        "subject": dict(image_subject),
                        "layers": statements,
                    },
                    "config has an invalid media type",
                ),
            ):
                invalid_attestation = self._json_blob(blobs, invalid_manifest)
                invalid_attestation["platform"] = {
                    "os": "unknown",
                    "architecture": "unknown",
                }
                invalid_attestation["annotations"] = {
                    "vnd.docker.reference.type": "attestation-manifest",
                    "vnd.docker.reference.digest": image_digest,
                }
                invalid_index = json.dumps(
                    {
                        "schemaVersion": 2,
                        "manifests": [image_without_platform, invalid_attestation],
                    },
                    sort_keys=True,
                ).encode()
                invalid = Path(temporary) / f"candidate-{label}-config.oci.tar"
                self._write_oci(invalid, invalid_index, blobs)
                with self.assertRaisesRegex(ValueError, expected_error):
                    OCI_VERIFY.verify(invalid, "1.3.0", "a" * 40)

            non_json_config_payload = b"not-json"
            non_json_config_digest = hashlib.sha256(
                non_json_config_payload
            ).hexdigest()
            non_json_config = {
                "mediaType": "application/vnd.oci.image.config.v1+json",
                "digest": f"sha256:{non_json_config_digest}",
                "size": len(non_json_config_payload),
            }
            blobs[f"blobs/sha256/{non_json_config_digest}"] = (
                non_json_config_payload
            )
            non_json_attestation = self._json_blob(
                blobs,
                {
                    "schemaVersion": 2,
                    "config": non_json_config,
                    "subject": dict(image_subject),
                    "layers": statements,
                },
            )
            non_json_attestation["platform"] = {
                "os": "unknown",
                "architecture": "unknown",
            }
            non_json_attestation["annotations"] = dict(
                attestation["annotations"]
            )
            non_json_index = json.dumps(
                {
                    "schemaVersion": 2,
                    "manifests": [image_without_platform, non_json_attestation],
                },
                sort_keys=True,
            ).encode()
            non_json_artifact = Path(temporary) / "non-json-config.oci.tar"
            self._write_oci(non_json_artifact, non_json_index, blobs)
            with self.assertRaises(ValueError):
                OCI_VERIFY.verify(non_json_artifact, "1.3.0", "a" * 40)

            for label, revision_value in (
                ("missing", None),
                ("mismatched", "b" * 40),
            ):
                provenance_statement = json.loads(
                    blobs[OCI_VERIFY._blob_path(statements[1]["digest"])]
                )
                external_parameters = provenance_statement["predicate"][
                    "buildDefinition"
                ]["externalParameters"]
                if revision_value is None:
                    external_parameters.clear()
                    provenance_statement["predicate"]["runDetails"][
                        "metadata"
                    ] = {"OPENQP_REVISION": "a" * 40}
                else:
                    external_parameters["args"][
                        "build-arg:OPENQP_REVISION"
                    ] = revision_value
                invalid_provenance = self._json_blob(
                    blobs,
                    provenance_statement,
                    "application/vnd.in-toto+json",
                )
                invalid_provenance_attestation = self._json_blob(
                    blobs,
                    {
                        "schemaVersion": 2,
                        "config": attestation_config,
                        "subject": dict(image_subject),
                        "layers": [statements[0], invalid_provenance],
                    },
                )
                invalid_provenance_attestation["platform"] = {
                    "os": "unknown",
                    "architecture": "unknown",
                }
                invalid_provenance_attestation["annotations"] = dict(
                    attestation["annotations"]
                )
                invalid_provenance_index = json.dumps(
                    {
                        "schemaVersion": 2,
                        "manifests": [
                            image_without_platform,
                            invalid_provenance_attestation,
                        ],
                    },
                    sort_keys=True,
                ).encode()
                invalid_provenance_artifact = (
                    Path(temporary) / f"{label}-provenance-revision.oci.tar"
                )
                self._write_oci(
                    invalid_provenance_artifact,
                    invalid_provenance_index,
                    blobs,
                )
                with self.assertRaisesRegex(ValueError, "provenance.*revision"):
                    OCI_VERIFY.verify(
                        invalid_provenance_artifact, "1.3.0", "a" * 40
                    )

            nested_index_descriptor = self._json_blob(
                blobs,
                {"schemaVersion": 2,
                 "manifests": [image_without_platform, attestation]},
                "application/vnd.oci.image.index.v1+json",
            )
            nested_index = json.dumps(
                {"schemaVersion": 2, "manifests": [nested_index_descriptor]},
                sort_keys=True,
            ).encode()
            nested = Path(temporary) / "candidate-nested-index.oci.tar"
            self._write_oci(nested, nested_index, blobs)
            OCI_VERIFY.verify(nested, "1.3.0", "a" * 40)

            image_with_wrong_platform = dict(image)
            image_with_wrong_platform["platform"] = {
                "os": "linux", "architecture": "arm64"
            }
            wrong_platform_index = json.dumps(
                {"schemaVersion": 2,
                 "manifests": [image_with_wrong_platform, attestation]},
                sort_keys=True,
            ).encode()
            wrong_platform = Path(temporary) / "wrong-platform.oci.tar"
            self._write_oci(wrong_platform, wrong_platform_index, blobs)
            with self.assertRaisesRegex(ValueError, "platform mismatch"):
                OCI_VERIFY.verify(wrong_platform, "1.3.0", "a" * 40)

            for empty_platform in ({}, None):
                image_with_empty_platform = dict(image)
                image_with_empty_platform["platform"] = empty_platform
                empty_platform_index = json.dumps(
                    {"schemaVersion": 2,
                     "manifests": [image_with_empty_platform, attestation]},
                    sort_keys=True,
                ).encode()
                empty_platform_artifact = (
                    Path(temporary) / f"empty-platform-{empty_platform is None}.oci.tar"
                )
                self._write_oci(
                    empty_platform_artifact, empty_platform_index, blobs
                )
                with self.assertRaisesRegex(ValueError, "platform mismatch"):
                    OCI_VERIFY.verify(
                        empty_platform_artifact, "1.3.0", "a" * 40
                    )

            tampered_blobs = dict(blobs)
            config_path = OCI_VERIFY._blob_path(config["digest"])
            tampered_blobs[config_path] = tampered_blobs[config_path].replace(
                b'"amd64"', b'"arm64"'
            )
            tampered = Path(temporary) / "tampered.oci.tar"
            self._write_oci(tampered, index, tampered_blobs)
            with self.assertRaisesRegex(ValueError, "descriptor digest mismatch"):
                OCI_VERIFY.verify(tampered, "1.3.0", "a" * 40)

            tampered_layer_blobs = dict(blobs)
            tampered_layer_payload = bytearray(layer_payload)
            tampered_layer_payload[0] ^= 1
            tampered_layer_blobs[OCI_VERIFY._blob_path(layer["digest"])] = bytes(
                tampered_layer_payload
            )
            self.assertEqual(
                len(tampered_layer_blobs[OCI_VERIFY._blob_path(layer["digest"])]),
                layer["size"],
            )
            tampered_layer = Path(temporary) / "tampered-layer.oci.tar"
            self._write_oci(tampered_layer, index, tampered_layer_blobs)
            with self.assertRaisesRegex(ValueError, "descriptor digest mismatch"):
                OCI_VERIFY.verify(tampered_layer, "1.3.0", "a" * 40)

            forged_spdx = self._json_blob(
                blobs,
                {
                    "_type": "https://in-toto.io/Statement/v1",
                    "subject": [{
                        "name": "openqp/openqp",
                        "digest": {"sha256": image_digest.split(":", 1)[1]},
                    }],
                    "predicateType": "https://example.test/not-spdx",
                    "predicate": {},
                },
                "application/vnd.in-toto+json",
            )
            forged_provenance = self._json_blob(
                blobs,
                {
                    "_type": "https://in-toto.io/Statement/v1",
                    "subject": [{
                        "name": "openqp/openqp",
                        "digest": {"sha256": image_digest.split(":", 1)[1]},
                    }],
                    "predicateType": "https://slsa.dev/provenance-forged",
                    "predicate": {},
                },
                "application/vnd.in-toto+json",
            )
            for name, layers, error in (
                ("forged-spdx", [forged_spdx, statements[1]], "no SPDX SBOM"),
                ("forged-slsa", [statements[0], forged_provenance], "no SLSA provenance"),
            ):
                forged_attestation = self._json_blob(
                    blobs,
                    {
                        "schemaVersion": 2,
                        "config": attestation_config,
                        "subject": dict(image_subject),
                        "layers": layers,
                    },
                )
                forged_attestation["platform"] = {
                    "os": "unknown", "architecture": "unknown"
                }
                forged_attestation["annotations"] = dict(attestation["annotations"])
                forged_index = json.dumps(
                    {"schemaVersion": 2,
                     "manifests": [image, forged_attestation]},
                    sort_keys=True,
                ).encode()
                forged_artifact = Path(temporary) / f"{name}.oci.tar"
                self._write_oci(forged_artifact, forged_index, blobs)
                with self.assertRaisesRegex(ValueError, error):
                    OCI_VERIFY.verify(
                        forged_artifact, "1.3.0", "a" * 40
                    )

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
                    "predicate": {
                        "spdxVersion": "SPDX-2.3",
                        "SPDXID": "SPDXRef-DOCUMENT",
                        "dataLicense": "CC0-1.0",
                        "name": "Unrelated image SBOM",
                        "documentNamespace": "https://openqp.dev/spdx/unrelated",
                        "creationInfo": {
                            "created": "2026-08-11T00:00:00Z",
                            "creators": ["Tool: BuildKit"],
                        },
                    },
                },
                "application/vnd.in-toto+json",
            )
            unrelated_attestation = self._json_blob(
                blobs,
                {
                    "schemaVersion": 2,
                    "config": attestation_config,
                    "subject": dict(image_subject),
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

            malformed_statements = (
                (
                    "wrong-type",
                    {
                        **json.loads(
                            blobs[OCI_VERIFY._blob_path(statements[0]["digest"])]
                        ),
                        "_type": "https://example.test/not-in-toto",
                    },
                    "unsupported in-toto statement type",
                ),
                (
                    "empty-spdx",
                    {
                        **json.loads(
                            blobs[OCI_VERIFY._blob_path(statements[0]["digest"])]
                        ),
                        "predicate": {},
                    },
                    "unsupported SPDX document version",
                ),
                (
                    "incomplete-spdx",
                    {
                        **json.loads(
                            blobs[OCI_VERIFY._blob_path(statements[0]["digest"])]
                        ),
                        "predicate": {
                            "spdxVersion": "SPDX-2.3",
                            "SPDXID": "SPDXRef-DOCUMENT",
                        },
                    },
                    "SPDX predicate lacks required document field",
                ),
                (
                    "empty-slsa",
                    {
                        **json.loads(
                            blobs[OCI_VERIFY._blob_path(statements[1]["digest"])]
                        ),
                        "predicate": {},
                    },
                    "SLSA v1 predicate lacks",
                ),
                (
                    "missing-slsa-parameters",
                    {
                        **json.loads(
                            blobs[OCI_VERIFY._blob_path(statements[1]["digest"])]
                        ),
                        "predicate": {
                            "buildDefinition": {
                                "buildType": "https://mobyproject.org/buildkit@v1"
                            },
                            "runDetails": {
                                "builder": {
                                    "id": "https://mobyproject.org/buildkit"
                                }
                            },
                        },
                    },
                    "buildDefinition.externalParameters",
                ),
                (
                    "bad-spdx-created",
                    {
                        **json.loads(
                            blobs[OCI_VERIFY._blob_path(statements[0]["digest"])]
                        ),
                        "predicate": {
                            **json.loads(
                                blobs[
                                    OCI_VERIFY._blob_path(statements[0]["digest"])
                                ]
                            )["predicate"],
                            "creationInfo": {
                                "created": "x",
                                "creators": ["Tool: BuildKit"],
                            },
                        },
                    },
                    "lacks valid creationInfo",
                ),
                (
                    "bad-spdx-creator",
                    {
                        **json.loads(
                            blobs[OCI_VERIFY._blob_path(statements[0]["digest"])]
                        ),
                        "predicate": {
                            **json.loads(
                                blobs[
                                    OCI_VERIFY._blob_path(statements[0]["digest"])
                                ]
                            )["predicate"],
                            "creationInfo": {
                                "created": "2026-08-11T00:00:00Z",
                                "creators": ["BuildKit"],
                            },
                        },
                    },
                    "lacks valid creationInfo",
                ),
            )
            for name, malformed_statement, expected_error in malformed_statements:
                malformed_layer = self._json_blob(
                    blobs, malformed_statement, "application/vnd.in-toto+json"
                )
                malformed_attestation = self._json_blob(
                    blobs,
                    {
                        "schemaVersion": 2,
                        "config": attestation_config,
                        "subject": dict(image_subject),
                        "layers": [malformed_layer, statements[1]],
                    },
                )
                malformed_attestation["platform"] = {
                    "os": "unknown", "architecture": "unknown"
                }
                malformed_attestation["annotations"] = dict(attestation["annotations"])
                malformed_index = json.dumps(
                    {"schemaVersion": 2, "manifests": [image, malformed_attestation]},
                    sort_keys=True,
                ).encode()
                malformed_artifact = Path(temporary) / f"{name}.oci.tar"
                self._write_oci(malformed_artifact, malformed_index, blobs)
                with self.assertRaisesRegex(ValueError, expected_error):
                    OCI_VERIFY.verify(
                        malformed_artifact, "1.3.0", "a" * 40
                    )
        self.assertEqual(summary["image_digest"], image_digest)
        self.assertEqual(summary["platform"], "linux/amd64")


if __name__ == "__main__":
    unittest.main()
