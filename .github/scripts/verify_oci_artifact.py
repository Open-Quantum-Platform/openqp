#!/usr/bin/env python3
"""Verify a build-only OpenQP OCI archive and its embedded attestations."""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from pathlib import Path
from typing import Any


def _blob_path(digest: str) -> str:
    algorithm, value = digest.split(":", 1)
    if algorithm != "sha256" or len(value) != 64:
        raise ValueError(f"unsupported OCI digest: {digest}")
    return f"blobs/{algorithm}/{value}"


def _read_json(archive: tarfile.TarFile, member_name: str) -> dict[str, Any]:
    member = archive.getmember(member_name)
    handle = archive.extractfile(member)
    if handle is None:
        raise ValueError(f"OCI member is not a regular file: {member_name}")
    return json.load(handle)


def _descriptor_bytes(
    archive: tarfile.TarFile, descriptor: dict[str, Any]
) -> bytes:
    digest = str(descriptor["digest"])
    member_name = _blob_path(digest)
    member = archive.getmember(member_name)
    handle = archive.extractfile(member)
    if handle is None or not member.isfile():
        raise ValueError(f"OCI blob is not a regular file: {member_name}")
    payload = handle.read()

    declared_size = descriptor.get("size")
    if not isinstance(declared_size, int) or declared_size < 0:
        raise ValueError(f"invalid OCI descriptor size for {digest}: {declared_size}")
    if len(payload) != declared_size:
        raise ValueError(
            f"OCI descriptor size mismatch for {digest}: "
            f"declared {declared_size}, actual {len(payload)}"
        )
    actual_digest = "sha256:" + hashlib.sha256(payload).hexdigest()
    if actual_digest != digest:
        raise ValueError(
            f"OCI descriptor digest mismatch: declared {digest}, actual {actual_digest}"
        )
    return payload


def _descriptor_json(
    archive: tarfile.TarFile, descriptor: dict[str, Any]
) -> dict[str, Any]:
    digest = str(descriptor["digest"])
    value = json.loads(_descriptor_bytes(archive, descriptor))
    if not isinstance(value, dict):
        raise ValueError(f"OCI descriptor does not contain a JSON object: {digest}")
    return value


def _archive_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify(
    artifact: Path, expected_version: str, expected_revision: str
) -> dict[str, Any]:
    with tarfile.open(artifact, mode="r:*") as archive:
        index = _read_json(archive, "index.json")
        descriptors = index.get("manifests", [])
        images = []
        attestations = []
        for descriptor in descriptors:
            annotations = descriptor.get("annotations", {})
            platform = descriptor.get("platform", {})
            if (
                annotations.get("vnd.docker.reference.type")
                == "attestation-manifest"
                or platform == {"architecture": "unknown", "os": "unknown"}
            ):
                attestations.append(descriptor)
            else:
                images.append(descriptor)

        if len(images) != 1:
            raise ValueError(f"expected one OCI image manifest, found {len(images)}")
        image_descriptor = images[0]
        image_manifest = _descriptor_json(archive, image_descriptor)
        config = _descriptor_json(archive, image_manifest["config"])
        config_platform = {
            "os": config.get("os"),
            "architecture": config.get("architecture"),
        }
        expected_platform = {"os": "linux", "architecture": "amd64"}
        if config_platform != expected_platform:
            raise ValueError(f"unexpected OCI image platform: {config_platform}")
        descriptor_platform = image_descriptor.get("platform")
        # OCI descriptors may omit the optional platform object for a
        # single-platform archive.  When present it must agree with the image
        # configuration, which is the authoritative runtime platform record.
        if (
            "platform" in image_descriptor
            and descriptor_platform != config_platform
        ):
            raise ValueError(
                "OCI descriptor/config platform mismatch: "
                f"{descriptor_platform} != {config_platform}"
            )
        layers = image_manifest.get("layers")
        if not isinstance(layers, list):
            raise ValueError("OCI image manifest layers must be a list")
        for layer in layers:
            if not isinstance(layer, dict):
                raise ValueError("OCI image layer descriptor must be an object")
            _descriptor_bytes(archive, layer)
        labels = config.get("config", {}).get("Labels", {})
        expected_labels = {
            "org.opencontainers.image.source": (
                "https://github.com/Open-Quantum-Platform/openqp"
            ),
            "org.opencontainers.image.version": expected_version,
            "org.opencontainers.image.revision": expected_revision,
            "org.opencontainers.image.licenses": (
                "LicenseRef-OpenQP-Research-1.0"
            ),
        }
        mismatches = {
            name: {"actual": labels.get(name), "expected": expected}
            for name, expected in expected_labels.items()
            if labels.get(name) != expected
        }
        if mismatches:
            raise ValueError(f"OCI labels do not match artifact inputs: {mismatches}")

        image_digest = str(image_descriptor["digest"])
        predicate_types: set[str] = set()
        image_predicate_types: set[str] = set()
        subjects: set[str] = set()
        if not attestations:
            raise ValueError("OCI archive contains no attestation manifest")
        for descriptor in attestations:
            manifest = _descriptor_json(archive, descriptor)
            subject = manifest.get("subject")
            if isinstance(subject, dict) and subject.get("digest"):
                subjects.add(str(subject["digest"]))
            annotation_subject = descriptor.get("annotations", {}).get(
                "vnd.docker.reference.digest"
            )
            if annotation_subject:
                subjects.add(str(annotation_subject))
            for layer in manifest.get("layers", []):
                statement = _descriptor_json(archive, layer)
                predicate_type = statement.get("predicateType")
                if predicate_type:
                    predicate_types.add(str(predicate_type))
                statement_subjects: set[str] = set()
                for statement_subject in statement.get("subject", []):
                    digest = statement_subject.get("digest", {}).get("sha256")
                    if digest:
                        subject_digest = f"sha256:{digest}"
                        statement_subjects.add(subject_digest)
                        subjects.add(subject_digest)
                if predicate_type and image_digest in statement_subjects:
                    image_predicate_types.add(str(predicate_type))

        if image_digest not in subjects:
            raise ValueError(
                f"attestations do not identify image {image_digest}; subjects={subjects}"
            )
        if not any("spdx" in value.lower() for value in image_predicate_types):
            raise ValueError(
                f"candidate image has no SPDX SBOM; "
                f"image_predicates={image_predicate_types}, all_predicates={predicate_types}"
            )
        if not any(
            "slsa.dev/provenance" in value for value in image_predicate_types
        ):
            raise ValueError(
                f"candidate image has no SLSA provenance; "
                f"image_predicates={image_predicate_types}, all_predicates={predicate_types}"
            )

    return {
        "archive": artifact.name,
        "archive_sha256": _archive_sha256(artifact),
        "image_digest": image_digest,
        "platform": "linux/amd64",
        "version": expected_version,
        "revision": expected_revision,
        "predicate_types": sorted(image_predicate_types),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact", type=Path)
    parser.add_argument("--version", required=True)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    summary = verify(args.artifact, args.version, args.revision)
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
