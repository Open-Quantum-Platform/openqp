#!/usr/bin/env python3
"""Verify a build-only OpenQP OCI archive and its embedded attestations."""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from pathlib import Path
from typing import Any


SPDX_PREDICATE_TYPES = {"https://spdx.dev/Document"}
SLSA_PROVENANCE_PREDICATE_TYPES = {
    "https://slsa.dev/provenance/v0.2",
    "https://slsa.dev/provenance/v1",
}
OCI_IMAGE_CONFIG_MEDIA_TYPE = "application/vnd.oci.image.config.v1+json"
OCI_INDEX_MEDIA_TYPES = {
    "application/vnd.oci.image.index.v1+json",
    "application/vnd.docker.distribution.manifest.list.v2+json",
}


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


def _leaf_manifests(
    archive: tarfile.TarFile,
    descriptors: Any,
    seen_indexes: set[str] | None = None,
) -> list[tuple[dict[str, Any], dict[str, Any]]]:
    """Authenticate and flatten nested OCI indexes into leaf manifests."""
    if not isinstance(descriptors, list):
        raise ValueError("OCI index manifests must be a list")
    if seen_indexes is None:
        seen_indexes = set()
    leaves = []
    for descriptor in descriptors:
        if not isinstance(descriptor, dict):
            raise ValueError("OCI index manifest descriptor must be an object")
        manifest = _descriptor_json(archive, descriptor)
        media_type = descriptor.get("mediaType")
        if media_type in OCI_INDEX_MEDIA_TYPES:
            digest = str(descriptor.get("digest"))
            if digest in seen_indexes:
                raise ValueError(f"cyclic OCI index descriptor: {digest}")
            nested = manifest.get("manifests")
            if not isinstance(nested, list):
                raise ValueError("OCI index descriptor does not contain manifests")
            seen_indexes.add(digest)
            leaves.extend(_leaf_manifests(archive, nested, seen_indexes))
            seen_indexes.remove(digest)
        else:
            if "manifests" in manifest:
                raise ValueError(
                    "OCI manifest contains a nested index under a non-index media type"
                )
            leaves.append((descriptor, manifest))
    return leaves


def verify(
    artifact: Path, expected_version: str, expected_revision: str
) -> dict[str, Any]:
    with tarfile.open(artifact, mode="r:*") as archive:
        index = _read_json(archive, "index.json")
        descriptors = index.get("manifests", [])
        images = []
        attestations = []
        for descriptor, manifest in _leaf_manifests(archive, descriptors):
            annotations = descriptor.get("annotations", {})
            platform = descriptor.get("platform", {})
            attestation_hint = (
                annotations.get("vnd.docker.reference.type")
                == "attestation-manifest"
                or platform == {"architecture": "unknown", "os": "unknown"}
            )
            config_descriptor = manifest.get("config")
            is_image = (
                isinstance(config_descriptor, dict)
                and config_descriptor.get("mediaType")
                == OCI_IMAGE_CONFIG_MEDIA_TYPE
                and isinstance(manifest.get("layers"), list)
            )
            is_attestation = (
                attestation_hint
                or "subject" in manifest
                or "artifactType" in manifest
                or isinstance(manifest.get("blobs"), list)
            )
            if is_image and is_attestation:
                raise ValueError("OCI manifest is ambiguously image and attestation")
            if is_image:
                images.append((descriptor, manifest))
            elif is_attestation:
                attestations.append((descriptor, manifest))
            else:
                raise ValueError("unrecognized OCI index manifest structure")

        if len(images) != 1:
            raise ValueError(f"expected one OCI image manifest, found {len(images)}")
        image_descriptor, image_manifest = images[0]
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
        for descriptor, manifest in attestations:
            subject = manifest.get("subject")
            if isinstance(subject, dict) and subject.get("digest"):
                subjects.add(str(subject["digest"]))
            annotation_subject = descriptor.get("annotations", {}).get(
                "vnd.docker.reference.digest"
            )
            if annotation_subject:
                subjects.add(str(annotation_subject))
            attestation_config = manifest.get("config")
            if isinstance(attestation_config, dict):
                _descriptor_bytes(archive, attestation_config)
            payload_descriptors = manifest.get("layers")
            if payload_descriptors is None:
                payload_descriptors = manifest.get("blobs")
            if not isinstance(payload_descriptors, list):
                raise ValueError("attestation manifest has no layer/blob descriptors")
            for layer in payload_descriptors:
                if not isinstance(layer, dict):
                    raise ValueError("attestation payload descriptor must be an object")
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
        if not image_predicate_types.intersection(SPDX_PREDICATE_TYPES):
            raise ValueError(
                f"candidate image has no SPDX SBOM; "
                f"image_predicates={image_predicate_types}, all_predicates={predicate_types}"
            )
        if not image_predicate_types.intersection(
            SLSA_PROVENANCE_PREDICATE_TYPES
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
