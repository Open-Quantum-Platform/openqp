#!/usr/bin/env python3
"""Verify a build-only OpenQP OCI archive and its embedded attestations."""

from __future__ import annotations

import argparse
import datetime
import gzip
import hashlib
import io
import json
import re
import tarfile
from pathlib import Path
from typing import Any


SPDX_PREDICATE_TYPES = {"https://spdx.dev/Document"}
SUPPORTED_SPDX_VERSIONS = {"SPDX-2.2", "SPDX-2.3"}
SLSA_PROVENANCE_PREDICATE_TYPES = {
    "https://slsa.dev/provenance/v0.2",
    "https://slsa.dev/provenance/v1",
}
IN_TOTO_STATEMENT_TYPES = {
    "https://in-toto.io/Statement/v0.1",
    "https://in-toto.io/Statement/v1",
}
OCI_IMAGE_LAYOUT_VERSION = "1.0.0"
OCI_IMAGE_CONFIG_MEDIA_TYPE = "application/vnd.oci.image.config.v1+json"
OCI_INDEX_MEDIA_TYPES = {
    "application/vnd.oci.image.index.v1+json",
    "application/vnd.docker.distribution.manifest.list.v2+json",
}
OCI_LEAF_MANIFEST_MEDIA_TYPES = {
    "application/vnd.oci.image.manifest.v1+json",
    "application/vnd.oci.artifact.manifest.v1+json",
    "application/vnd.docker.distribution.manifest.v2+json",
}
UNCOMPRESSED_LAYER_MEDIA_TYPES = {
    "application/vnd.oci.image.layer.v1.tar",
    "application/vnd.oci.image.layer.nondistributable.v1.tar",
}
GZIP_LAYER_MEDIA_TYPES = {
    "application/vnd.oci.image.layer.v1.tar+gzip",
    "application/vnd.oci.image.layer.nondistributable.v1.tar+gzip",
    "application/vnd.docker.image.rootfs.diff.tar.gzip",
    "application/vnd.docker.image.rootfs.foreign.diff.tar.gzip",
}


def _blob_path(digest: str) -> str:
    algorithm, value = digest.split(":", 1)
    if algorithm != "sha256" or len(value) != 64:
        raise ValueError(f"unsupported OCI digest: {digest}")
    return f"blobs/{algorithm}/{value}"


def _read_json(archive: tarfile.TarFile, member_name: str) -> dict[str, Any]:
    try:
        member = archive.getmember(member_name)
    except KeyError as exc:
        raise ValueError(f"OCI archive is missing {member_name}") from exc
    handle = archive.extractfile(member)
    if handle is None or not member.isfile():
        raise ValueError(f"OCI member is not a regular file: {member_name}")
    value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"OCI member does not contain a JSON object: {member_name}")
    return value


def _statement_subjects(statement: dict[str, Any]) -> tuple[str, set[str]]:
    """Validate an in-toto statement and return its typed SHA-256 subjects."""
    statement_type = statement.get("_type")
    if statement_type not in IN_TOTO_STATEMENT_TYPES:
        raise ValueError(f"unsupported in-toto statement type: {statement_type}")

    predicate_type = statement.get("predicateType")
    if not isinstance(predicate_type, str) or not predicate_type:
        raise ValueError("in-toto statement has no predicateType URI")
    predicate = statement.get("predicate")
    if not isinstance(predicate, dict):
        raise ValueError(f"{predicate_type} predicate must be an object")

    if predicate_type in SPDX_PREDICATE_TYPES:
        spdx_version = predicate.get("spdxVersion")
        if spdx_version not in SUPPORTED_SPDX_VERSIONS:
            raise ValueError(f"unsupported SPDX document version: {spdx_version}")
        if predicate.get("SPDXID") != "SPDXRef-DOCUMENT":
            raise ValueError("SPDX predicate is not an SPDX document")
        for field in ("dataLicense", "name", "documentNamespace"):
            if (
                not isinstance(predicate.get(field), str)
                or not predicate[field].strip()
            ):
                raise ValueError(f"SPDX predicate lacks required document field: {field}")
        if predicate["dataLicense"] != "CC0-1.0":
            raise ValueError("SPDX document dataLicense must be CC0-1.0")
        creation_info = predicate.get("creationInfo")
        creators = creation_info.get("creators") if isinstance(creation_info, dict) else None
        created = creation_info.get("created") if isinstance(creation_info, dict) else None
        valid_creators = (
            isinstance(creators, list)
            and bool(creators)
            and all(
                isinstance(creator, str)
                and any(
                    creator.startswith(f"{kind}: ")
                    and bool(creator.removeprefix(f"{kind}: ").strip())
                    for kind in ("Person", "Organization", "Tool")
                )
                for creator in creators
            )
        )
        if (
            not isinstance(creation_info, dict)
            or not isinstance(created, str)
            or not re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", created)
            or not valid_creators
        ):
            raise ValueError("SPDX predicate lacks valid creationInfo")
        try:
            datetime.datetime.strptime(created, "%Y-%m-%dT%H:%M:%SZ")
        except ValueError as exc:
            raise ValueError("SPDX predicate has an invalid creation timestamp") from exc
    elif predicate_type == "https://slsa.dev/provenance/v0.2":
        builder = predicate.get("builder")
        invocation = predicate.get("invocation")
        metadata = predicate.get("metadata")
        materials = predicate.get("materials")
        if (
            not isinstance(builder, dict)
            or not isinstance(builder.get("id"), str)
            or not builder["id"]
            or not isinstance(predicate.get("buildType"), str)
            or not predicate["buildType"]
            or not isinstance(invocation, dict)
            or not isinstance(metadata, dict)
            or not isinstance(materials, list)
            or not all(isinstance(material, dict) for material in materials)
        ):
            raise ValueError(
                "SLSA v0.2 predicate lacks builder.id, buildType, invocation, "
                "metadata, or materials"
            )
    elif predicate_type == "https://slsa.dev/provenance/v1":
        build_definition = predicate.get("buildDefinition")
        run_details = predicate.get("runDetails")
        builder = run_details.get("builder") if isinstance(run_details, dict) else None
        external_parameters = (
            build_definition.get("externalParameters")
            if isinstance(build_definition, dict)
            else None
        )
        if (
            not isinstance(build_definition, dict)
            or not isinstance(build_definition.get("buildType"), str)
            or not build_definition["buildType"]
            or not isinstance(external_parameters, dict)
            or not isinstance(builder, dict)
            or not isinstance(builder.get("id"), str)
            or not builder["id"]
        ):
            raise ValueError(
                "SLSA v1 predicate lacks buildDefinition.buildType, "
                "buildDefinition.externalParameters, or "
                "runDetails.builder.id"
            )
        if (
            "internalParameters" in build_definition
            and not isinstance(build_definition["internalParameters"], dict)
        ):
            raise ValueError(
                "SLSA v1 buildDefinition.internalParameters must be an object"
            )
        if (
            "resolvedDependencies" in build_definition
            and not isinstance(build_definition["resolvedDependencies"], list)
        ):
            raise ValueError(
                "SLSA v1 buildDefinition.resolvedDependencies must be a list"
            )
        if "metadata" in run_details and not isinstance(run_details["metadata"], dict):
            raise ValueError("SLSA v1 runDetails.metadata must be an object")

    raw_subjects = statement.get("subject")
    if not isinstance(raw_subjects, list) or not raw_subjects:
        raise ValueError("in-toto statement subject must be a non-empty list")
    subjects: set[str] = set()
    for subject in raw_subjects:
        if not isinstance(subject, dict) or not isinstance(subject.get("digest"), dict):
            raise ValueError("in-toto statement subject lacks a digest object")
        digest = subject["digest"].get("sha256")
        if digest is not None:
            if (
                not isinstance(digest, str)
                or len(digest) != 64
                or any(character not in "0123456789abcdef" for character in digest)
            ):
                raise ValueError(f"invalid in-toto subject SHA-256: {digest}")
            subjects.add(f"sha256:{digest}")
    return predicate_type, subjects


def _descriptor_bytes(
    archive: tarfile.TarFile, descriptor: dict[str, Any]
) -> bytes:
    media_type = descriptor.get("mediaType")
    if not isinstance(media_type, str) or not media_type:
        raise ValueError("OCI descriptor has no mediaType")
    digest = str(descriptor["digest"])
    member_name = _blob_path(digest)
    member = archive.getmember(member_name)
    handle = archive.extractfile(member)
    if handle is None or not member.isfile():
        raise ValueError(f"OCI blob is not a regular file: {member_name}")
    payload = handle.read()

    declared_size = descriptor.get("size")
    if (
        not isinstance(declared_size, int)
        or isinstance(declared_size, bool)
        or declared_size < 0
    ):
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


def _require_schema_version(document: dict[str, Any], label: str) -> None:
    if document.get("schemaVersion") != 2:
        raise ValueError(f"{label} schemaVersion must be 2")


def _layer_diff_id(payload: bytes, media_type: Any) -> str:
    if media_type in UNCOMPRESSED_LAYER_MEDIA_TYPES:
        uncompressed = payload
    elif media_type in GZIP_LAYER_MEDIA_TYPES:
        try:
            uncompressed = gzip.decompress(payload)
        except (EOFError, OSError) as exc:
            raise ValueError(f"invalid gzip OCI image layer: {media_type}") from exc
    else:
        raise ValueError(f"unsupported OCI image layer media type: {media_type}")
    try:
        with tarfile.open(fileobj=io.BytesIO(uncompressed), mode="r:") as layer:
            # Reading the complete member table validates every tar header and
            # catches truncated member payloads without extracting anything.
            layer.getmembers()
    except (OSError, tarfile.TarError) as exc:
        raise ValueError(f"OCI image layer is not a valid tar archive: {media_type}") from exc
    return "sha256:" + hashlib.sha256(uncompressed).hexdigest()


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
        media_type = descriptor.get("mediaType")
        if media_type not in OCI_INDEX_MEDIA_TYPES | OCI_LEAF_MANIFEST_MEDIA_TYPES:
            raise ValueError(
                f"unsupported OCI manifest descriptor media type: {media_type}"
            )
        manifest = _descriptor_json(archive, descriptor)
        _require_schema_version(manifest, "OCI manifest")
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
        layout = _read_json(archive, "oci-layout")
        layout_version = layout.get("imageLayoutVersion")
        if layout_version != OCI_IMAGE_LAYOUT_VERSION:
            raise ValueError(
                f"unsupported OCI image layout version: {layout_version}"
            )
        index = _read_json(archive, "index.json")
        _require_schema_version(index, "OCI index")
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
            # BuildKit exports attestation manifests in OCI-image shape: they
            # may carry both an image-config descriptor and a layers list.
            # Explicit attestation markers must therefore take precedence;
            # the attestation config and payload descriptors are still
            # authenticated below before their contents are trusted.
            if is_attestation:
                attestations.append((descriptor, manifest))
            elif is_image:
                images.append((descriptor, manifest))
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
        rootfs = config.get("rootfs")
        diff_ids = rootfs.get("diff_ids") if isinstance(rootfs, dict) else None
        if (
            not isinstance(rootfs, dict)
            or rootfs.get("type") != "layers"
            or not isinstance(diff_ids, list)
            or len(diff_ids) != len(layers)
        ):
            raise ValueError(
                "OCI image config rootfs must contain one diff_id per manifest layer"
            )
        for layer in layers:
            if not isinstance(layer, dict):
                raise ValueError("OCI image layer descriptor must be an object")
        actual_diff_ids = [
            _layer_diff_id(
                _descriptor_bytes(archive, layer), layer.get("mediaType")
            )
            for layer in layers
        ]
        if diff_ids != actual_diff_ids:
            raise ValueError(
                "OCI image config rootfs diff_ids do not match manifest layers: "
                f"declared {diff_ids}, actual {actual_diff_ids}"
            )
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
            has_layers = "layers" in manifest
            attestation_config = manifest.get("config")
            if "config" in manifest:
                if not isinstance(attestation_config, dict):
                    raise ValueError(
                        "attestation config descriptor must be an object"
                    )
                _descriptor_bytes(archive, attestation_config)
            elif has_layers:
                raise ValueError(
                    "OCI image-shaped attestation manifest has no config descriptor"
                )
            if has_layers:
                payload_descriptors = manifest.get("layers")
            else:
                payload_descriptors = manifest.get("blobs")
            if not isinstance(payload_descriptors, list):
                raise ValueError("attestation manifest has no layer/blob descriptors")
            for layer in payload_descriptors:
                if not isinstance(layer, dict):
                    raise ValueError("attestation payload descriptor must be an object")
                statement = _descriptor_json(archive, layer)
                predicate_type, statement_subjects = _statement_subjects(statement)
                predicate_types.add(predicate_type)
                subjects.update(statement_subjects)
                if image_digest in statement_subjects:
                    image_predicate_types.add(predicate_type)

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
