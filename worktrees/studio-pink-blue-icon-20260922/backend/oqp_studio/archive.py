"""Extract release archives without allowing files to escape the destination."""

from __future__ import annotations

import os
import re
import shutil
import stat
import tarfile
import zipfile
from pathlib import Path, PurePosixPath


class UnsafeArchiveError(ValueError):
    """The archive contains a path or file type that is unsafe to extract."""


def _destination(root: Path, member_name: str) -> Path:
    """Resolve one archive member below ``root`` or reject it."""
    normalized = member_name.replace("\\", "/")
    path = PurePosixPath(normalized)
    if (
        not normalized
        or path.is_absolute()
        or ".." in path.parts
        or re.match(r"^[A-Za-z]:", normalized)
    ):
        raise UnsafeArchiveError(f"unsafe archive path: {member_name}")
    root = root.resolve()
    destination = (root / Path(*path.parts)).resolve()
    if destination != root and root not in destination.parents:
        raise UnsafeArchiveError(f"archive path escapes its destination: {member_name}")
    return destination


def _link_target(root: Path, destination: Path, link_name: str) -> str:
    """Validate a relative archive link and return its normalized target."""
    normalized = link_name.replace("\\", "/")
    path = PurePosixPath(normalized)
    if (
        not normalized
        or path.is_absolute()
        or re.match(r"^[A-Za-z]:", normalized)
    ):
        raise UnsafeArchiveError(f"unsafe archive link target: {link_name}")
    root = root.resolve()
    target = (destination.parent / Path(*path.parts)).resolve()
    if target != root and root not in target.parents:
        raise UnsafeArchiveError(f"archive link escapes its destination: {link_name}")
    return normalized


def validate_links(root: Path) -> None:
    """Reject links that escape an extracted tree after it has been rearranged."""
    root = root.resolve()

    def reject_walk_error(error: OSError) -> None:
        raise UnsafeArchiveError(
            f"could not inspect archive directory: {error.filename or root}"
        ) from error

    for current, directory_names, file_names in os.walk(
        root, followlinks=False, onerror=reject_walk_error
    ):
        current_path = Path(current)
        for name in directory_names + file_names:
            candidate = current_path / name
            if candidate.is_symlink():
                _link_target(root, candidate, os.readlink(candidate))


def _create_link(destination: Path, link_name: str) -> None:
    """Create an archive link with the correct Windows link type when known."""
    target = destination.parent / Path(*PurePosixPath(link_name).parts)
    destination.symlink_to(link_name, target_is_directory=target.is_dir())


def extract_zip(archive: zipfile.ZipFile, target: Path) -> None:
    """Extract files, directories, and safe relative links from a ZIP archive."""
    members: list[tuple[zipfile.ZipInfo, Path, str | None]] = []
    for member in archive.infolist():
        destination = _destination(target, member.filename)
        mode = member.external_attr >> 16
        link_name = None
        if stat.S_ISLNK(mode):
            try:
                link_name = archive.read(member).decode("utf-8")
            except UnicodeDecodeError as error:
                raise UnsafeArchiveError(f"archive link is not UTF-8: {member.filename}") from error
            _link_target(target, destination, link_name)
        members.append((member, destination, link_name))

    for member, destination, link_name in members:
        if link_name is not None:
            continue
        if member.is_dir():
            destination.mkdir(parents=True, exist_ok=True)
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        with archive.open(member) as source, destination.open("wb") as output:
            shutil.copyfileobj(source, output)
        permissions = (member.external_attr >> 16) & 0o777
        if permissions:
            destination.chmod(permissions)

    for member, destination, link_name in members:
        if link_name is None:
            continue
        destination = _destination(target, member.filename)
        _link_target(target, destination, link_name)
        destination.parent.mkdir(parents=True, exist_ok=True)
        _create_link(destination, link_name)

    validate_links(target)


def extract_tar(archive: tarfile.TarFile, target: Path) -> None:
    """Extract files, directories, and safe relative links from a TAR archive."""
    members: list[tuple[tarfile.TarInfo, Path, str | None]] = []
    for member in archive.getmembers():
        destination = _destination(target, member.name)
        link_name = None
        if member.issym():
            link_name = _link_target(target, destination, member.linkname)
        elif not (member.isfile() or member.isdir()):
            raise UnsafeArchiveError(
                f"archive devices and hard links are not allowed: {member.name}"
            )
        members.append((member, destination, link_name))

    for member, destination, link_name in members:
        if link_name is not None:
            continue
        if member.isdir():
            destination.mkdir(parents=True, exist_ok=True)
            continue
        source = archive.extractfile(member)
        if source is None:
            raise UnsafeArchiveError(f"could not read archive member: {member.name}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with source, destination.open("wb") as output:
            shutil.copyfileobj(source, output)
        destination.chmod(member.mode & 0o777)

    # Links come last so no regular archive member can be written through one.
    for member, destination, link_name in members:
        if link_name is None:
            continue
        destination = _destination(target, member.name)
        _link_target(target, destination, link_name)
        destination.parent.mkdir(parents=True, exist_ok=True)
        _create_link(destination, link_name)

    validate_links(target)

    # Restrictive directory modes must come last or child extraction can lose
    # write/traverse access. Deepest-first also handles parents without +x.
    directories = [(member, destination) for member, destination, _ in members if member.isdir()]
    for member, destination in sorted(
        directories, key=lambda item: len(item[1].parts), reverse=True
    ):
        destination.chmod(member.mode & 0o777)
