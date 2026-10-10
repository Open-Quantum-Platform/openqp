import io
import os
import stat
import tarfile
import zipfile

import pytest

from oqp_studio.archive import UnsafeArchiveError, extract_tar, extract_zip, validate_links
from oqp_studio.engine import _strip_single_root


def test_zip_extracts_regular_file_and_preserves_executable_bit(tmp_path):
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w") as archive:
        member = zipfile.ZipInfo("openqp/openqp")
        member.external_attr = 0o755 << 16
        archive.writestr(member, b"engine")

    with zipfile.ZipFile(payload) as archive:
        extract_zip(archive, tmp_path)

    executable = tmp_path / "openqp" / "openqp"
    assert executable.read_bytes() == b"engine"
    if os.name != "nt":
        # Windows has no POSIX executable bit; engine.install handles its .exe
        # directly and only Unix-like platforms need the restored archive mode.
        assert executable.stat().st_mode & stat.S_IXUSR


@pytest.mark.parametrize("name", ["../outside", "/absolute", "C:\\outside", "C:outside"])
def test_zip_rejects_paths_outside_destination(tmp_path, name):
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w") as archive:
        archive.writestr(name, b"escape")

    with zipfile.ZipFile(payload) as archive, pytest.raises(UnsafeArchiveError):
        extract_zip(archive, tmp_path)

    assert not (tmp_path.parent / "outside").exists()


def test_zip_extracts_safe_relative_symlinks(tmp_path):
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w") as archive:
        archive.writestr("app/lib/real.dylib", b"library")
        link = zipfile.ZipInfo("app/lib/alias.dylib")
        link.external_attr = (stat.S_IFLNK | 0o777) << 16
        archive.writestr(link, b"real.dylib")

    with zipfile.ZipFile(payload) as archive:
        extract_zip(archive, tmp_path)

    link = tmp_path / "app" / "lib" / "alias.dylib"
    assert link.is_symlink()
    assert link.read_bytes() == b"library"


@pytest.mark.skipif(
    os.name == "nt",
    reason="Windows normalizes '..' before following directory symlinks",
)
def test_zip_revalidates_the_completed_symlink_graph(tmp_path):
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w") as archive:
        first = zipfile.ZipInfo("app/openqp")
        first.external_attr = (stat.S_IFLNK | 0o777) << 16
        archive.writestr(first, b"sub/link/../../outside")
        later = zipfile.ZipInfo("app/sub/link")
        later.external_attr = (stat.S_IFLNK | 0o777) << 16
        archive.writestr(later, b"..")

    with zipfile.ZipFile(payload) as archive, pytest.raises(UnsafeArchiveError):
        extract_zip(archive, tmp_path)


def test_tar_rejects_symlinks_that_escape_the_destination(tmp_path):
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w") as archive:
        member = tarfile.TarInfo("openqp/link")
        member.type = tarfile.SYMTYPE
        member.linkname = "../../outside"
        archive.addfile(member)
    payload.seek(0)

    with tarfile.open(fileobj=payload) as archive, pytest.raises(UnsafeArchiveError):
        extract_tar(archive, tmp_path)


def test_tar_extracts_safe_relative_symlinks(tmp_path):
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w") as archive:
        target = tarfile.TarInfo("app/lib/real.dylib")
        target.size = len(b"library")
        archive.addfile(target, io.BytesIO(b"library"))
        link = tarfile.TarInfo("app/lib/alias.dylib")
        link.type = tarfile.SYMTYPE
        link.linkname = "real.dylib"
        archive.addfile(link)
    payload.seek(0)

    with tarfile.open(fileobj=payload) as archive:
        extract_tar(archive, tmp_path)

    link = tmp_path / "app" / "lib" / "alias.dylib"
    assert link.is_symlink()
    assert link.read_bytes() == b"library"


@pytest.mark.skipif(
    os.name == "nt",
    reason="Windows normalizes '..' before following directory symlinks",
)
def test_tar_revalidates_the_completed_symlink_graph(tmp_path):
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w") as archive:
        first = tarfile.TarInfo("app/openqp")
        first.type = tarfile.SYMTYPE
        first.linkname = "sub/link/../../outside"
        archive.addfile(first)
        later = tarfile.TarInfo("app/sub/link")
        later.type = tarfile.SYMTYPE
        later.linkname = ".."
        archive.addfile(later)
    payload.seek(0)

    with tarfile.open(fileobj=payload) as archive, pytest.raises(UnsafeArchiveError):
        extract_tar(archive, tmp_path)


def test_tar_extracts_regular_files(tmp_path):
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w") as archive:
        member = tarfile.TarInfo("openqp/README.txt")
        member.size = len(b"OpenQP")
        member.mode = 0o644
        archive.addfile(member, io.BytesIO(b"OpenQP"))
    payload.seek(0)

    with tarfile.open(fileobj=payload) as archive:
        extract_tar(archive, tmp_path)

    assert (tmp_path / "openqp" / "README.txt").read_text() == "OpenQP"


def test_tar_defers_restrictive_directory_modes_until_children_exist(tmp_path):
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w") as archive:
        directory = tarfile.TarInfo("openqp/lib")
        directory.type = tarfile.DIRTYPE
        directory.mode = 0o555
        archive.addfile(directory)
        child = tarfile.TarInfo("openqp/lib/library.dat")
        child.size = len(b"library")
        archive.addfile(child, io.BytesIO(b"library"))
    payload.seek(0)

    with tarfile.open(fileobj=payload) as archive:
        extract_tar(archive, tmp_path)

    directory = tmp_path / "openqp" / "lib"
    assert (directory / "library.dat").read_bytes() == b"library"
    if os.name != "nt":
        assert stat.S_IMODE(directory.stat().st_mode) == 0o555


@pytest.mark.skipif(os.name == "nt", reason="Windows does not enforce POSIX directory modes")
def test_link_validation_rejects_an_unreadable_directory_after_flattening(tmp_path):
    staging = tmp_path / "staging"
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w") as archive:
        directory = tarfile.TarInfo("wrapper/locked")
        directory.type = tarfile.DIRTYPE
        directory.mode = 0o300
        archive.addfile(directory)
        link = tarfile.TarInfo("wrapper/locked/link")
        link.type = tarfile.SYMTYPE
        link.linkname = "../../settings.json"
        archive.addfile(link)
    payload.seek(0)

    with tarfile.open(fileobj=payload) as archive:
        extract_tar(archive, staging)
    _strip_single_root(staging)

    locked = staging / "locked"
    try:
        with pytest.raises(UnsafeArchiveError):
            validate_links(staging)
    finally:
        locked.chmod(0o700)


def test_single_symlink_wrapper_is_not_flattened(tmp_path):
    staging = tmp_path / "staging"
    staging.mkdir()
    wrapper = staging / "wrapper"
    wrapper.symlink_to(".", target_is_directory=True)

    _strip_single_root(staging)

    assert wrapper.is_symlink()
    assert not (tmp_path / "staging.unpack").exists()


def test_links_are_revalidated_after_archive_root_is_flattened(tmp_path):
    staging = tmp_path / "staging"
    (tmp_path / "settings.json").write_text("outside")
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w") as archive:
        link = tarfile.TarInfo("wrapper/openqp")
        link.type = tarfile.SYMTYPE
        link.linkname = "../settings.json"
        archive.addfile(link)
    payload.seek(0)

    with tarfile.open(fileobj=payload) as archive:
        extract_tar(archive, staging)
    validate_links(staging)

    _strip_single_root(staging)
    assert (staging / "openqp").is_symlink()
    with pytest.raises(UnsafeArchiveError):
        validate_links(staging)
