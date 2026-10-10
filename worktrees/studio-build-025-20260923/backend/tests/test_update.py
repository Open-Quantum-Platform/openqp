import io
import tarfile

import pytest

from oqp_studio import update
from oqp_studio.archive import UnsafeArchiveError


def _write_tarball(path, members):
    with tarfile.open(path, mode="w:gz") as archive:
        for member, content in members:
            archive.addfile(member, io.BytesIO(content) if content is not None else None)


def test_macos_update_revalidates_links_against_the_app_bundle(tmp_path, monkeypatch):
    tarball = tmp_path / "update.tar.gz"
    payload = tarfile.TarInfo("payload")
    payload.size = len(b"outside")
    link = tarfile.TarInfo("OQP Studio.app/Contents/MacOS/openqp")
    link.type = tarfile.SYMTYPE
    link.linkname = "../../../payload"
    _write_tarball(tarball, [(payload, b"outside"), (link, None)])

    staging = tmp_path / "staging"
    staging.mkdir()
    monkeypatch.setattr(update, "installed_app_path", lambda: tmp_path / "installed.app")
    monkeypatch.setattr(update.tempfile, "mkdtemp", lambda **_kwargs: str(staging))

    with pytest.raises(UnsafeArchiveError):
        update.install_macos(tarball)


def test_macos_update_rejects_a_symlink_app_bundle(tmp_path, monkeypatch):
    tarball = tmp_path / "update.tar.gz"
    payload = tarfile.TarInfo("payload")
    payload.type = tarfile.DIRTYPE
    application = tarfile.TarInfo("OQP Studio.app")
    application.type = tarfile.SYMTYPE
    application.linkname = "payload"
    _write_tarball(tarball, [(payload, None), (application, None)])

    staging = tmp_path / "staging"
    staging.mkdir()
    monkeypatch.setattr(update, "installed_app_path", lambda: tmp_path / "installed.app")
    monkeypatch.setattr(update.tempfile, "mkdtemp", lambda **_kwargs: str(staging))

    with pytest.raises(RuntimeError, match="one real application bundle"):
        update.install_macos(tarball)
