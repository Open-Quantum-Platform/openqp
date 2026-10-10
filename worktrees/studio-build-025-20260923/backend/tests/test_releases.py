import copy
import hashlib
import io

import pytest

from oqp_studio import releases, update


def manifest():
    name = "OQP-Studio-0.2.4-linux-x86_64-with-engine.deb"
    return {"schema_version": 1, "version": "0.2.4", "channel": "stable",
            "studio_commit": "a" * 40, "engine": {"commit": "b" * 40},
            "assets": [{"name": name, "url": releases.ASSETS + "v0.2.4/" + name,
                        "sha256": hashlib.sha256(b"installer").hexdigest(), "size": 9,
                        "kind": "installer", "variant": "with-engine"}]}


def test_linux_update_preserves_offline_variant(monkeypatch):
    monkeypatch.setattr(update, "_platform_asset", lambda: ("linux-x86_64", ".AppImage"))
    assets = manifest()["assets"] + [{"name": "OQP-Studio-0.2.4-linux-x86_64.AppImage"}]
    assert update.pick_asset(assets, True)["name"].endswith("with-engine.deb")
    assert update.pick_asset(assets, False)["name"].endswith(".AppImage")
    assert update.pick_asset(assets[1:], True) is None


@pytest.mark.parametrize("field,value", [
    ("name", "../../payload"), ("sha256", ""), ("size", 0),
    ("url", "https://github.com/other/repository/releases/download/v0.2.4/payload"),
])
def test_rejects_unbound_or_unverifiable_asset(field, value):
    data = manifest()
    data["assets"][0][field] = value
    with pytest.raises(ValueError):
        releases.validate(data, "v0.2.4")


def test_rejects_manifest_from_different_studio_version():
    with pytest.raises(ValueError):
        releases.validate(manifest(), "v0.2.5")


def test_engine_download_uses_installed_version_not_latest(monkeypatch):
    requested = []
    monkeypatch.setattr(releases, "__version__", "0.2.4")
    monkeypatch.setattr(releases, "read_json", lambda url: requested.append(url) or {})
    monkeypatch.setattr(releases, "manifest_for", lambda r: manifest())
    releases.installed_manifest()
    assert requested == [releases.API + "/tags/v0.2.4"]


def test_corrupt_download_is_deleted_before_installation(tmp_path, monkeypatch):
    monkeypatch.setattr(releases.urllib.request, "urlopen", lambda *a, **k: io.BytesIO(b"corrupted"))
    target = tmp_path / "installer"
    with pytest.raises(ValueError, match="checksum"):
        releases.download(manifest()["assets"][0], target)
    assert not target.exists()


def test_download_validates_size_and_sha256(tmp_path, monkeypatch):
    monkeypatch.setattr(releases.urllib.request, "urlopen", lambda *a, **k: io.BytesIO(b"installer"))
    assert releases.download(manifest()["assets"][0], tmp_path / "ok").read_bytes() == b"installer"


def test_stable_is_default_and_preview_requires_explicit_choice(tmp_path, monkeypatch):
    monkeypatch.setattr(releases.network, "settings_path", lambda: tmp_path / "network.json")
    assert releases.channel() == "stable"
    releases.set_channel("preview")
    assert releases.channel() == "preview"
    with pytest.raises(ValueError):
        releases.set_channel("nightly")


def test_preview_cannot_be_claimed_as_stable(monkeypatch):
    data = copy.deepcopy(manifest())
    data["channel"] = "preview"
    monkeypatch.setattr(releases, "read_json", lambda url: data)
    with pytest.raises(ValueError, match="channel"):
        releases.manifest_for({"tag_name": "v0.2.4", "prerelease": False,
                              "assets": [{"browser_download_url": releases.ASSETS + "v0.2.4/manifest.json"}]})


@pytest.mark.parametrize('current,newer', [
    ('0.2.5-rc.1', '0.2.5-rc.2'), ('0.2.5-rc.9', '0.2.5-rc.10'),
    ('0.2.5-rc.10', '0.2.5'), ('0.2.4', '0.2.5-rc.1'),
])
def test_preview_and_stable_version_order(current, newer):
    assert releases.version_key(newer) > releases.version_key(current)


@pytest.mark.parametrize('text,expected', [
    ('nmr(gauge=giao,acid=true)', True),
    ('nmr(gauge=giao, acid=False)', False),
    ('[properties]\nscf_prop=nmr,acid', True),
    ('# nmr(gauge=giao,acid=true)\nenergy', False),
])
def test_acid_detection_covers_ui_and_legacy_input(text, expected):
    from oqp_studio import engine
    assert engine.requires_acid(text) == expected


def test_acid_rejects_unverified_remote_engine():
    from oqp_studio import engine
    with pytest.raises(ValueError, match='qualified Studio engine'):
        engine.check_capabilities('nmr(gauge=giao,acid=true)', 'ssh')
