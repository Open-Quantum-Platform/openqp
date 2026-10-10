"""TLS settings for the app's outbound requests.

Campus and corporate networks routinely inspect TLS by re-signing traffic
with a private root. Three things can make that work, in descending order of
preference, and this module implements all three:

1. Verify against the operating system's trust store, where the network's
   root is already installed (otherwise browsers would fail too).
2. Verify against a CA bundle the user points at, for a root that was handed
   out as a file but never installed system-wide.
3. Skip verification, which the user has to switch on deliberately.
"""

from __future__ import annotations

import json
import os
import ssl
import sys
from pathlib import Path


def settings_path() -> Path:
    """Where the settings live, following each platform's convention."""
    override = os.environ.get("OQP_STUDIO_CONFIG")
    if override:
        return Path(override)
    if sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support" / "OQP Studio"
    elif os.name == "nt":
        base = Path(os.environ.get("APPDATA", Path.home())) / "OQP Studio"
    else:
        base = Path(os.environ.get("XDG_CONFIG_HOME",
                                   Path.home() / ".config")) / "oqp-studio"
    return base / "network.json"


DEFAULTS = {"ca_bundle": "", "insecure": False}


def load() -> dict:
    """Stored settings, with environment variables taking precedence."""
    values = dict(DEFAULTS)
    try:
        stored = json.loads(settings_path().read_text())
        if isinstance(stored, dict):
            values["ca_bundle"] = str(stored.get("ca_bundle") or "")
            values["insecure"] = bool(stored.get("insecure"))
    except (OSError, ValueError):
        pass
    # An environment variable wins, so a site can configure this centrally.
    env_bundle = os.environ.get("OQP_STUDIO_CA_BUNDLE") or os.environ.get("SSL_CERT_FILE")
    if env_bundle:
        values["ca_bundle"] = env_bundle
    if os.environ.get("OQP_STUDIO_INSECURE_SSL") == "1":
        values["insecure"] = True
    return values


def save(ca_bundle: str = "", insecure: bool = False) -> dict:
    path = settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"ca_bundle": ca_bundle, "insecure": bool(insecure)},
                               indent=2))
    return load()


def use_system_trust_store() -> bool:
    """Verify TLS against the OS trust store. True when that succeeded.

    A frozen build ships certifi's CA list, which does not include private
    roots; the OS keychain does.
    """
    try:
        import truststore
    except ImportError:
        return False
    try:
        truststore.inject_into_ssl()
    except Exception:  # noqa: BLE001 — never let this stop the server starting
        return False
    return True


def context() -> ssl.SSLContext | None:
    """The SSL context for outbound requests, or None for the default one."""
    values = load()
    if values["insecure"]:
        # Deliberately unverified: the user asked for this after their network
        # could not be verified any other way.
        unverified = ssl._create_unverified_context()
        return unverified
    bundle = values["ca_bundle"]
    if bundle and Path(bundle).is_file():
        return ssl.create_default_context(cafile=bundle)
    return None


def status() -> dict:
    """What the app is doing about TLS right now, for the settings dialog."""
    values = load()
    bundle = values["ca_bundle"]
    return {
        "system_trust_store": _system_trust_active,
        "ca_bundle": bundle,
        "ca_bundle_found": bool(bundle) and Path(bundle).is_file(),
        "insecure": values["insecure"],
        "settings_path": str(settings_path()),
    }


_system_trust_active = False


def activate() -> None:
    global _system_trust_active
    _system_trust_active = use_system_trust_store()
