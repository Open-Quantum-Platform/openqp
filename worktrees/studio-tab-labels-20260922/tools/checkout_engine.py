#!/usr/bin/env python3
"""Resolve the public gateway once; check out the same engine SHA on every OS."""
from __future__ import annotations

import argparse
import base64
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
GATEWAY = "https://qchemlab.knu.ac.kr/open-quantum-platform/openqp.git"
SSH_GATEWAY = "gitlab-qchemlab:open-quantum-platform/openqp.git"


def configuration() -> dict:
    config = json.loads((ROOT / "engine-upstream.json").read_text())
    if config != {"repository": GATEWAY, "tracking_branch": "main"}:
        raise ValueError("engine upstream must be the public OpenQP gateway main")
    return config


def git_env(ssh: bool = False) -> dict:
    env = os.environ.copy()
    env["GIT_TERMINAL_PROMPT"] = "0"
    # No credentials in argv, output, persistent git config, or the manifest.
    token = env.get("CI_JOB_TOKEN") or env.get("OPENQP_GATEWAY_TOKEN")
    if token and not ssh:
        user = "gitlab-ci-token" if env.get("CI_JOB_TOKEN") else "oauth2"
        encoded = base64.b64encode(f"{user}:{token}".encode()).decode()
        index = int(env.get("GIT_CONFIG_COUNT", "0"))
        env["GIT_CONFIG_COUNT"] = str(index + 1)
        env[f"GIT_CONFIG_KEY_{index}"] = "http.https://qchemlab.knu.ac.kr/.extraHeader"
        env[f"GIT_CONFIG_VALUE_{index}"] = f"Authorization: Basic {encoded}"
    return env


def checkout(destination: Path, ref: str, *, ssh: bool = False) -> dict:
    config = configuration()
    # Tags, branch names, or full SHAs, never git options or revision expressions.
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._/-]*", ref) or ".." in ref:
        raise ValueError("invalid engine ref")
    if destination.exists():
        raise ValueError("engine destination must not already exist")
    destination.mkdir(parents=True)
    env = git_env(ssh)

    def git(*args: str) -> str:
        return subprocess.check_output(
            ["git", "-C", str(destination), *args], env=env, text=True
        ).strip()

    git("init", "--quiet")
    git("remote", "add", "origin", SSH_GATEWAY if ssh else config["repository"])
    # Full commit ancestry allows checking that an explicit ref is public main
    # history. Blob filtering avoids downloading every historical engine tree.
    git("fetch", "--quiet", "--filter=blob:none", "--no-tags", "origin",
        "+refs/heads/main:refs/remotes/origin/main")
    if ref == "main":
        commit = git("rev-parse", "refs/remotes/origin/main^{commit}")
    else:
        git("fetch", "--quiet", "--filter=blob:none", "--no-tags", "origin", ref)
        commit = git("rev-parse", "FETCH_HEAD^{commit}")
    git("merge-base", "--is-ancestor", commit, "refs/remotes/origin/main")
    git("checkout", "--quiet", "--detach", commit)
    manifest = {**config, "requested_ref": ref, "commit": commit}
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref", default="main")
    parser.add_argument("--destination", type=Path, default=Path("openqp"))
    parser.add_argument("--manifest", type=Path, default=Path("engine-source.json"))
    parser.add_argument("--ssh", action="store_true", help="Use Ultra's verified gitlab-qchemlab alias")
    args = parser.parse_args()
    manifest = checkout(args.destination, args.ref or "main", ssh=args.ssh)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    if output := os.environ.get("GITHUB_OUTPUT"):
        with open(output, "a") as stream:
            stream.write(f"sha={manifest['commit']}\ntag={manifest['commit']}\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, subprocess.CalledProcessError) as error:
        print(f"Engine checkout refused: {error}", file=sys.stderr)
        sys.exit(1)
