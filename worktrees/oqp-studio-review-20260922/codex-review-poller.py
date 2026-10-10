#!/usr/bin/python3
from __future__ import annotations

import fcntl
import json
import os
from pathlib import Path
import pwd
import re
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request


GITLAB_URL = "https://qchemlab.knu.ac.kr"
PROJECTS = {
    17: "open-quantum-platform/openqp",
    19: "open-quantum-platform/internal/openqp",
    18: "open-quantum-platform/oqp-studio",
}
TOKEN_FILE = Path("/var/lib/gitlab-review/gitlab.token")
MIRROR_ROOT = Path("/var/lib/gitlab-review/mirror")
WORK_ROOT = Path("/var/lib/codex-review/work")
STATE_FILE = Path("/var/lib/codex-review/state.json")
LOCK_FILE = Path("/run/codex-review-poller.lock")
CODEX = "/var/lib/codex-review/.local/bin/codex"
SAFE_REF = re.compile(r"[A-Za-z0-9._/-]+")


def api(method: str, path: str, data: dict[str, str] | None = None):
    token = TOKEN_FILE.read_text(encoding="utf-8").strip()
    body = None
    headers = {"PRIVATE-TOKEN": token, "User-Agent": "openqp-codex-review-runner/1"}
    if data is not None:
        body = urllib.parse.urlencode(data).encode("utf-8")
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    request = urllib.request.Request(GITLAB_URL + path, data=body, headers=headers, method=method)
    with urllib.request.urlopen(request, timeout=30) as response:
        payload = response.read()
    return json.loads(payload) if payload else None


def run(
    command: list[str],
    *,
    timeout: int = 600,
    capture: bool = False,
    user: str | None = None,
) -> subprocess.CompletedProcess[str]:
    drop_privileges = None
    if user is not None:
        account = pwd.getpwnam(user)

        def drop_privileges() -> None:
            os.initgroups(user, account.pw_gid)
            os.setgid(account.pw_gid)
            os.setuid(account.pw_uid)

    return subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
        timeout=timeout,
        preexec_fn=drop_privileges,
    )


def as_user(user: str, command: list[str], *, timeout: int = 600, capture: bool = False):
    home = "/var/lib/gitlab-review" if user == "gitlab-review" else "/var/lib/codex-review"
    return run(
        ["/usr/bin/env", f"HOME={home}", *command],
        timeout=timeout,
        capture=capture,
        user=user,
    )


def load_state() -> dict[str, str]:
    if not STATE_FILE.exists():
        return {}
    return json.loads(STATE_FILE.read_text(encoding="utf-8"))


def save_state(state: dict[str, str]) -> None:
    temporary = STATE_FILE.with_suffix(".tmp")
    temporary.write_text(json.dumps(state, sort_keys=True) + "\n", encoding="utf-8")
    os.chmod(temporary, 0o600)
    os.replace(temporary, STATE_FILE)


def existing_review(project_id: int, iid: int, sha: str) -> bool:
    marker = f"<!-- codex-review-runner sha:{sha} -->"
    notes = api("GET", f"/api/v4/projects/{project_id}/merge_requests/{iid}/notes?per_page=100&sort=desc")
    return any(marker in str(note.get("body", "")) for note in notes)


def add_eyes(project_id: int, iid: int) -> None:
    try:
        api("POST", f"/api/v4/projects/{project_id}/merge_requests/{iid}/award_emoji", {"name": "eyes"})
    except urllib.error.HTTPError as error:
        if error.code not in {404, 409}:
            raise


def mirror_path(project_id: int) -> Path:
    return MIRROR_ROOT / f"project-{project_id}.git"


def ensure_mirror(project_id: int, project_path: str) -> Path:
    mirror = mirror_path(project_id)
    if mirror.exists():
        return mirror
    as_user("gitlab-review", ["/usr/bin/git", "init", "--bare", str(mirror)])
    as_user(
        "gitlab-review",
        [
            "/usr/bin/git",
            "-C",
            str(mirror),
            "remote",
            "add",
            "origin",
            f"https://codex-review-bot@qchemlab.knu.ac.kr/{project_path}.git",
        ],
    )
    return mirror


def fetch_and_checkout(
    project_id: int,
    project_path: str,
    iid: int,
    sha: str,
    target_branch: str,
) -> tuple[Path, str]:
    if not SAFE_REF.fullmatch(target_branch) or ".." in target_branch:
        raise ValueError(f"unsafe target branch: {target_branch!r}")
    mirror = ensure_mirror(project_id, project_path)
    target_ref = f"target-{iid}"
    source_ref = f"mr-{iid}"
    git_env = [
        "/usr/bin/env",
        "GIT_ASKPASS=/usr/local/libexec/codex-review-askpass",
        "GIT_TERMINAL_PROMPT=0",
        "/usr/bin/git",
        "-C",
        str(mirror),
        "fetch",
        "--force",
        "--prune",
        "origin",
        f"+refs/heads/{target_branch}:refs/heads/{target_ref}",
        f"+refs/merge-requests/{iid}/head:refs/heads/{source_ref}",
    ]
    as_user("gitlab-review", git_env, timeout=600)

    work = WORK_ROOT / f"project-{project_id}-mr-{iid}"
    if work.exists():
        as_user("codex-review", ["/usr/bin/rm", "-rf", "--", str(work)], timeout=120)
    as_user(
        "codex-review",
        ["/usr/bin/git", "-c", f"safe.directory={mirror}", "clone", "--shared", "--no-checkout", str(mirror), str(work)],
        timeout=600,
    )
    as_user(
        "codex-review",
        ["/usr/bin/git", "-C", str(work), "checkout", "--detach", sha],
        timeout=120,
    )
    return work, f"origin/{target_ref}"


def review(work: Path, base: str, project_path: str, iid: int, title: str) -> str:
    result = as_user(
        "codex-review",
        [
            "/usr/bin/env",
            "CODEX_HOME=/var/lib/codex-review/.codex",
            CODEX,
            "exec",
            "--ephemeral",
            "--ignore-user-config",
            "--sandbox",
            "read-only",
            "--color",
            "never",
            "-C",
            str(work),
            "review",
            "--base",
            base,
            "--title",
            f"GitLab {project_path} MR !{iid}: {title}",
        ],
        timeout=1800,
        capture=True,
    )
    output = (result.stdout or "").strip()
    if not output:
        raise RuntimeError(f"Codex produced no review output; stderr={result.stderr[-2000:]}")
    return output[-60000:]


def post_review(project_id: int, mr: dict, output: str) -> None:
    iid = int(mr["iid"])
    sha = str(mr["sha"])
    body = (
        "## Codex automated review\n\n"
        f"Reviewed commit `{sha[:12]}` with the official Codex CLI.\n\n"
        f"{output}\n\n"
        f"<!-- codex-review-runner sha:{sha} -->"
    )
    api("POST", f"/api/v4/projects/{project_id}/merge_requests/{iid}/notes", {"body": body})


def main() -> int:
    LOCK_FILE.parent.mkdir(parents=True, exist_ok=True)
    with LOCK_FILE.open("w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print("another review poll is active")
            return 0

        state = load_state()
        for project_id, project_path in PROJECTS.items():
            merge_requests = api(
                "GET",
                f"/api/v4/projects/{project_id}/merge_requests?state=opened&per_page=100",
            )
            for mr in merge_requests:
                iid = int(mr["iid"])
                sha = str(mr["sha"])
                key = f"{project_id}:{iid}"
                if state.get(key) == sha or existing_review(project_id, iid, sha):
                    state[key] = sha
                    continue
                print(f"reviewing {project_path}!{iid} {sha}", flush=True)
                add_eyes(project_id, iid)
                work, base = fetch_and_checkout(
                    project_id,
                    project_path,
                    iid,
                    sha,
                    str(mr["target_branch"]),
                )
                output = review(work, base, project_path, iid, str(mr["title"]))
                post_review(project_id, mr, output)
                state[key] = sha
                save_state(state)
                print(f"posted review for {project_path}!{iid} {sha}", flush=True)
        save_state(state)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"codex review poll failed: {type(error).__name__}: {error}", file=sys.stderr)
        raise
