"""Run the production sync shell against an offline GitLab/Git simulator."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github/workflows/sync-gitlab.yml"
SHA = "a" * 40

CURL = r'''#!/usr/bin/env python3
import json, os, pathlib, sys
args = sys.argv[1:]
def option(name, default=None):
    return args[args.index(name) + 1] if name in args else default
state_path = pathlib.Path(os.environ["SIM_STATE"])
state = json.loads(state_path.read_text()) if state_path.exists() else {}
case = os.environ["SIM_CASE"]
url = args[-1]
method = option("--request", "GET")
code, body, transport_exit = 200, {}, 0
if "/repository/branches/" in url:
    key = "branch_get"
    count = state.get(key, 0) + 1
    if case == "branch_transport_exhausted": transport_exit = 28
    elif case.startswith("branch_transport_"):
        transport_exit = int(case.rsplit("_", 1)[1]) if count < 3 else 0
    elif case == "branch_certificate": transport_exit = 60
    if case == "branch_401": code = 401
    elif case == "branch_missing" or (case == "branch_delay" and count < 3): code = 404
    else: body = {"commit": {"id": "b" * 40 if case == "wrong_sha" else "a" * 40}}
elif method == "POST":
    key = "post"
    count = state.get(key, 0) + 1
    if case in ("mr_missing", "mr_delay", "lost_response") and (case != "mr_delay" or count == 1):
        code, body = 400, {"message": {"source_branch": ["does not exist"]}}
    elif case == "invalid_400": code, body = 400, {"message": {"target_branch": ["does not exist"]}}
    else: code, body = 201, {"iid": 39}
elif method == "PUT": key = "merge"
else:
    key = "list"
    body = [{"iid": 39}] if case == "lost_response" and state.get("post", 0) else []
state[key] = state.get(key, 0) + 1
state_path.write_text(json.dumps(state))
if transport_exit:
    if "--write-out" in args: print("000", end="")
    print("simulated curl transport failure", file=sys.stderr)
    sys.exit(transport_exit)
pathlib.Path(option("--output")).write_text(json.dumps(body))
if "--write-out" in args: print(code, end="")
'''

class GitLabSyncRecoveryTests(unittest.TestCase):
    def run_sync(self, case, workflow=WORKFLOW):
        # Exercise the actual workflow commands, not a copied retry algorithm.
        text = workflow.read_text()
        script = text.split("        run: |\n", 1)[1]
        script = "\n".join(line[10:] if line.startswith("          ") else line
                           for line in script.splitlines())
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, body in {
                "curl": CURL,
                "git": '#!/bin/sh\ncase " $* " in *" merge-base "*) exit 1;; esac\nexit 0\n',
                "sleep": '#!/bin/sh\nexit 0\n',
            }.items():
                path = root / name
                path.write_text(body)
                path.chmod(0o755)
            # jq is part of the GitHub runner contract; use the real parser.
            self.assertIsNotNone(shutil.which("jq"))
            state = root / "state.json"
            env = dict(os.environ, PATH=str(root) + os.pathsep + os.environ["PATH"],
                       GITHUB_SHA=SHA, GITLAB_TOKEN="offline-test-token",
                       GITLAB_PROJECT_URL="https://invalid.example/openqp.git",
                       GITLAB_API_URL="https://invalid.example/api/v4",
                       GITLAB_PROJECT_ID="17", SIM_CASE=case, SIM_STATE=str(state))
            # macOS base64 lacks -w0; emulate only this setup command in the
            # simulator. The production workflow runs on Ubuntu.
            script = script.replace('base64 -w0', 'base64')
            result = subprocess.run(["bash", "-c", script], env=env,
                                    capture_output=True, text=True, timeout=20)
            return result, json.loads(state.read_text())

    def test_branch_api_visibility_delay(self):
        result, calls = self.run_sync("branch_delay")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(calls, {"branch_get": 3, "list": 1, "post": 1, "merge": 1})

    def test_branch_transport_errors_recover_without_corrupting_http_status(self):
        for code in (6, 7, 28, 35):
            with self.subTest(curl_exit=code):
                result, calls = self.run_sync(f"branch_transport_{code}")
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(calls, {"branch_get": 3, "list": 1, "post": 1, "merge": 1})

    def test_branch_transport_exhaustion_and_certificate_failure_stop_before_mr(self):
        for case, expected in (("branch_transport_exhausted", 3), ("branch_certificate", 1)):
            with self.subTest(case=case):
                result, calls = self.run_sync(case)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(calls, {"branch_get": expected})

    def test_mr_creation_visibility_delay(self):
        result, calls = self.run_sync("mr_delay")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(calls["post"], 2)
        self.assertEqual(calls["list"], 2)
        self.assertEqual(calls["merge"], 1)

    def test_recovery_finds_existing_mr_before_reposting(self):
        result, calls = self.run_sync("lost_response")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(calls["post"], 1)
        self.assertEqual(calls["merge"], 1)

    def test_permanent_branch_errors_do_not_create_mr(self):
        for case, expected in (("branch_401", 1), ("wrong_sha", 1), ("branch_missing", 3)):
            with self.subTest(case=case):
                result, calls = self.run_sync(case)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(calls, {"branch_get": expected})

    def test_other_400_is_not_retried(self):
        result, calls = self.run_sync("invalid_400")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(calls["post"], 1)
        self.assertNotIn("merge", calls)

    def test_visibility_retries_are_bounded(self):
        result, calls = self.run_sync("mr_missing")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(calls["post"], 3)
        self.assertEqual(calls["list"], 4)
        self.assertNotIn("merge", calls)

if __name__ == "__main__":
    unittest.main()
