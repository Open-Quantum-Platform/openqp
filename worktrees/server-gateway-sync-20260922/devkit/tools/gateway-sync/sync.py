#!/usr/bin/python3
"""Inbound-only GitLab sync; never checks out or executes fetched repository code."""
from __future__ import annotations
import argparse
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from openqp_sync_policy import sync_gateway, INITIAL_SHA, INITIAL_DOC, INITIAL_DOC_BLOB

BASE = 'https://qchemlab.knu.ac.kr'
SOURCE = 17
TARGET = 19
PATHS = {17: 'open-quantum-platform/openqp', 19: 'open-quantum-platform/internal/openqp'}
PREFIX = 'gateway-sync/'
ADOPT_IID = 4
ADOPT_BRANCH = 'codex/gateway-internal-sync-20260922'
SHA = re.compile(r'^[0-9a-f]{40}$')
REQUIRED_JOBS = {'source-policy', 'linux-openblas64-full', 'linux-openblas64-mpi',
                 'macos-intel-accelerate-ilp64', 'macos-arm64-accelerate-ilp64',
                 'windows-intel-mkl-ilp64', 'windows-source-regression'}

def pipeline_gate(mr, jobs, project=TARGET):
    p = mr.get('head_pipeline') or {}
    if (p.get('sha') != mr['sha'] or p.get('project_id') != project
            or p.get('source') != 'merge_request_event' or p.get('status') != 'success'):
        return 'waiting_successful_ci'
    if not REQUIRED_JOBS.issubset({j['name'] for j in jobs}):
        return 'blocked_missing_ci_jobs'
    if any(j['status'] != 'success' for j in jobs):
        return 'blocked_incomplete_ci_jobs'
    return 'clean'


class Sync:
    def __init__(self, state: Path, dry_run=False):
        self.state = state
        self.repo = state / 'repo.git'
        self.token_file = state / 'gitlab.token'
        self.dry_run = dry_run
        self.env = dict(os.environ, GIT_ASKPASS='/usr/local/libexec/openqp-inbound-askpass',
                        GIT_TERMINAL_PROMPT='0', GIT_CONFIG_NOSYSTEM='1',
                        GIT_CONFIG_GLOBAL='/dev/null', GIT_AUTHOR_NAME='OpenQP Inbound Sync',
                        GIT_AUTHOR_EMAIL='openqp-inbound-sync@noreply.qchemlab.knu.ac.kr',
                        GIT_COMMITTER_NAME='OpenQP Inbound Sync',
                        GIT_COMMITTER_EMAIL='openqp-inbound-sync@noreply.qchemlab.knu.ac.kr')

    def api(self, method, path, data=None, token_project=TARGET):
        # Writes are constrained in code to the private target and its token rotation.
        if method != 'GET' and not (path.startswith(f'/projects/{TARGET}/') or
                (token_project == SOURCE and path == f'/projects/{SOURCE}/access_tokens/self/rotate')):
            raise ValueError('Attempted write outside private project')
        request = urllib.request.Request(BASE + '/api/v4' + path, method=method,
            headers={'PRIVATE-TOKEN': (self.state / ('gateway.token' if token_project == SOURCE
                                                   else 'gitlab.token')).read_text().strip(),
                     'Content-Type': 'application/json'},
            data=json.dumps(data).encode() if data is not None else None)
        with urllib.request.urlopen(request, timeout=45) as response:
            payload = response.read()
            return json.loads(payload) if payload else None

    def pages(self, path, token_project=TARGET):
        items = []
        sep = '&' if '?' in path else '?'
        for page in range(1, 101):
            batch = self.api('GET', f'{path}{sep}per_page=100&page={page}', token_project=token_project)
            items.extend(batch)
            if len(batch) < 100:
                return items
        raise RuntimeError('Pagination limit reached; refusing partial evidence')

    def status(self, status, **fields):
        data = {'time': dt.datetime.now(dt.timezone.utc).isoformat(), 'status': status, **fields}
        tmp = self.state / 'status.tmp'
        tmp.write_text(json.dumps(data, sort_keys=True) + '\n')
        tmp.replace(self.state / 'status.json')
        print(json.dumps(data, sort_keys=True), flush=True)
        return data

    def git(self, *args, check=True, stdin=None, gateway=False):
        result = subprocess.run(['/usr/bin/git', '-c', 'core.hooksPath=/dev/null',
            '-c', 'protocol.file.allow=never', '--git-dir', str(self.repo), *args],
            env=dict(self.env, OPENQP_SYNC_TOKEN=str(self.state /
                ('gateway.token' if gateway else 'gitlab.token'))),
            text=True, input=stdin, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=300)
        if check and result.returncode:
            # Do not log git stderr, URLs with secrets, credential helpers or payloads.
            raise RuntimeError(f'git {args[0]} failed with exit {result.returncode}')
        return result

    def fetch(self, project, branch, ref):
        self.git('fetch', '--no-tags', f'{BASE}/{PATHS[project]}.git',
                 f'refs/heads/{branch}:{ref}', gateway=project == SOURCE)
        sha = self.git('rev-parse', ref).stdout.strip()
        if not SHA.fullmatch(sha):
            raise ValueError('Invalid commit SHA')
        return sha

    def ancestor(self, older, newer):
        result = self.git('merge-base', '--is-ancestor', older, newer, check=False)
        if result.returncode not in (0, 1):
            raise RuntimeError('Ancestry check failed')
        return result.returncode == 0

    def merged_commit(self, left, right):
        result = self.git('merge-tree', '--write-tree', left, right, check=False)
        if result.returncode == 1:
            return None
        if result.returncode:
            raise RuntimeError('Merge-tree failed')
        tree = result.stdout.splitlines()[0]
        if not SHA.fullmatch(tree):
            raise RuntimeError('Invalid merge tree')
        return self.git('commit-tree', tree, '-p', left, '-p', right,
            stdin=f'Merge gateway update into private OpenQP\n\nParents: {left} {right}\n').stdout.strip()

    def push_branch(self, commit, branch):
        if not branch.startswith(PREFIX):
            raise ValueError('Refusing to push a non-service branch')
        self.git('push', f'{BASE}/{PATHS[TARGET]}.git', f'{commit}:refs/heads/{branch}')
        remote = self.git('ls-remote', f'{BASE}/{PATHS[TARGET]}.git', f'refs/heads/{branch}').stdout.split()
        if not remote or remote[0] != commit:
            raise RuntimeError('Published SHA verification failed')

    def rotate_token(self):
        for project in (SOURCE, TARGET):
            self.rotate_project_token(project)

    def rotate_project_token(self, project):
        token_file = self.state / ('gateway.token' if project == SOURCE else 'gitlab.token')
        details = self.api('GET', '/personal_access_tokens/self', token_project=project)
        expires = details.get('expires_at')
        if not expires:
            return
        remaining = (dt.date.fromisoformat(expires) - dt.date.today()).days
        if remaining > 30 or self.dry_run:
            return
        rotated = self.api('POST', f'/projects/{project}/access_tokens/self/rotate',
                           {'expires_at': (dt.date.today() + dt.timedelta(days=365)).isoformat()},
                           token_project=project)
        token = rotated['token']
        tmp = token_file.with_suffix('.token.new')
        fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
        with os.fdopen(fd, 'w') as stream:
            stream.write(token + '\n')
            stream.flush()
            os.fsync(stream.fileno())
        tmp.replace(token_file)

    def upstream_gate(self, gateway):
        pipelines = self.pages(f'/projects/{SOURCE}/pipelines?sha={gateway}', token_project=SOURCE)
        if not pipelines:
            return 'waiting_upstream_ci'
        pipeline = max(pipelines, key=lambda p: p['id'])
        details = self.api('GET', f'/projects/{SOURCE}/pipelines/{pipeline["id"]}', token_project=SOURCE)
        jobs = self.pages(f'/projects/{SOURCE}/pipelines/{pipeline["id"]}/jobs', token_project=SOURCE)
        verdict = pipeline_gate({'sha': gateway, 'head_pipeline': details}, jobs, SOURCE)
        return 'clean' if verdict == 'clean' else 'waiting_upstream_ci'

    def pure_sync(self, mr, source, gateway, target):
        pinned = sync_gateway(TARGET, mr)
        if not pinned or not self.ancestor(pinned, gateway) or not self.ancestor(pinned, source):
            return None
        expected = self.git('merge-tree', '--write-tree', target, pinned, check=False)
        if expected.returncode:
            return None
        tree = expected.stdout.splitlines()[0]
        if self.git('rev-parse', source + '^{tree}').stdout.strip() == tree:
            return pinned
        # The initial approved MR also contains one documented operations file.
        # Pin its full source SHA and blob; never exempt an arbitrary docs subtree.
        if source == INITIAL_SHA:
            changes = self.git('diff', '--name-only', tree, source).stdout.splitlines()
            blob = self.git('rev-parse', source + ':' + INITIAL_DOC).stdout.strip()
            if changes == [INITIAL_DOC] and blob == INITIAL_DOC_BLOB:
                return pinned
        return None

    def handle_mr(self, mr, gateway, target):
        iid = mr['iid']
        prefix = f'/projects/{TARGET}/merge_requests/{iid}'
        mr = self.api('GET', prefix)
        if (mr['source_project_id'] != TARGET or mr['target_project_id'] != TARGET
                or mr['target_branch'] != 'main' or mr['state'] != 'opened'):
            return self.status('blocked_mr_scope', mr=iid)
        source = self.fetch(TARGET, mr['source_branch'], f'refs/sync/mr-{iid}')
        if source != mr['sha']:
            return self.status('waiting_source_changed', mr=iid)
        # Never mutate an adopted human branch, or merge an untested new target.
        if not self.ancestor(target, source):
            if self.dry_run or not mr['source_branch'].startswith(PREFIX):
                return self.status('blocked_target_changed', mr=iid)
            merged = self.merged_commit(source, target)
            if merged is None:
                return self.status('blocked_target_conflict', mr=iid)
            self.push_branch(merged, mr['source_branch'])
            return self.status('refreshed_target', mr=iid, sha=merged)
        pinned = self.pure_sync(mr, source, gateway, target)
        if pinned is None:
            return self.status('blocked_non_sync_changes', mr=iid, sha=source)
        verdict = self.upstream_gate(pinned)
        if verdict != 'clean':
            return self.status(verdict, mr=iid, gateway=pinned)
        pipeline = mr.get('head_pipeline') or {}
        jobs = (self.pages(f'/projects/{TARGET}/pipelines/{pipeline["id"]}/jobs')
                if pipeline.get('id') else [])
        verdict = pipeline_gate(mr, jobs)
        if verdict != 'clean':
            return self.status(verdict, mr=iid, sha=source)
        discussions = self.pages(prefix + '/discussions')
        if any(n.get('resolvable') and not n.get('resolved')
               for d in discussions for n in d.get('notes', [])):
            return self.status('blocked_discussions', mr=iid)
        if mr.get('draft') or mr.get('detailed_merge_status') != 'mergeable':
            return self.status('waiting_mergeable', mr=iid)
        live_target = self.api('GET', f'/projects/{TARGET}/repository/branches/main')['commit']['id']
        fresh = self.api('GET', prefix)
        if live_target != target or fresh['sha'] != source:
            return self.status('waiting_refs_changed', mr=iid)
        if self.dry_run:
            return self.status('would_merge', mr=iid, sha=source)
        # GitLab FF-only policy atomically rejects a newly diverging target.
        policy = self.api('GET', f'/projects/{TARGET}')
        if policy.get('merge_method') != 'ff' or policy.get('squash_option') == 'always':
            return self.status('blocked_project_policy', mr=iid)
        self.api('PUT', prefix + '/merge', {'sha': source, 'squash': False, 'should_remove_source_branch': False})
        fresh = self.api('GET', prefix)
        if fresh['state'] != 'merged':
            return self.status('waiting_merge_completion', mr=iid)
        main = self.fetch(TARGET, 'main', 'refs/sync/internal')
        if not self.ancestor(source, main):
            raise RuntimeError('Merged source missing from internal main')
        return self.status('merged', mr=iid, sha=source, main=main,
                           latest_gateway_included=self.ancestor(gateway, main))

    def run(self):
        self.rotate_token()
        project = self.api('GET', f'/projects/{TARGET}')
        if (project['path_with_namespace'] != PATHS[TARGET]
                or project.get('merge_method') != 'ff'
                or project.get('squash_option') == 'always'
                or not project.get('only_allow_merge_if_pipeline_succeeds')
                or not project.get('only_allow_merge_if_all_discussions_are_resolved')):
            return self.status('blocked_project_policy')
        if not self.repo.exists():
            self.git('init', '--bare')
        gateway = self.fetch(SOURCE, 'main', 'refs/sync/gateway')
        target = self.fetch(TARGET, 'main', 'refs/sync/internal')
        if self.ancestor(gateway, target):
            return self.status('up_to_date', gateway=gateway, internal=target)
        bot_id = self.api('GET', '/user')['id']
        mrs = self.pages(f'/projects/{TARGET}/merge_requests?state=opened&target_branch=main')
        owned = [m for m in mrs if m['source_project_id'] == TARGET and
                 ((m['source_branch'].startswith(PREFIX) and m['author']['id'] == bot_id)
                  or (m['iid'] == ADOPT_IID and m['source_branch'] == ADOPT_BRANCH))]
        if len(owned) > 1:
            return self.status('blocked_multiple_sync_mrs', mrs=[m['iid'] for m in owned])
        if owned:
            return self.handle_mr(owned[0], gateway, target)
        # A manually closed initial MR must not be silently recreated.
        initial = self.api('GET', f'/projects/{TARGET}/merge_requests/{ADOPT_IID}')
        if initial['state'] == 'closed' and initial['source_branch'] == ADOPT_BRANCH:
            return self.status('blocked_closed_initial_mr', mr=ADOPT_IID)
        branch = f'{PREFIX}{gateway}/{target}'
        previous = self.pages(f'/projects/{TARGET}/merge_requests?state=all&source_branch={urllib.parse.quote(branch, safe="")}')
        if previous:
            return self.status('blocked_existing_sync_record', mrs=[m['iid'] for m in previous])
        result = self.git('ls-remote', f'{BASE}/{PATHS[TARGET]}.git', f'refs/heads/{branch}').stdout.split()
        if result:
            candidate = self.fetch(TARGET, branch, 'refs/sync/recovery')
            parents = self.git('show', '-s', '--format=%P', candidate).stdout.strip().split()
            expected_tree = self.git('merge-tree', '--write-tree', target, gateway).stdout.splitlines()[0]
            if parents != [target, gateway] or self.git('rev-parse', f'{candidate}^{{tree}}').stdout.strip() != expected_tree:
                return self.status('blocked_unexpected_existing_branch', branch=branch)
        else:
            candidate = self.merged_commit(target, gateway)
            if candidate is None:
                return self.status('blocked_merge_conflict', gateway=gateway, internal=target)
            if self.dry_run:
                return self.status('would_create_mr', gateway=gateway, internal=target)
            self.push_branch(candidate, branch)
        if self.dry_run:
            return self.status('would_recover_mr', sha=candidate)
        mr = self.api('POST', f'/projects/{TARGET}/merge_requests', {
            'source_branch': branch, 'target_branch': 'main', 'target_project_id': TARGET,
            'title': f'Sync gateway main {gateway[:12]} into private OpenQP',
            'description': f'Inbound-only integration of gateway `{gateway}` into internal `{target}`. '
                'Upstream verification is reused for unchanged inbound code; full internal integration CI gates merging. '
                f'Source: {BASE}/{PATHS[SOURCE]}/-/commit/{gateway}',
            'remove_source_branch': False, 'squash': False})
        if mr['source_project_id'] != TARGET or mr['target_project_id'] != TARGET:
            raise RuntimeError('Unexpected MR destination')
        return self.status('created_mr', mr=mr['iid'], sha=candidate, gateway=gateway)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--state', type=Path, default=Path('/var/lib/openqp-inbound-sync'))
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    args.state.mkdir(parents=True, exist_ok=True)
    with (args.state / 'service.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return 0
        sync = Sync(args.state, args.dry_run)
        try:
            sync.run()
        except urllib.error.HTTPError as exc:
            sync.status('error', kind='http', code=exc.code)
            return 1
        except Exception as exc:
            sync.status('error', kind=type(exc).__name__)
            return 1
    return 0

if __name__ == '__main__':
    sys.exit(main())
