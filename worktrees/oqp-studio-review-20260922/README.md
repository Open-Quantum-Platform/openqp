# OpenQP Codex review runner

Dedicated review worker for these GitLab projects:

- 17: `open-quantum-platform/openqp`
- 19: `open-quantum-platform/internal/openqp`
- 18: `open-quantum-platform/oqp-studio`

Each open merge request is reviewed once per source commit. New source commits
are picked up by the existing timer. Review notes are posted by
`codex-review-bot`; this service does not merge merge requests.

The bot inherits Developer access from the private group. For each project,
the `codex-review` user's Git configuration must trust the exact root-managed
mirror path `/var/lib/gitlab-review/mirror/project-<id>.git` so Git can clone
the mirror across the two service identities. Do not use a wildcard trust path.

- VM: `codex-review-runner.qc.lab` (VMID 108)
- Poll interval: one minute, with up to ten seconds randomized delay
- GitLab identity: `codex-review-bot`
- Codex CLI: 0.155.1, authenticated with ChatGPT under `codex-review`
- GitLab credential: `/var/lib/gitlab-review/gitlab.token` (never copy into this source directory)
- Codex credential: `/var/lib/codex-review/.codex/auth.json` (never copy into this source directory)

The root-owned poller splits privileges between `gitlab-review` for GitLab fetch/API access and `codex-review` for Codex execution. Codex runs with its read-only sandbox. The Ubuntu AppArmor profile grants user-namespace creation only to the pinned Codex bwrap executable.

Active files:

- `/usr/local/sbin/codex-review-poller`
- `/usr/local/libexec/codex-review-askpass`
- `/etc/systemd/system/codex-review-poller.service`
- `/etc/systemd/system/codex-review-poller.timer`
- `/etc/apparmor.d/codex-review-bwrap`

After changing a source file, install it to the corresponding active path, run `systemctl daemon-reload` for unit changes, test with `systemctl start codex-review-poller.service`, and verify the GitLab note author and commit marker. Keep MR source content untrusted and do not expose either credential to the other OS user.

If Codex CLI is upgraded, update the exact bwrap path in the AppArmor profile and verify a sandboxed review before enabling the timer.
