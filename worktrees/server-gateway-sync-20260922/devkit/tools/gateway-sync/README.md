# OpenQP inbound sync service

Runs on the existing `codex-review-runner` server (192.168.1.46), independently
of Ultra, Codex desktop, or a user login. The systemd timer checks five minutes
after each completed run. The previous Codex heartbeat is paused.

Direction: official GitHub main → existing GitLab gateway inbound workflow →
private `open-quantum-platform/internal/openqp` (project 19).
This service implements only the second arrow and never writes to project 17
or GitHub. It does not update filesystem development checkouts.

Separate dedicated bots provide project-19 Maintainer access and project-17
Reporter access. Their tokens live only in the private server state directory.
Tokens are rotated through the self-rotation API when 30 days remain; rotation
failure stops the run. No credentials belong in this repository.

The worker fetches into a bare repository, combines histories using merge-tree
and commit-tree without checking out or executing source, and opens one internal
MR at a time. Private files/history are retained; conflicts stop integration.
An interrupted push can recover its MR only after parent/tree verification.
It adopts initial internal MR !4 without modifying that human-created branch.
A changed target on an adopted branch requires manual integration and retesting.

Merging reuses upstream review for unchanged inbound code. The pinned gateway
commit must be accepted on gateway main with all seven upstream CI jobs passing.
The candidate must exactly match a conflict-free merge of target and gateway;
the initial MR has a fixed SHA/blob exception for its single operations document.
New integration edits stop automatic merging and require ordinary development
review. All seven internal CI jobs must pass on the exact candidate SHA, with resolved discussions,
matching live refs, and GitLab mergeability. Missing upstream verification, new integration code, missing jobs, failed tests,
conflicts, or changed refs stop merging. Final merging uses
the GitLab MR API with the expected source SHA and the project's native gates.
Internal project 19 must use fast-forward-only merging. This makes GitLab reject
a target that diverges after the last check instead of merging an untested tree.
Contributors can merge current main into their branch and rerun CI; history
rewriting is not required. Squashing is explicitly disabled to preserve ancestry. Unknown new CI formats require an operator to review this implementation.

## Server operations

Installed code: `/usr/local/libexec/openqp-inbound-sync.py` and
`/usr/local/libexec/openqp-inbound-askpass` (root-owned).
State: `/var/lib/openqp-inbound-sync` (0700, `openqp-sync` user).
Latest result: `status.json`; execution log: systemd journal.
A process lock excludes overlapping runs. No Ultra registry calls occur at runtime.

```sh
sudo systemctl status openqp-inbound-sync.timer
sudo journalctl -u openqp-inbound-sync.service -n 20 --no-pager
sudo cat /var/lib/openqp-inbound-sync/status.json
sudo -u openqp-sync /usr/bin/python3 /usr/local/libexec/openqp-inbound-sync.py --dry-run
sudo systemctl start openqp-inbound-sync.service
```

Pause/rollback: `sudo systemctl disable --now openqp-inbound-sync.timer`.
This preserves all Git refs and MRs. Never solve a sync failure by force-pushing
or enabling an outbound private mirror. Do not enable the old desktop heartbeat
alongside this timer. Deployment is tracked in CHC-20260922-608F70.

Validation: `python3 -m unittest discover -s devkit/tools/gateway-sync -p 'test_*.py' -v`.

The review poller skips only project-19 MRs from the dedicated sync bot (64)
with fully pinned branch names, plus the exact approved initial MR !4 SHA.
Other OpenQP/Studio development MRs retain automatic Codex review. The shared
root-owned `openqp_sync_policy.py` defines this narrow identity check. Existing
upstream findings remain recorded in MR !4 and are tracked separately; they
are not erased or silently marked resolved by the sync worker.
