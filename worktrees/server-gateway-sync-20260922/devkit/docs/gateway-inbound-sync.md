# Inbound gateway updates into the private development repository

The user-authorized direction is:

`GitHub Open-Quantum-Platform/openqp:main` →
`GitLab open-quantum-platform/openqp:main` (project 17) →
`GitLab open-quantum-platform/internal/openqp:main` (project 19).

The first step is the existing GitHub inbound-sync workflow. The second step
runs through `openqp-inbound-sync.timer` on the existing review server
(`codex-review-runner`, 192.168.1.46), five minutes after each completed run.
Ultra and the Codex desktop app are not runtime dependencies. The old desktop
heartbeat is paused. It does not update filesystem development checkouts.
It is not an outbound mirror and never copies private history to the gateway.

For each incoming update:

1. Fetch the live gateway and internal main refs. If the gateway tip is already
   an ancestor of internal main, do nothing. Check existing sync MRs first.
2. Create a unique sync branch from internal main in the service bare repository. Merge
   the exact gateway tip into it without rebasing, resetting, or force-pushing.
   Preserve private history and devkit content. If there are conflicts, stop
   automatic integration and report the affected files; never choose a side
   globally or discard private changes.
3. Inspect the integration diff, run applicable source checks, publish only to
   project 19, and create a normal MR to its main branch through the API with
   `target_project_id=19` explicitly set. Do not use Git push options to create
   the MR: fork defaults can choose upstream project 17 as the target. Verify
   both source and target project IDs are 19. Include both source
   SHAs and the gateway commit link. Serialize sync MRs so newer gateway commits
   do not cause simultaneous competing merges. Never merge unrelated MRs.
4. Reuse upstream verification for unchanged incoming code. Verify that the pinned
   gateway commit is accepted on gateway main and its seven upstream CI jobs
   passed. Require the candidate tree to exactly equal the normal merge of
   current internal main and the pinned gateway commit; unexpected integration
   edits stop automatic merging and require ordinary development review. Initial
   MR !4 has an exact source-SHA/document-blob exception for this operations note.
   Require all seven internal CI jobs to pass for the exact candidate SHA, resolved
   discussions, matching refs, and GitLab mergeability. Fast-forward-only merging
   rejects a concurrently diverging target, and squashing is disabled. Pure sync
   MRs do not require or trigger duplicate Codex review. Ordinary development MRs
   retain full automatic review. Existing upstream findings are tracked separately.
5. Verify remote internal main contains the gateway tip. Record the result in
   `/var/lib/openqp-inbound-sync/status.json` and the server journal. Deployment
   provenance is tracked separately in research operations task
   `CHC-20260922-608F70`; runtime synchronization makes no registry calls.

Server implementation and operation instructions are in internal MR !5, under
`devkit/tools/gateway-sync/`. The service adopts initial MR !4 without changing
its source branch; failed CI, conflicts and unreviewed integration edits block merging.

Initial integration: gateway `255aef4ac87cae7c18f70d987ebad9b052f76481` into
internal `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`. The eight incoming commits
cover already published SCF/TRAH and QM/MM active-atom changes. Their integration
is conflict-free; existing private development content is retained.

This automation does not alter the conditional Studio release plan, global
skills, local working branches, or on-demand-only GitHub outbound mirrors.

MR !4 was merged by the server worker on 2026-09-22 after upstream pipeline
2680 and internal pipeline 2702 passed. Internal main contains gateway
`255aef4ac87cae7c18f70d987ebad9b052f76481`; private history is retained.
