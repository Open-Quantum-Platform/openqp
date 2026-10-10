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
4. Wait for the exact source SHA's GitLab CI and automated Codex review. Address
   actionable integration issues without relaxing tests or branch protections.
   Only merge after CI succeeds, review issues are resolved, and the source SHA
   still matches. Re-evaluate if target main has moved; integrate and retest it
   when required. Do not interpret pending, missing, skipped or failed CI as
   success. Enable GitLab auto-merge only after review is clear and required CI
   enforcement is active.
5. Verify remote internal main contains the gateway tip. Record the result in
   `/var/lib/openqp-inbound-sync/status.json` and the server journal. Deployment
   provenance is tracked separately in research operations task
   `CHC-20260922-608F70`; runtime synchronization makes no registry calls.

Server implementation and operation instructions are in internal MR !5, under
`devkit/tools/gateway-sync/`. The service adopts initial MR !4 without changing
its source branch; CI failures and review findings continue to block merging.

Initial integration: gateway `255aef4ac87cae7c18f70d987ebad9b052f76481` into
internal `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`. The eight incoming commits
cover already published SCF/TRAH and QM/MM active-atom changes. Their integration
is conflict-free; existing private development content is retained.

This automation does not alter the conditional Studio release plan, global
skills, local working branches, or on-demand-only GitHub outbound mirrors.
