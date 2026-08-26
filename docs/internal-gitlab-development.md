# Internal GitLab development and release boundary

The private CHC GitLab project is the system of record for OpenQP development.
Merge requests, review comments, issue discussion, CI traces, intermediate
branches, and unreleased source changes remain in GitLab.

GitHub is a publication target, not a development remote. Normal development
must not push branches, comments, pull requests, or CI state to GitHub. A GitHub
write is allowed only for an explicitly approved release after the release
commit, intended ref, artifacts, checksums, licensing gates, and live GitHub ref
have been verified. Plain force push is prohibited. Any exceptional
force-with-lease operation requires separate approval for the exact ref and
expected remote SHA.

Imports from `Open-Quantum-Platform/openqp` are read-only. Synchronization adds
new refs or fast-forwards existing refs and must not delete or overwrite
divergent internal branches. GitLab CI must not commit or push source code.
