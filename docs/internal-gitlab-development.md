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

The private CI acceptance matrix is deliberately platform-native:

- Linux x86_64 builds the package with OpenBLAS ILP64, runs the complete example
  and Python regression suites, and separately verifies a two-rank MPI build.
- macOS x86_64 on Zeus and macOS arm64 on macmaster build with Python 3.11,
  GCC 15, and the default Apple Accelerate ILP64 policy, then run the complete
  example suite. Each host has a runner-only persistent external cache.
- Windows x86_64 on Winserver1 builds with Intel oneAPI ifx/icx and MKL ILP64,
  runs the complete example suite, and also runs a fast portable source gate.

Runner authentication is restricted to the private OpenQP GitLab group. CI
jobs never receive GitHub write credentials. Platform caches are isolated by
runner and toolchain; jobs on a given runner execute serially so two builds do
not write the same external cache concurrently.
