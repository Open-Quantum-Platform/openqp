# GitLab access for cheol/openqp

- Use the `openqp-gitlab-development` skill for all OpenQP repository work.
- On this Mac, use `glab --hostname qchemlab.knu.ac.kr` to read and update
  merge-request discussions in `cheol/openqp`. Its Developer-level project
  token is stored in the macOS Keychain and expires on 2026-12-20; never ask
  the user to paste it, print it, or copy it into a repository or settings file.
- A mode-0600 compatibility copy for tools that require a token file is at
  `~/.gitlab-qchemlab-token`; read it only for GitLab authentication and never
  display its contents.
- Use the configured GitLab SSH remote for Git fetch and push. The API token is
  limited to `cheol/openqp`; it is not an upstream or instance-wide credential.
