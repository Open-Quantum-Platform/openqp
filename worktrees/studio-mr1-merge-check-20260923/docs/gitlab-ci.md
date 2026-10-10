# Studio continuous integration

GitLab merge requests and `main` run the same platform matrix as
`.github/workflows/ci.yml`: native Linux, macOS and Windows jobs for both
backend and frontend. A pipeline must pass before merging.

| Check | Runtime | Commands |
| --- | --- | --- |
| Backend, each OS | Python 3.12 | Install `backend[dev]` in a job-local venv, `ruff check .`, `pytest -v` |
| Frontend, each OS | Node 24 | `npm ci`, responsiveness and ACID-control regressions, `npm audit --audit-level=moderate`, `npm run build` |
| Engine source | Linux / Python 3.12 | Engine checkout unit tests and public-gateway checkout with recorded SHA |

Linux uses the `small-linux` Docker runner with official Python/Node images.
macOS uses the existing `macmaster,macos,shell` runner; Windows uses the existing
`windows,winserver1,shell` runner. The scripts check the actual operating system:
a Linux container is not accepted as a Windows or macOS test.

Shell runners keep their existing OpenQP Python installations. When necessary,
`tools/ci/backend.py` obtains Python 3.12 through a job-local uv installation.
The macOS frontend bootstrap downloads Node 24.19.0 and verifies the published
SHA-256 checksum before extraction. All downloaded tools and package caches
are confined to the job checkout. Windows must provide Node 24 and a Python
interpreter capable of creating a venv. A missing prerequisite fails the job.

Backend results are retained as JUnit reports; each successful frontend build
is retained as an artifact for one week. Jobs fail on lint, test, dependency
audit or build failure, and required platform checks are not manual or allowed
to fail. A runner being unavailable leaves its job pending and blocks merging.

The installer and standalone-engine packaging workflows in `.github/workflows/`
are separate release operations, not part of the GitHub PR CI matrix. This CI
change does not publish a release or run those packaging workflows.
