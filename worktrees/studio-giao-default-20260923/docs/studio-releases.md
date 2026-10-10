# Studio releases and the bundled OpenQP engine

## Engine upstream (implemented)

Studio develops in the private GitLab `open-quantum-platform/oqp-studio`
project. Its distributable engine comes exclusively from
`https://qchemlab.knu.ac.kr/open-quantum-platform/openqp`, the public export
gateway linked to the official GitHub OpenQP repository. The gateway itself
is access-controlled; "public" describes the content selected for external
publication. `internal/openqp` is not an engine source for Studio distribution.

`engine-upstream.json` declares the gateway and tracking branch. Run:

```sh
python3 tools/checkout_engine.py --ssh
```

On Ultra this uses the verified `gitlab-qchemlab` SSH route, checks out the
current gateway `main` into a new `openqp/` directory, and writes
`engine-source.json` with the full commit. Existing destinations are refused.
To reproduce a candidate, pass `--ref <full-commit-sha>` and a new destination.
A ref outside gateway main's ancestry is refused. There is no automatic
fallback to internal OpenQP or GitHub when gateway access fails.

CI uses HTTPS with ephemeral `CI_JOB_TOKEN` authentication. The gateway must
allow inbound job tokens from Studio, and the initiating user must have
permission to read the gateway. Credentials never enter the checkout URL or
manifest. A legacy external runner needs a separately provisioned read-only
`OPENQP_GATEWAY_TOKEN`; do not re-enable GitHub Actions just to build Studio.
The historical GitHub release workflow now resolves main once and supplies
that exact SHA to every engine platform job. This also removes the old
hard-coded engine SHA and the latest-GitHub-release default.

Resolving latest main is candidate selection, **not** release qualification.
Gateway CI success alone does not verify Studio compatibility. Engine/native
builds, installer tests, and the feature checks below remain required.

## Findings on 2026-09-22

- Studio GitLab/GitHub main is `2cb4fe93eeabcd3c04d1c48a354df305757e0a4b`.
  GitHub source is private and its Actions are disabled.
- OpenQP GitHub release `v1.3.1` carries Studio `0.2.3` executables as well as
  the independent OpenQP engine archives. Older downloads must remain valid.
- The historical Studio release build pins engine commit
  `49583799238d1ff02387463b33683dbe853d4d80`; standalone workflows instead
  default to GitHub's latest engine release.
- The app's update API and Releases menu point to the private Studio GitHub
  repository. Ordinary users cannot use that endpoint. Replacing it with
  OpenQP's `/releases/latest` would also be wrong: that tag names the engine,
  not Studio, and a later engine release need not include Studio installers.
- Backend and Tauri versions are 0.2.3; frontend package metadata remains
  0.1.0. Engine runtime labels use package version, without the commit.
- The legacy publisher includes wheel/sdist and optional PyPI publishing.
  That predates the decision to keep Studio source private. It must not be
  used as the new public publisher.

## Selected distribution: public Studio Pages (activation pending)

The user selected this distribution architecture on 2026-09-23. Keep the existing
GitHub Pages documentation site and use its top-level **OQP Studio** menu
with **Download**, **Release notes**, **Engine compatibility**, and **Manual**.
Use `https://open-quantum-platform.github.io/openqp-docs/studio/` as the landing
page and `/studio/download/` as the stable download URL. This avoids another
hosting service/domain while giving Studio a separate product entry point.
The existing `/studio/packages/` page can link to the new download page.

Create a separate **public binary-distribution repository**, named
`Open-Quantum-Platform/oqp-studio-releases`. It contains release notes,
checksums, manifest and installers only. Keep Studio development/source in
private GitLab. Do not mirror Studio source, Python wheels/sdists, source maps,
private branch history, or build credentials into this distribution repo.
GitHub's automatically generated source archives then contain only the small
public distribution repository, not Studio source.

Use Studio tags such as `v0.2.4`, independently of OpenQP tags. OpenQP's release
page should link to Studio downloads; new Studio installers should stop being
attached to engine releases. Retain the old v1.3.1 assets and URLs for existing
users. Public publishing remains an explicit on-demand promotion from GitLab.

### Two versions, plus an exact engine commit

Every release manifest and About dialog should identify:

- Studio version and source commit;
- bundled engine package version **and full gateway commit**, plus upstream
  release tag only if the commit is exactly that tag;
- channel, platform/architecture, asset SHA-256 and required feature set;
  retain private build-pipeline provenance in GitLab;
- the actually selected runtime engine, separately from the bundled one, when
  a user chooses a local installation or remote execution.

For an engine beyond its latest release, display for example
`OpenQP 1.3.1 + gateway abcdef012345 (snapshot)` rather than claiming it is the
unaltered 1.3.1 release. This is a display identity; it must not rewrite
upstream package metadata. Studio and engine version numbers need not match.

An engine-only bundle refresh is still a new Studio patch/preview release;
never replace an already released installer under the same version.

### Stable and Preview

**Stable** is a tested Studio/engine pair. Its engine may be a qualified gateway
snapshot; it does not have to wait for the next formal OpenQP release.
**Preview** follows newer gateway main commits and is explicitly opt-in.

On a gateway change, prepare a candidate with one resolved engine SHA. Build
all platforms from that SHA. Test representative Studio workflows, input
keywords, result formats, native execution and offline all-in-one startup.
New features must declare the required engine capability/commit and stay
unavailable with an older user-selected engine. Promote only the exact tested
artifacts, without rebuilding between qualification and publication.
Automating candidate creation is a later CI step, not an automatic public
release or silent replacement of an installed engine.

### Migration order

1. Adopt the gateway-only checkout and source tests in this MR.
2. Port the existing platform build jobs to GitLab runners. Preserve the
   existing native build/cache/ILP64 rules. Add compatibility and binary-only
   payload gates; unify Studio version metadata.
3. Add the distribution repo, release manifest, and docs download menu.
   Before public release, explicitly synchronize the selected public-safe engine
   commit to the official GitHub gateway mirror and verify its source link.
   Publish a separately authorized pilot candidate as a draft, check anonymous
   download access and checksums, then publish it.
4. Switch both frontend Releases and backend update/engine-download endpoints
   to the public Studio feed in the same release. Keep Stable/Preview feeds
   distinct and bind engine downloads to the selected Studio manifest. Verify
   installer variant preservation (with-engine vs slim), including Linux.
5. Test upgrade from the currently distributed 0.2.3 installers. If their
   embedded updater still uses the now-private URL, the first upgrade needs
   the new public installer/download page; changing the server alone cannot
   repair a hard-coded URL in an already installed app.
6. After one successful migration release, add Studio links to engine release
   notes. Leave historical assets in place.

### CI responsibilities

The public Pages site is the product entry point. It does not execute Studio
or build native installers. Split the work into three explicit operations:

1. **MR validation in private GitLab.** Run backend lint/tests and frontend
   install/audit/build on native Windows, macOS and Linux runners for the exact
   MR head. This is the six-job matrix carried over from GitHub CI. These
   checks do not qualify a downloadable installer by themselves.
2. **Candidate packaging in private GitLab.** On an explicit candidate request,
   resolve one gateway engine SHA, build the standard and with-engine desktop
   installers, smoke-test installation/startup and required NMR/ACID workflows,
   and retain the artifacts plus checksums. Reuse the platform build scripts;
   do not use the old GitHub source-repository release publisher or attach new
   installers to an OpenQP engine release. Qualify every advertised architecture,
   including both Mac architectures if both are offered.
3. **Public promotion and Pages update.** Promote those exact tested artifacts,
   without rebuilding, to `Open-Quantum-Platform/oqp-studio-releases`. Publish
   release notes, `SHA256SUMS` and a manifest alongside them. Then update the
   public docs download page with the immutable release/asset URLs. The existing
   docs Pages workflow publishes the website; no Studio build credential is
   available to that workflow.

The first Pages revision keeps verified Studio 0.2.3 download links as an
explicitly historical release. It must not present development features as
released or link to a distribution repository that has not been created yet.
Switch app updates and engine downloads only together with the first qualified
independent release and its public feed. Do not repoint an installed app at an
empty endpoint.

The candidate manifest must contain the Studio version and source SHA, engine
version and full gateway SHA, channel, each asset's platform/architecture and
variant, immutable download URL, size and SHA-256. The public manifest must not
contain private CI URLs, repository credentials or source archives. Keep the
private CI provenance with the candidate artifacts instead.

### Release notes and version comparison

The public `studio/whats-new/` page is the release-note archive. Each entry
must name the previous Studio release and contain a before/after comparison,
features, fixes, known limitations and upgrade instructions, together with the
Studio version, date, channel and exact bundled-engine identity. Link that
entry from the download page and public release description.

Keep three states distinct: published releases, changes implemented in the
development source, and work still under review. Promote only the notes for
changes contained in the exact qualified installers. Never infer missing
historical release notes or assign a future version/date before qualification.

### Current activation state

The Pages content and this distribution decision are prepared on isolated
branches. The binary-only repository, updater migration and native installer
promotion are not activated by this documentation change. The existing NMR
merge prerequisite and full platform qualification remain release gates.
Historical public installers are retained unchanged.

## References

- [GitLab job-token access and private repository cloning](https://docs.gitlab.com/ci/jobs/ci_job_token/)
- [GitHub release assets and automatically generated archives](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases)


## Studio 0.2.4 candidate implementation

`release-candidate.json` pins the public engine commit and the Studio version.
`tools/release/freeze_engine.py` freezes an already-built isolated engine;
`qualify_engine.py` runs the bundled NMR/ACID example without development tools
on PATH and checks its energy, shielding record and four finite cube grids.
`package.py` builds both installer variants and checks the frozen sidecar.
`manifest.py` requires all four advertised native architectures, both variants,
matching exact source commits, successful qualification and unchanged hashes.
Public promotion copies these artifacts without rebuilding them.

The app reads the public binary-distribution manifest, verifies each download's
size and SHA-256 before installation, and preserves the installed variant.
Linux offline updates select the with-engine DEB. An absent offline variant
blocks the update. Engine downloads use the installed Studio tag, independently
of whichever Studio version is newest. Stable is the default; Preview is an
explicit About-menu selection. The About panel identifies the release engine
and separately reports the selected local/bundled runtime when it is known.
ACID submissions require an engine manifest declaring the qualified capability;
older or unverified local/remote engines are not silently accepted for ACID.

The 0.2.3 updater has a hard-coded private endpoint. Its first migration requires
a manual install from the public download page. Historical assets stay intact.

Mac release packaging preserves PyInstaller BUNDLE framework/resource links,
signs the enclosing Tauri app with an ad-hoc identity, and verifies every Mach-O
minimum deployment target plus the complete signature before archiving it.
Ad-hoc signatures are not Apple notarization; the normal first-launch approval
in macOS Privacy & Security can still be required. No paid certificate or
notarization service is provisioned by this release process.
