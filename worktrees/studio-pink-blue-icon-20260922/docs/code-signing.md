# Code signing and notarization

Unsigned installers still work, but each OS shows a warning the first time:

- **macOS** refuses to open the app ("damaged" / "unidentified developer")
  until the quarantine flag is cleared with
  `xattr -cr "/Applications/OQP Studio.app"`, or the user approves it under
  System Settings → Privacy & Security.
- **Windows** SmartScreen shows "Windows protected your PC"; the user must
  click *More info → Run anyway*.
- **Linux** has no such gate; deb and AppImage need no signing.

The release workflow already supports signing. It is entirely driven by
repository secrets: with none configured the build produces unsigned
installers exactly as it does today, and adding them switches signing on with
no further code change.

## Does it have to be a paid certificate?

**macOS: yes for signing — but the fee can be waived.** Gatekeeper only trusts
Developer ID certificates issued through the Apple Developer Program; a
self-signed or ad-hoc signature is not trusted and does not remove the
warning. However, Apple waives the USD 99 fee for
[accredited educational institutions](https://developer.apple.com/help/account/membership/fee-waivers/),
nonprofits, and government entities, and South Korea is an eligible region.
A university enrolls as an organization (D-U-N-S number required) and selects
the fee-waiver option during enrollment, so signing and notarization end up
costing nothing. The waiver excludes entities selling paid apps or digital
goods, which does not apply here.

**Windows: yes, with no waiver.** Authenticode trust requires a certificate
from a commercial CA; self-signed certificates are not trusted.

**Neither is required to avoid the warning entirely** — see the tarball route
below.

## Installing without any certificate

macOS applies the quarantine flag that triggers Gatekeeper only when a file
arrives through an app that opts into it — browsers, mail clients, Messages.
Files fetched with `curl` are not quarantined, so an unsigned app installed
from the terminal opens normally:

```bash
curl -L -o oqp-studio.tar.gz <tarball URL>
tar xzf oqp-studio.tar.gz -C /Applications
open "/Applications/OQP Studio.app"
```

The release workflow publishes `OQP-Studio-*.app.tar.gz` alongside the DMG for
exactly this. The same idea applies on Windows: files downloaded by PowerShell
or a package manager (Scoop, winget) carry no mark-of-the-web, so SmartScreen
stays quiet.

## macOS — Developer ID + notarization

Requires an **Apple Developer Program** membership (USD 99/year, or free with
the educational fee waiver above). Notarization is not available on a free
Apple ID.

1. In the Apple Developer portal create a **Developer ID Application**
   certificate, install it in Keychain Access, then export it as a `.p12`
   with a password.
2. Base64-encode it: `base64 -i certificate.p12 | pbcopy`.
3. Create an app-specific password for your Apple ID at appleid.apple.com
   (Sign-In and Security → App-Specific Passwords).
4. Add these repository secrets (Settings → Secrets and variables → Actions):

   | Secret | Value |
   | --- | --- |
   | `APPLE_CERTIFICATE` | base64 of the `.p12` |
   | `APPLE_CERTIFICATE_PASSWORD` | the `.p12` export password |
   | `APPLE_SIGNING_IDENTITY` | e.g. `Developer ID Application: Your Name (TEAMID)` |
   | `APPLE_ID` | the Apple ID email |
   | `APPLE_APP_PASSWORD` | the app-specific password |
   | `APPLE_TEAM_ID` | the 10-character team ID |

Tauri then signs the bundle with a hardened runtime and submits the DMG for
notarization during `tauri build`; the resulting DMG opens with no warning and
no `xattr` step.

An App Store Connect API key (`APPLE_API_ISSUER` / `APPLE_API_KEY` /
`APPLE_API_KEY_PATH`) can replace the Apple ID pair; the workflow passes
whichever variables exist.

## Windows — Authenticode

Since June 2023 OV code-signing keys must live on hardware tokens or an HSM,
so a `.pfx` file in CI is no longer possible. Use a cloud signing service and
set one secret, `WINDOWS_SIGN_COMMAND`, to a command that signs the file
passed as `%1`:

- **Azure Trusted Signing** — about USD 10/month, the cheapest current option;
  individuals need a verifiable three-year identity history.
  Example: `trusted-signing-cli -e https://eus.codesigning.azure.net -a ACCOUNT -c PROFILE %1`
- **DigiCert KeyLocker**, **SSL.com eSigner** — roughly USD 200–600/year.

EV certificates additionally clear SmartScreen reputation immediately; OV
certificates build reputation over time.

Note that SmartScreen warnings fade as a signed binary accumulates downloads,
so signing is worth it mainly for wide distribution.

## Without signing

For a research group, the practical alternatives are:

- publish the `xattr -cr` one-liner next to the download (current approach);
- distribute through a package manager that handles the gate — a Homebrew
  cask installed with `--no-quarantine`, or conda-forge;
- have users build from source, where no gate applies.
