# Desktop shell (Tauri 2)

Tauri 2 shell wrapping the frontend and the frozen backend. It builds the
published installers: MSI/NSIS on Windows, DMG and `.app` tarball on macOS,
AppImage and `.deb` on Linux — each in a plain and a `-with-engine` variant.

## Layout

```
src-tauri/
  Cargo.toml         Rust crate (tauri 2 + shell plugin)
  tauri.conf.json    window, bundled frontend dist, bundling
  src/main.rs        app entry
  capabilities/      window permission set
  binaries/          where the frozen backend is staged (empty in git)
  icons/             placeholder icons (replace with real OQP artwork)
```

## Developing

Requires Rust (stable) and, on Linux, the Tauri 2 prerequisites
(webkit2gtk-4.1, libappindicator — see the tauri.app docs).

```bash
cargo install tauri-cli --version "^2"
cd shell/src-tauri
cargo run            # opens the bundled frontend in a native window
cargo tauri build    # bundles installers for the host OS
```

`cargo run` needs a frozen backend staged in `binaries/` first, because
`externalBin` declares one; the repository ships that directory empty. See
"Working on the desktop shell" in the top-level README for the one-time build,
and for why most work does not belong in this loop at all.

The desktop window serves `frontend/dist` from Tauri's internal origin. It
does not use Vite or a loopback HTTP origin.

## What the shell does

1. Spawns the PyInstaller-frozen backend as a sidecar with `--stdio`, passing
   `OQP_STUDIO_RESOURCES` so it can find an engine bundled by the installer —
   the resource directory sits next to the executable on Windows, in
   `Contents/Resources` on macOS, and under `/usr/lib` on Linux, and only
   Tauri knows which.
2. Routes frontend `/api` requests through a Tauri command as newline-delimited
   JSON. No TCP port is allocated.
3. Keeps the backend child process until the desktop window exits.

## Remaining

- Real icons (`icon.ico`, `icon.icns`) and file associations for `.oqp`,
  `.molden`, `.cube`, `.hess.json`.
- Code signing and notarization; see [../docs/code-signing.md](../docs/code-signing.md).
- Forward sidecar stderr to `backend.log` so startup remains observable in the
  field.
