# OQP Studio logo concept 01

Status: original approved concept, superseded by the user-selected pink/blue
palette in `studio-logo-pink-blue.md`. Retained as design history.
Installer publication remains deferred until the relevant NMR PR/MR merges.
Generated using the built-in image generation tool on 2026-09-22.
Image: `studio-logo-concept-01.png`.

The abstract orbital-inspired Q uses two contrasting colors; it is a brand
symbol rather than a literal orbital surface. Desktop icon variants were exported without redesign using Tauri CLI 2.11.5.
The desktop shell uses PNG, Windows ICO and macOS ICNS; the Studio web
development frontend uses the 32px PNG favicon. The embedded OpenqpView viewer
retains its separate identity. No installed application bundle was modified.

## Exact generation prompt

Use case: logo-brand. Create one polished, original app icon concept for OQP Studio, a professional molecular quantum-chemistry visualization desktop application. A single square 1024x1024 image, icon artwork only, no presentation board, no wordmark or text. Design a distinctive sculptural Q emblem abstracted from smooth molecular orbital lobes: two broad opposing curved lobes enclose a precise negative-space aperture, with a short elegant diagonal lower-right stroke making the Q unmistakable. The lobes evoke the two phases of an orbital without pretending to be a literal scientific plot. Strong simple silhouette, expertly balanced geometry, generous margins, excellent legibility at small app-icon sizes. Rich restrained color, luminous cool blue-teal opposed to warm copper-orange, subtle satin depth rather than glossy plastic. Set on a beautifully refined deep charcoal rounded-square app tile with just a little tonal depth; straight-on, perfectly centered, no perspective mockup. Make the mark sophisticated, memorable and suitable for a serious scientific application, with crisp vector-like contours. Avoid generic atom symbols with orbiting electrons, balls-and-sticks, tiny particles, busy details, neon cyberpunk, lens flares, excessive glow, gradients that muddy the edges, stock infinity symbols, text, labels, watermarks and any existing company's logo. This is a fresh concept, not an edit of the existing icon.

## Reproduce the platform exports

```sh
npm exec --yes --cache .cache/npm --package @tauri-apps/cli@2.11.5 -- \
  tauri icon docs/design/studio-logo-concept-01.png --output .cache/generated-icons
```

Copy `32x32.png`, `128x128.png`, `128x128@2x.png`, `icon.png`, `icon.ico`,
and `icon.icns` to `shell/src-tauri/icons/`; use `32x32.png` for
`frontend/public/studio-icon.png`. Only desktop outputs are checked in.
