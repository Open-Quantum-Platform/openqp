# OQP Studio pink/blue icon

Status: final palette approved by the user and explicitly requested for application.
Source: `studio-logo-pink-blue.png`, generated with the built-in image tool from
the approved Q icon. The pink and blue were deepened through user review.

The approved source is preserved without further recoloring. Tauri CLI 2.11.5
exports the desktop PNG, Windows ICO and macOS ICNS resources. The 32px PNG is
the favicon; the 128px PNG supplies the 24px header logo, replacing the former
CSS blue dot. The embedded OpenqpView viewer retains its separate identity.
Installed app bundles are not patched. Include these assets in the next Studio
build; public installer publication retains the agreed NMR merge prerequisite.

## Reproduce the platform exports

```sh
npm exec --yes --cache .cache/npm --package @tauri-apps/cli@2.11.5 -- \
  tauri icon docs/design/studio-logo-pink-blue.png --output .cache/generated-icons
```

Copy `32x32.png`, `128x128.png`, `128x128@2x.png`, `icon.png`, `icon.ico`, and
`icon.icns` to `shell/src-tauri/icons/`. Copy `32x32.png` to
`frontend/public/studio-icon.png` and `128x128.png` to
`frontend/public/studio-logo.png`.

## Final image edit prompt

Use case: precise-object-edit. Edit this exact OQP Studio Q app icon. The user asks AGAIN for stronger/deeper colors after a previous too-subtle increase. Make BOTH lobes noticeably richer and more saturated than this reference, not a barely perceptible change. Left lobe: clear medium rich rose pink, main body midtone around #DD6B9B with deeper rosy shadows. Right lobe including Q tail: clear medium azure blue, main body midtone around #559BD7 with deeper blue shadows. Reduce the broad washed-out almost-white highlights so pink and blue remain apparent over most of each surface, while retaining small soft highlights and satin 3D volume. Keep the same pink/blue hues, not red, purple or navy. Change ONLY surface colors and their matching tonal shading. Preserve the Q shape exactly, both organic lobes, outlines, gap, tail, placement, size, smooth satin texture, light direction, dark charcoal rounded-square tile, bevels, square framing, and exterior transparency. No redesign, no added elements, no text or labels. One standalone app icon preview.
