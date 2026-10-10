# Wiring openqp-gpu into OpenQP

Two files here are the OpenQP-side seam that diverts the MRSF-TDDFT Davidson
solver to this GPU library:

- `routec_sig.F90` — the bridge module (`use routec_sig`). Dual-mode: a
  compile-time link (`-DOQP_GPU_LINKED`) or a runtime `dlopen`.
- `tdhf_mrsf_energy.F90` — OpenQP's MRSF energy driver, patched to call the
  bridge inside the Davidson loop (`routec_sig_available/begin/apply/end`).

`openqp_gpu.cmake` is the build glue: an `OPENQP_WITH_GPU` option that
auto-downloads, builds, and links the library.

## Steps

1. Copy `routec_sig.F90` into OpenQP's `source/` and build it into `liboqp`.
   Apply the `tdhf_mrsf_energy.F90` changes (the `routec_sig` `use` and the
   `routec_sig_available()` branch around the sigma triple; see the file here).

2. In OpenQP's top-level `CMakeLists.txt`:

   ```cmake
   include(${CMAKE_SOURCE_DIR}/source/openqp_gpu.cmake)   # wherever you put it
   # ... after the liboqp target exists:
   openqp_gpu_attach(oqp)                                 # oqp = the liboqp target
   ```

3. Build openqp-gpu **separately** (its own repo/build) and install it, then
   build OpenQP pointing at it:

   ```
   # in the openqp-gpu checkout, once:
   cmake -B build -DCMAKE_INSTALL_PREFIX=$HOME/opt/openqp-gpu
   cmake --build build && cmake --install build

   # then in OpenQP:
   cmake -B build -DOPENQP_WITH_GPU=ON \
         -DCMAKE_PREFIX_PATH=$HOME/opt/openqp-gpu \
         -DCMAKE_CUDA_ARCHITECTURES=80
   cmake --build build
   ```

## What the option does

openqp-gpu is a **separate** project; OpenQP does not vendor or download it, so
the public OpenQP tree has no dependency on that repository.

| `OPENQP_WITH_GPU` | Behaviour |
|-------------------|-----------|
| `ON`  | Links the separately built openqp-gpu -- found via `find_package(openqp_gpu)` on `CMAKE_PREFIX_PATH`, or built in-tree from `-DOPENQP_GPU_SOURCE_DIR=<checkout>`. Compiles the seams with `-DOQP_GPU_LINKED`; the GPU path runs with **no runtime env var**. Needs a CUDA toolkit. Errors with instructions if the library is not found. |
| `OFF` (default) | Nothing is linked; the seams keep their runtime behaviour (inert unless the `OQP_ROUTEC_*` env vars name a dylib). OpenQP builds exactly as before. |

So the GPU backend is fully separable: whether openqp-gpu is public or private,
OpenQP itself never references it. Those with the library plug it in (installed
package, local source, or runtime dlopen); everyone else builds stock OpenQP.

## Runtime control

- **Linked build**: on by default. Opt out per run with `OQP_ROUTEC_SIG=0`
  (also `off` / `none`) to fall back to OpenQP's native exact-integral path.
- **dlopen build**: off by default. Point `OQP_ROUTEC_SIG` at
  `libopenqp_gpu.so` to enable it.

Either way, when the session is active you get one stderr line
`[routec_sig] MRSF Davidson diverted to GPU sigma-session`, and the density
tensor is read from `OQP_ROUTEC_B` (see the main README for the pyscf-free
builder that produces it).

## Validated

- `routec_sig.F90` compiles in both modes; in the linked build the symbols
  resolve at link time and `routec_sig_available()` is true with no env var; the
  `OQP_ROUTEC_SIG=0` opt-out disables it.
- `openqp_gpu.cmake` with `OPENQP_WITH_GPU=ON` configures the subproject, links
  `libopenqp_gpu`, and sets `-DOQP_GPU_LINKED`; `OFF` is a clean no-op.
- End-to-end through the live OpenQP MRSF Davidson (dlopen build): excitation
  energies match native to 2e-6 at (H2O)16 (see `RESULTS.md`).
