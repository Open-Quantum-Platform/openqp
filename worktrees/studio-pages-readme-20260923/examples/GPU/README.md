# Optional GPU backend (`OPENQP_WITH_GPU`)

`OPENQP_WITH_GPU` is a **build-time** option (not a runtime input flag) that
links the optional CUDA library `openqp-gpu` and compiles the `routec_sig` /
`routec_bridge` seams in their linked path. It is **off by default**; a stock
build is unaffected.

The seams are **transparent and self-falling-back**: the input decks are
identical whether or not the GPU is present, and any decline (missing library,
unsupported case, or a device error) reverts to the native path. So the GPU
backend is exercised by the *existing* examples — it just runs their heavy
two-electron work on the device:

| GPU seam | Exercised by |
|----------|--------------|
| MRSF-TDDFT sigma (energy) | `examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.{inp,json}` |
| MRSF-TDDFT 2e gradient | `examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_GRADIENT.{inp,json}` |
| SCF J/K + DFT Vxc | `examples/DFT/*` and `examples/HF/*` |
| ground-state 2e gradient | `examples/OPT/*`, `examples/HF/*` gradient decks |

## Coverage classification

This is a **run-only optional backend** (in the sense of AGENTS.md rule 2, like
other optional backends): CI runners have no GPU, so the device path cannot be
exercised in CI. The committed `.json` references for the examples above are
validated by CI on the **native** path (GPU off). On a machine with the GPU
backend built, the same decks reproduce those references to density-fitting
accuracy (~1e-6 Ha) — that is the GPU regression check.

To run an example on the GPU:

```bash
# build once with the backend (needs a CUDA toolkit)
cmake -B build -DOPENQP_WITH_GPU=ON -DCMAKE_CUDA_ARCHITECTURES=80
cmake --build build
# then run any MRSF/HF/DFT/gradient example as usual; it diverts to the device
openqp examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.inp
# per-run opt-out (fall back to native): OQP_ROUTEC_SIG=0 (and OQP_ROUTEC_*_LIB=0)
```
