# Performance options and the `perf` preset

OpenQP's performance knobs are ordinary **input keys** backed by the shared control
struct — there are no environment variables involved. A single opt-in preset
`[input] perf` bundles them into one accuracy↔speed dial.

```
[input]
perf = 1        # 0, 1, 2, or 3 (default 1; set perf=-1 to disable the preset)
```

## The four levels

| `perf` | Use it for | What it does | Accuracy |
|---|---|---|---|
| **0** | strict reference / reproducibility | every accelerator off, every cutoff tightest | bit-reference (fixed thread count) |
| **1** | **recommended production** | only exact, proven-helpful knobs: MRSF response cutoff `1e-8`, z-vector warm-start (+ always-on Fock digestion) | ≈ reference (≤ µEh) |
| **2** | faster, tiny degradation | `perf=1` + coarse-to-fine XC grid + gradient Schwarz cutoff `1e-8` | gradients within ~5×10⁻⁷ a.u.; SCF exact at convergence |
| **3** | aggressive, **degradation allowed** | looser cutoffs traded for speed: gradient `1e-7` (~10⁻⁵ a.u.), response `1e-6` (~few µeV) | small, controlled — verify for your system |

**`perf=1` is the default** (recommended production, exact). Set `perf=-1` to disable the preset
and leave every knob at its control default. All levels were calibrated
against a CPU benchmark (MKL + Apple Accelerate) across HF/DFT/TDDFT/MRSF energies and MRSF
gradients.

> **Why these and not the rest.** On the CPU builds benchmarked, the reliable speedups are
> the `perf=1` knobs (response cutoff, warm-start, Fock digestion) and, at higher levels,
> *looser integral cutoffs*. The XC Φ-cache, IncDFT, MRSF FP32 and progressive screening
> (`pscreen`) were performance-neutral to *negative* on CPU (FP32 ~2× slower under MKL), so
> **no preset enables them** — they remain available as explicit input keys for regimes the
> CPU benchmark doesn't cover (e.g. GPU XC). `perf=1` is the default; raise it only if your
> own timing shows a gain.

## Individual input keys

Every knob is also a direct input key. Each defaults to `auto` (defer to the preset);
an explicit value **overrides** the preset.

| Section | Key | Meaning | Values |
|---|---|---|---|
| `[scf]` | `xc_c2f` | coarse-to-fine XC grid during SCF descent | `on`/`off`/`auto` |
| `[scf]` | `xc_phi_cache` | cache collocation Φ across SCF iters (exact; opt-in) | `on`/`off`/`auto` |
| `[scf]` | `xc_incdft` | incremental DFT (experimental; opt-in) | `on`/`off`/`auto` |
| `[scf]` | `pscreen` | progressive integral/grid screening (+ `pscreen_cap`) | `on`/`off` |
| `[scf]` | `grad_cutoff` | Schwarz cutoff for the 2e-derivative gradient build | a number, e.g. `1.0d-8`, or `auto` |
| `[tdhf]` | `resp_cutoff` | 2e cutoff for the MRSF response build | a number, e.g. `1e-8`, or `auto` |
| `[tdhf]` | `fp32` | single-precision MRSF response digestion (opt-in) | `on`/`off`/`auto` |
| `[tdhf]` | `zv_warmstart` | reuse the previous step's z-vector as the CPHF seed (exact) | `on`/`off`/`auto` |

Precedence (low → high): control default → `perf` preset → explicit input key. Example —
production speed but keep the full XC grid and tighten the response on one job:

```
[input]
perf = 2
[scf]
xc_c2f = off
[tdhf]
resp_cutoff = 5e-11
```

The resolved settings (and any warnings) are printed near the top of the run log:

```
   Performance settings (perf = 2)
     scf.xc_c2f        = off       (input)
     scf.grad_cutoff   = 1.0d-8    (preset)
     tdhf.resp_cutoff  = 5e-11     (input)
     tdhf.zv_warmstart = on        (preset)
     ...
```

## Safety notes

- `perf=0` is the bitwise reference **for a fixed thread count**. The threaded Fock build
  uses a dynamic-schedule reduction, so results are not bit-identical across thread counts
  at any level — pin `OMP_NUM_THREADS` for strict reproducibility.
- `perf=3` deliberately trades a little accuracy for speed; do not use it for reference
  numbers or tight 2nd-order properties (IR/Raman/NMR) — use `perf ≤ 1` there.
- `fp32` (opt-in) is non-reproducible, can flip near-degenerate excited states, and was
  net-slower than fp64 on CPU; only enable it where you have measured a gain.
- `xc_incdft` (opt-in) is experimental and typically *slows* SCF convergence.
