# openqp-gpu

GPU-accelerated density-fitting HF, DFT, and MRSF-TDDFT (energy + gradient) for
[OpenQP](https://github.com/Open-Quantum-Platform/openqp).

A self-contained CUDA library. Given the density-fitting tensor, the one-electron
matrices, and a starting guess, it runs the entire calculation — Coulomb and
exchange builds, exchange-correlation, Fock assembly, diagonalization, DIIS, the
MRSF excited-state sigma operator, and the nuclear gradient — resident on a single
GPU. The SCF loop depends only on cuBLAS and cuSOLVER; there is no external
quantum-chemistry package in the loop.

## Status

| Path | State |
|------|-------|
| Restricted / hybrid density-fitting SCF (`routec_scf_solve`) | **working, validated** — energies match reference to 1e-10 |
| DFT exchange-correlation (RKS/UKS, `routec_vxc`) | **working, validated** — exact to 1e-14 on the XC energy |
| MRSF-TDDFT excitation energies (`routec_sig_*`, the OpenQP `OQP_ROUTEC_SIG` seam) | **working, validated** — matches native OpenQP to 2e-6 (density-fitting error) end-to-end |
| MRSF-TDDFT nuclear gradient (`routec_grad2_mrsf`) | **working, validated** — matches the CPU reference to 1e-11 relative |
| Ground-state nuclear gradient, 2e part (`routec_grad2`) | **working, validated** — matches the CPU reference to 1e-11 relative |
| Native DF-tensor builder (`openqp_gpu_build_df`) | **working, validated** — fully pyscf-free tensor drives the MRSF seam to 6e-6 (density-fitting error) end-to-end |

## Performance

All numbers on 1× A100-PCIE-40GB, cc-pVDZ, density fitting. Two baselines are
reported and labelled: **gpu4pyscf** (pyscf's GPU code — a GPU-vs-GPU comparison)
for the ground-state SCF and DFT, and **native OpenQP** (the shipping CPU code,
16 threads) for MRSF, which has no GPU baseline to compare against.

### Ground-state SCF and DFT — vs gpu4pyscf 1.7.1 (GPU vs GPU)

Energies identical to the pyscf result to 1e-10. Wall-to-wall SCF time.

| Calculation                     | openqp-gpu | gpu4pyscf | speedup |
|---------------------------------|-----------:|----------:|--------:|
| Hartree–Fock, (H2O)16           |   0.19 s   |   0.61 s  | **3.2x** |
| Hartree–Fock, (H2O)24           |   0.66 s   |   3.04 s  | **4.6x** |
| DFT, pure GGA (BLYP)            |   0.565 s  |   0.697 s | **1.2x** |
| DFT, hybrid (BH&HLYP)          |   0.473 s  |   0.854 s | **1.8x** |

The exchange-fraction gradient is expected: gpu4pyscf pays ~50% more for the
exact-exchange in a hybrid, while our occupied-orbital exchange makes it nearly
free — so the advantage widens from the pure functional (1.2x) to the hybrid
(1.8x) to pure Hartree–Fock (3–5x).

### MRSF-TDDFT — vs native OpenQP (GPU vs CPU)

The excited-state solver's per-iteration work (the sigma operator) and the MRSF
gradient's two-electron assembly are the expensive parts; both move to the GPU.

| Calculation                                        | native CPU | openqp-gpu | speedup |
|----------------------------------------------------|-----------:|-----------:|--------:|
| MRSF sigma, per Davidson iteration, (H2O)8         |   2.41 s   |   0.13 s   | **18x** |
| MRSF sigma, per Davidson iteration, (H2O)16        |  16.00 s   |   0.68 s   | **23x** |
| MRSF sigma, per Davidson iteration, (H2O)24        |  40.56 s   |   3.05 s   | **13x** |
| MRSF gradient, total, (H2O)8                       |  11.4 s    |   2.3 s    | **5x**  |
| — its two-electron exchange assembly alone         |   9.4 s    |   0.16 s   | **59x** |

MRSF excitation energies through the live OpenQP Davidson reproduce the native
exact-integral result to density-fitting accuracy at every size (3.9e-6, 2.0e-6,
6.7e-6 Ha for (H2O)8/16/24), using a fully pyscf-free tensor.

MRSF excitation energies through the live OpenQP Davidson solver reproduce the
native exact-integral result to 2e-6 Ha (the density-fitting error) at (H2O)16.

### Big-system memory — compressed density fitting (CDF)

The density-fitting tensor is what limits system size: it grows fast and fills
GPU memory. With `OQP_CDF_ONDEV` the library holds it compressed on the device
and rebuilds small slices on the fly, so the full dense tensor never resides.
Two levers, both opt-in, both bit-identical by default, both kept below the
density-fitting error:

- **pair compaction** — drop the tensor rows that are exactly zero (lossless);
- **precision tiers** (`ROUTEC_CDF_LOWPREC`) — store each auxiliary column in
  64/32/16/8-bit by its magnitude. 32-bit is the safe workhorse (about 2x
  smaller with the energy far below the fitting error); 16/8-bit help further
  where the magnitude spread is wider.

Peak GPU memory, Hartree–Fock, same geometry family:

| System  |   dense | + compaction | + 32-bit tier | vs dense |
|---------|--------:|-------------:|--------------:|---------:|
| (H2O)16 | 5281 MB |      2205 MB |       1855 MB | **2.9x** |
| (H2O)24 | 15773 MB |     4043 MB |       3171 MB | **5.0x** |

Energies stay within the density-fitting error (at (H2O)24 the precision tiers
add 3.5e-6 Ha, against the ~1e-3 fitting error). The saving grows with system
size, so a 40 GB card that runs out near 50 waters on the dense tensor reaches
substantially larger clusters.

### Ground-state SCF/DFT through the live OpenQP seam

Running HF and DFT inside OpenQP with the GPU seam active (its Fock J/K and DFT
Vxc built on the device from a pyscf-free tensor) reproduces the correct
density-fitting energy and replaces OpenQP's native exact-integral SCF:

| Calculation, (H2O)8 | native SCF (exact ERI) | GPU seam (DF) | speedup | dE |
|---------------------|-----------------------:|--------------:|--------:|---:|
| Hartree–Fock        | 8.59 s | 0.89 s | **9.7x** | 4.0e-4 Ha (DF) |
| DFT (BH&HLYP)       | 9.27 s | 1.48 s | **6.3x** | 9.1e-5 Ha (DF) |

The GPU HF energy matches an independent pyscf DF-RHF reference to 2.5e-5 Ha
(the two implementations' DF-fit difference), confirming the seam is correct at
density-fitting accuracy; the difference from native OpenQP is the DF error
itself. (This native-exact vs GPU-DF comparison reflects the total gain from
enabling the seam — switching to density fitting *and* the GPU; the 3–5x above
is the pure GPU-vs-GPU number at fixed method.)

## Independence (no pyscf)

The entire path — library, tensor builder, and prep — is pyscf-free, validated
end to end: a tensor built with no pyscf drives the live OpenQP MRSF Davidson to
**6e-6 Ha** (the density-fitting error) at (H2O)8, identical quality to the
pyscf-built tensor.

- **Integrals.** The builder computes the 3-center integrals (μν|P) *and* the
  2-center auxiliary metric (P|Q) on the GPU with its own rotation kernels
  (`ROUTEC_NATIVE_METRIC`). The metric reconstructs identical electron-repulsion
  integrals to 4.3e-11 against the pyscf metric.
- **Frame.** Orbital shells come from OpenQP's own `get_basis()` (its exact
  internal AO convention), so the tensor is guaranteed to be in the same frame
  as the wavefunction — removing the pyscf-frame calibration bridge that caused
  earlier stale-tensor errors. The cartesian d-component reordering between the
  builder and OpenQP is fixed and baked into the prep (`tools/dfbuild/native_prep.py`).
- **Auxiliary basis** comes from `basis_set_exchange` (not pyscf).

Workflow: `native_prep.py NX NY NZ <dir>` (OpenQP run → shells + operands), then
`ROUTEC_CART=1 ROUTEC_NATIVE_METRIC=1 ROUTEC_B_SIGMA=1 openqp_gpu_build_df ...`
writes a tensor the MRSF/SCF loaders read directly.

## Layout

- `src/scf.cu` — density-fitting HF/hybrid SCF + J/K engine. cuBLAS + cuSOLVER.
- `src/xc.cu` — exchange-correlation kernels (`routec_vxc`). Self-contained.
- `src/sigma.cu` — MRSF-TDDFT sigma-session (the `OQP_ROUTEC_SIG` seam engine).
- `src/grad.cu` (+ `routec_grad_gen*`) — two-electron nuclear gradient (ground
  state and MRSF), Gamma/gamma assembly on the device.
- `tools/dfbuild/` — the native density-fitting tensor builder
  (`openqp_gpu_build_df`): 3-center integrals + 2-center metric on the GPU.
- `seam/mrsf/` — the OpenQP-side seam files (`routec_sig.F90`,
  `tdhf_mrsf_energy.F90`) that divert the MRSF Davidson solver to this library.
- `data/routec_tables.bin` — tabulated radial/Gaunt data for the builder.

The SCF loop and the XC kernels are rotation-free; the rotation codegen appears
only in the gradient's derivative integrals and the tensor builder — both
one-time or optional, both isolated, and both where the rotation method is
genuinely fastest (high angular momentum).

## Build

```
cmake -B build -DCMAKE_CUDA_ARCHITECTURES=80
cmake --build build
```

Requires CUDA (tested 12.6), cuBLAS, cuSOLVER; a host LAPACK/BLAS for the
gradient and the native metric. A100 = `sm_80`.

## Provenance

Extracted and rebuilt from the OpenQP GPU research driver (June 2026): the
device-resident DF-SCF engine with occupied-orbital exchange and on-device CDIIS,
the exchange-correlation kernels, the MRSF sigma-session, and the two-electron
gradient.

## Consolidated direct-integral GPU backend

The audited OpenQP/private GPU branches are consolidated here, with a separate
`openqp_gpu_metc` C/CUDA library for direct-integral contraction, persistent
workspaces, Python planning, and numerical regression tests. The existing
DF-SCF/gradient engine remains independent. See
[the integration record](docs/integration/README.md) for source branches,
component locations, A100 validation, and remaining host/experimental boundaries.
