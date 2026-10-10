# R3 e2e seam — wire the GPU sigma-session into the live MRSF Davidson — RESULTS

Lane: chc4 (build + GPU A100 GPU1) + local Mac (CPU-stub design/port). 2026-06-17.
No commit. Followed DESIGN.md (sessions/20260617_r3_e2e_design/DESIGN.md).

## Status: DONE. G-sigma3 PASS (both CPU-stub and real GPU v3), both multiplicities.

The `OQP_ROUTEC_SIG` seam is implemented, built into liboqp, and validated
end-to-end: the live OpenQP MRSF Davidson runs entirely through the device
sigma-session (whole 6a/6b/6c per-vector triple replaced) and reproduces the
native exact-ERI excitation energies to ~1.7e-11 Ha.

## Headline numbers (stock H2O BH&HLYP/6-31G* MRSF, ROHF mult-3, nstate=3)
System is SPHERICAL: nbf=19, nocca=6, noccb=4, ntrial = nocca*(nbf-noccb) = 90.
B = exact-Cholesky `~/blane_accept/B_exact_dense.bin` (naux=165=19*20/2, dense),
which reproduces native exact-ERI physics (so this is the TIGHT gate, not DF).

| leg                              | singlet exc (Ha)                          | triplet exc (Ha)                       |
|----------------------------------|-------------------------------------------|----------------------------------------|
| native (no env, exact ERI)       | -0.283514008  0.045299858  0.103564966    | 0.030213436  0.094801444  0.100958330  |
| CPU-stub sig-session (OQP_ROUTEC_SIG) | -0.283514007 0.045299858 0.103564966  | 0.030213436  0.094801444  0.100958330  |
| GPU v3 sig-session (A100, GPU1)  | -0.283514007  0.045299858  0.103564966    | 0.030213436  0.094801444  0.100958330  |

- **G-sigma3 (GPU-sig vs native):** S 1.69e-11 Ha, T 1.70e-11 Ha — gate <=1e-6 -> PASS (~5 orders margin).
- **G-sigma3 tight (GPU-sig vs CPU-stub, same B/same math):** T 3.9e-15 (machine), S 1.69e-11 (shared convergence-level offset, not an engine delta) — gate <=1e-10 -> PASS.
- **Native bit-match:** native exc reproduce the frozen 2026-06-12 reference vector exactly (last digit of singlet state-2 differs by 1e-9 = tdconv noise). Proves the seam-built liboqp is correct on the INERT path.
- **Fallback discipline:** OQP_ROUTEC_SIG pointing at a missing .so logs the dlopen
  failure and falls cleanly to native (exact native exc) — seam is inert-safe.

## G-sigma1 (per-iteration plumbing) — demonstrated via the e2e gate + the stub port
The CPU stub (libroutec_sig_cpustub.so) is a verbatim C++ port of sig_ref.py's
validated Stage-0 pipeline (mrsfcbc/mrsfmntoia/mrsfesum + DF-B J/K). Standalone
port check (check_cpustub.py vs sig_ref.build_sigma_stage0): rel 5e-15..1e-14,
S and T, nbf in {24,25,40}. Because the exact-B sig leg matches native to
1.7e-11, the per-iteration amo handed back through the seam equals native amo to
that level — G-sigma1 is closed transitively (any per-iter amo error would
propagate to the excitation energies, which it does not).

## What was implemented
1. **`source/routec_sig.F90`** (new module, clone of routec_bridge dlopen/dlsym):
   `routec_sig_available / _begin / _apply / _end`; iso_c_binding iface to the v3
   ABI (init/set_scale/iter/free, all scalars by pointer; kind=merge(3,1,mrst==3)).
   Reads $OQP_ROUTEC_SIG. Globbed into liboqp automatically (source/*.F90).
2. **`source/modules/tdhf_mrsf_energy.F90`** (2 edits + use + decls):
   - pre-loop gate (after fa/fb built, after scale_exch): `use_sig` decided ONCE,
     active only for mrst in {1,3}, .not.umrsf, all spc_*==HFscale, no CAM, and
     `routec_sig_begin(nbf, mo_a, mo_b, fa, fb, nocca, noccb, mrst, scale_exch)==0`.
   - in-loop branch: `if (use_sig)` -> `routec_sig_apply(bvec_mo(:,ist:iend), nv,
     amo(:,ist:iend))` (mid-run decline = hard abort); `else` = the verbatim
     native 6a/6b/6c block. Subspace solver (rparedms/rpaeig/.../rpanewb) untouched.
   - post-loop `if (use_sig) call routec_sig_end()` next to int2_driver%clean().
3. **`routec_sig_cpustub.cpp`** (CPU validation stub, the v3 ABI via exact Stage-0
   math + DF-B from the same B-file format the GPU engine reads). The one porting
   bug found+fixed: the ball CV-block contraction `Ca_o @ tmp.T` had a/b swapped.

## MO-Fock handoff (the HIGH risk) — CONFIRMED correct
The seam hands the unpacked-square MO Fock `fa`/`fb` (from
orthogonal_transform_sym(fock_*,mo_*)+unpack_matrix, driver ~369-388) to
sig_init's fmo_a/fmo_b. If this were the wrong matrix (e.g. raw packed fock_a),
mrsfesum's orbital-energy diagonal terms would be wrong and excitations grossly
off. The 1.7e-11 match proves fa/fb are handed correctly. (ixcore level-shift:
flagship has no ixcore, so fa/fb are the clean MO Fock; guarded/irrelevant here.)

## Files
Local (this session wd/):
- routec_sig.F90              the new bridge module
- tdhf_mrsf_energy.F90        the patched driver (2 edits)
- routec_sig_cpustub.cpp      CPU validation stub (v3 ABI, exact Stage-0)
- check_cpustub.py            stub-vs-sig_ref port check (PASS 5e-15)
- run_sig.py                  the e2e run harness (native / sig legs, mult 1/3)
- cmp.py, cmp2.py             delta tabulators
- libroutec_sig_cpustub.so    locally-built stub (Mac)
chc4:
- ~/openqp_routec/source/routec_sig.F90  (built into liboqp; backup: none needed, new)
- ~/openqp_routec/source/modules/tdhf_mrsf_energy.F90  (backup .pre_sig kept)
- ~/openqp_routec/build_chc4/source/liboqp.so  (rebuilt WITH seam; 4 routec_sig syms)
- ~/r3_e2e_seam/{libroutec_sig_cpustub.so, run_sig.py, env.sh, oqp_root_lane/, *.npz, *.log}
- B: ~/blane_accept/B_exact_dense.bin (nbf=19 exact-Cholesky, reused)
- GPU engine: ~/r3_sigsession/libroutec_oqp_gpu2_v3.so

## Run protocol (chc4)
```
source ~/r3_e2e_seam/env.sh                 # OPENQP_ROOT=oqp_root_lane (seam liboqp)
cd ~/r3_e2e_seam
$DRIVER run_sig.py native 1 --nstate 3      # gold reference
CUDA_VISIBLE_DEVICES=1 $DRIVER run_sig.py sig 1 --nstate 3 \
   --so ~/r3_sigsession/libroutec_oqp_gpu2_v3.so --b ~/blane_accept/B_exact_dense.bin
$DRIVER cmp2.py                              # tabulate deltas
```

## Restore (no commit; lane is on COPIES with backups)
`cp ~/openqp_routec/source/modules/tdhf_mrsf_energy.F90.pre_sig
    ~/openqp_routec/source/modules/tdhf_mrsf_energy.F90`
and `rm ~/openqp_routec/source/routec_sig.F90`, then rebuild to revert.

## Honest notes / open items
- The machine-zero gate used an EXACT-Cholesky B (B_exact_dense.bin), so the
  ~1.7e-11 IS the tight gate (DF-fit error is zeroed). A coarse DF-B
  (def2-universal-jkfit) leg would land at the ~1e-5 DF level (G-sigma4); not run
  here because pyscf is banned (the only on-hand DF-B exporter used pyscf) and
  the exact-Cholesky B already gives the stronger result.
- Validated on the spherical flagship (nbf=19). Cartesian (nbf=25) would need a
  matching B; the seam dims (ntrial via iatogen) are basis-agnostic.
- spc!=HFscale / CAM / quintet / UMRSF are guarded to native (engine unsupported);
  the pre-loop gate enforces this. Not exercised (flagship needs none).
- The seam edits sit in the chc4 working tree on COPIES with backups; NOT committed.

## 2026-07-04 — revalidated against the rebuilt engine in openqp-gpu

The lost June GPU sigma dylib was re-implemented as src/sigma.cu in this
repository (low-rank cuBLAS J/K; machine-precision gate vs the CPU stub,
worst rel 8.7e-17). The live OpenQP MRSF Davidson (chc4 seam build,
openqp_routec + routec_sig.F90) pointed at libopenqp_gpu.so with
B_exact_dense.bin reproduces native exact-ERI excitations:

    mult=1 max|dE| = 1.687e-11 Ha   dE_scf = 5.7e-14
    mult=3 max|dE| = 1.698e-11 Ha   dE_scf = 4.3e-14

identical margins to June. MRSF energy through openqp-gpu: CORRECT.

## 2026-07-05 — production-scale MRSF energy: correct AND 21x

Engine correctness at scale (tools/sigval/sig_gpu_validate.py vs the CPU stub
referee, real (H2O)16 dims nbf=400 nocca=81 noccb=79): worst rel 1.14e-16 —
machine precision. Warm GPU sigma iteration 72 ms.

Live OpenQP MRSF Davidson through the OQP_ROUTEC_SIG seam, (H2O)16 cc-pVDZ
BH&HLYP, ROHF mult-3 ref, singlet nstate=6, A100:

    per-iter sigma:  native 16.94 s/it  vs  GPU-session 0.79 s/it  = 21x
    SCF identical (-1221.967448033); excitations vs native = 2.06e-6 Ha (DF err)

Fresh DF-B built by stage2_gpu_inloop_e2e.py (grid 2 2 4, def2-universal-jkfit,
naux=1808). Validated at (H2O)2 (1.0e-5), (H2O)8 (5.8e-6), (H2O)16 (2.1e-6).

CAUTION: the June mrsf_e2e/B_oqp.bin (Jun 11) is STALE — driving the session
with it gives ~1 Ha-wrong excitations. Its frame calibration had degraded;
always rebuild B with the current exporter. The clean long-term fix is the
native openqp-gpu build_df tool fed an OpenQP shell dump (no pyscf frame bridge).
