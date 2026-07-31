# OpenQP — Analytical MRSF-TDDFT NAC — Handoff

**Status: WORKING & VALIDATED.** OpenQP's analytical MRSF NAC reproduces KNU-GAMESS's
numerical MRSF NAC on all 6 benchmark pairs, |cos| ≥ 0.99995, norms within 1%.

---

## Address of the final working code

| | |
|---|---|
| **Local path** | `/bighome/alireza/openqp-nac` |
| **Remote (push/pull)** | `git@github.com-alireza-chc2:Alireza-Lashkaripour/openqp-dev-private.git` |
| **Browse** | https://github.com/Alireza-Lashkaripour/openqp-dev-private/tree/nac |
| **Branch** | `nac` |
| **Commit** | `e98d08d` |
| **Paper** | Overleaf `6a29063d942f99866eab0cb4`, working file `main_v2.tex` |

> Standing rule: push ONLY to `devprivate` (this private fork), branch `nac`.
> `upstream` (Open-Quantum-Platform/openqp) is `no-push-allowed`; never touch it or `origin`/main.

---

## What works (the deliverable)

The derivative coupling is `d_IJ = d_ov − d_amp`:
- **`d_ov`** (orbital/Pulay term) — fully closed-form (transition-density kernel + Z-vector).
- **`d_amp`** (amplitude term) — finite difference of the *analytic* MRSF matrix–vector
  product, transported into a fixed reference-orbital frame (a difference of the
  analytic OPERATOR, not of energies/states). Cheap (~490×) because the operator is
  applied only to the transported kets, never built.

Honest label: **semi-analytic** (one FD of the analytic operator remains in `d_amp`).

**Validation (MRSF-BHHLYP/6-31G*, |cos| ≥ 0.99995, ratio within 1%):**
H2O S0/S1, S0/S2, S1/S2; and three S0/S1 MECIs — twisted-pyramidal ethylene (0.107 eV
gap), ethylidene, PSB3. Matches component-by-component, not just in norm.

---

## Key files / entry points

| What | Where |
|---|---|
| Analytical NAC driver | `pyoqp/oqp/library/single_point.py` → `analytical_nac()` (~L2009) |
| The transported-matvec `d_amp` (the working term) | same file → `_compute_amp_damp()` (~L1819) |
| Numerical NAC (in-code) | same file → `numerical_nac()` (~L2203) |
| MRSF matvec (operator A) | `source/modules/tdhf_mrsf_energy.F90` → `mrsf_matvec_apply` |
| Closed-form pieces (esum etc.) | `source/modules/tdhf_mrsf_gradient.F90` |
| Full status write-up | `tools/nac_validation/NAC_STATUS_WRITEUP.md` |
| Validation inputs | `/bighome/alireza/gamess/h2nac/ana/{h2o,eth,ed,psb}_ana.inp` |

---

## Build

MRSF is **ROHF triplet reference + a DFT functional** — there is no "HF-MRSF"; never
run an MRSF input without a functional.

```bash
# LOGIN node (compute nodes lack /opt/conda for the cffi step)
bash tools/nac_validation/build_oqp.sh
# = PATH GCCcore-12.3/bin ; CFLAGS=-g0 ; ninja -C build install ; -> lib/liboqp.so
```

## Run / reproduce the validation

Always on SLURM, never heavy on the login node.

```bash
cd /bighome/alireza/gamess/h2nac
# run_pyoqp.sh sets PYTHONPATH to THIS repo's pyoqp + the correct libs/preload
sbatch -p def -c 16 run_pyoqp.sh <script.py> <input.inp>
```

- The shipped analytical vs numerical comparison uses `ana/*_ana.inp` with
  `[nac] type=analytical`.
- Closed-form diagnostics live in `tools/nac_validation/p15..p38*.py`; the esum FD
  self-tests run under env `NAC_ESUM_FDTEST` (writes `/tmp/nac_esum_{fdtest,xcfd,1efd}.out`).

---

## What is still open (the closed-form refinement)

Removing the last FD from `d_amp` (making it fully closed-form). Not required for the
result above; it does not affect any validated number.

Progress:
- **Proven:** the closed form must be built in the response-operator (matvec) frame.
  The excited-state-gradient route is excluded — a homogeneity test shows the gradient
  is degree-2, so polarization returns `X_I^T M X_J` with M the gradient chain's
  Hessian; M matches the true operator on the eigenvector diagonals but NOT off-diagonal
  (H2O S0/S2 polarization = exactly ½ the correct value).
- **Done + FD-validated (1e-9..1e-5):** the explicit `esum` term
  `Tr[P^IJ · dF^skel/dR]` (1e + 2e + XC), plus the transport term (validated as
  necessary on C1 ethylene: (2,3) −0.54→+0.95, (1,3) +0.24→+0.78).
- **Open:** the interstate orbital response. A free 4-parameter fit over the terms in
  hand leaves 8–88% residuals ⇒ the target is NOT in their span ⇒ a structurally new
  term is needed (the Li–Liu / Fatehi–Subotnik Lagrangian construction, unpublished
  for MRSF). This is a derivation problem, not debugging.

### Landmines (each cost a wrong result; all caught only by an independent FD check)
1. `fock_deriv_contract_os` is `intent(out)` and zeroes `gx` → it OVERWRITES; sum the
   two spin contributions through a scratch array.
2. `utddft_xc_gradient` overwrites, scales `da/db/pa/pb` by basis norms IN PLACE (pass
   copies), and adds a probe-independent ground-state XC gradient → subtract a
   zero-probe baseline.
3. `grd1`'s `logtol` is a LOG threshold; passing `log(10)*int2e_cutoff` with the
   harness value `1e-20` silently screens terms to zero.
4. On a symmetric molecule (H2O, C2v) a 1-D irrep block makes `|cos|=1` automatic and
   meaningless — judge by ratio and test on a C1 system (ethylene).

---

*Deep detail + full derivations: `tools/nac_validation/NAC_STATUS_WRITEUP.md`.*
