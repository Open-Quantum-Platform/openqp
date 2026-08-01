# ROUTE A: the derivative-sigma amp-channel engine — implementation spec

GOAL: production `analytical_nac()` at the certified formula's accuracy:

    d_IJ = antisym[ G_met . dXt_J  +  gamma_IJ : (Sk + U) ]
         = antisym[ ytil_IJ . w_J  +  gamma-channel ]           (adjoint)

The gamma-channel is DONE (7.32/7.35: z-seam with direct-injected gamma
+ V-mask S^x contractions + gamma:Sk from NAC_DUMP_DS). The remaining
object is w_J^c = (dA/dx_c) X_J contracted with ytil — to be built as a
WHOLE (7.39/7.41: no M:U decomposition).

## The quantity

ytil.w^c = d/dx_c [ ytil^T A(x) X_J ]  at frozen amplitudes, where
A(x) is the MRSF-TDA response operator with FULL implicit dependence:
AO integrals (1e, 2e, XC), C(x) (CPHF), F_AO[D(x)].

This is EXACTLY an interstate excited-state-gradient-type contraction
of the amplitude PAIR (ytil, X_J) — same architecture as the existing
state gradient, with the state's (X,X) bilinears replaced by the
symmetrized (ytil, X_J) bilinears, and the state-specific eliminations
REDERIVED for the bilinear (7.30's chain-polarization failure showed
the existing chain's eliminations assume the eigenpair; do not reuse
them blindly).

## Structure (each term has a certified referee already frozen)

1. SKELETON (DONE, reuse): mrsf_nac_amp + mrsf_nac_esum with ytil
   slot-injected == ytil.w_skel to 1e-5 (gate G-A, 7.30(a)).
2. ORBITAL-RESPONSE: one z-solve per pair with the DIRECT-INJECTED
   RHS X = Mt + gamma (v7i-certified seam == pack(X).U to 1e-3), PLUS
   the elimination terms (V-mask) — CAUTION: this reproduces only the
   J2-level accuracy; the missing Delta (7.39) must come from item 3.
3. THE DELTA CHANNEL (the actual new work): implement the interstate
   W^IJ / relaxed-density contraction FROM SCRATCH per Sec. 4 of the
   derivation, following the STATE gradient's own architecture
   (tdhf_mrsf_z_vector's density builds + mrsfrowcal W + the grd2
   contraction) but with every density rebuilt for the (ytil, X_J)
   bilinear WITHOUT eigen/normalization assumptions:
     - tij/tab -> symmetrized interstate difference densities of
       (ytil, X_J)  [mrsf_interstate_tden already does this]
     - hxa/hxb, ab1_mo terms -> bilinear versions (polarize the
       EXISTING builders at the DENSITY level, not the gradient level)
     - W^IJ -> mrsfrowcal generalization with the bilinear inputs and
       WITHOUT the eigenvector shortcuts; derive each term against the
       per-coordinate referee below.
   REFEREE (frozen, per coordinate, per pair): Delta^c =
   [ytil.(w_ref - w_skel) - Mt:Ux]^c from H2O_energy_tlf0_v7o.npz and
   ETH_energy_v7o.npz(+ctx). A term-by-term Fortran build must
   reproduce Delta; when it does, items 2+3 merge into the standard
   [z.B^x - Tr(W^IJ S^x)] form and the J2 assembly closes at the
   theory level (1e-4..1e-6).
4. REWIRE analytical_nac() (single_point.py): per pair:
   (i) G_met (closed form from the kernel — replace the unit sweep),
   (ii) MINRES ytil on the matvec (certified 2.6e-9),
   (iii) skeleton engines (slot injection),
   (iv) the item-2/3 response contraction,
   (v) gamma-channel (z_gamma seam + V-mask + gamma:Sk),
   (vi) antisymmetrize (I,J)/(J,I); NEVER stamp.
5. VALIDATION LADDER: H2O + ethylene vs frozen d_num (machine target);
   sum rules (translational: proven clean for every certified channel);
   then Acrolein NAMD (ndtlf=2 trajectories, see HANDOFF related-
   session map) as the production shakedown.

## Landmine checklist for the implementer (all documented in 7.24-7.41)
- transpose every 2-D Fortran<->numpy boundary; Python-created arrays
  keep dims verbatim, Fortran-created are reversed.
- in-process phases only; FD sweeps never cross processes (7.36).
- the seam hook antisymmetrizes L internally; seam(e_pq) = -U_full.
- sfrorhs adds 2FT terms absent from the matvec (NAC_ZERO_2FT to
  compare); do not polarize the gradient chain (7.30).
- C2v H2O degeneracies: judge structure only on ethylene.
- int2e_cutoff=1e-20 for matvec linearity.
- G[P] builds: DM push + hf_energy(control.maxit=1); certified
  linear/self-adjoint/null at 1e-10.
