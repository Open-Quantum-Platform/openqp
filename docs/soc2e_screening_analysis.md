# 2e SOC screening bug: diagnosis and fix

Reference: Komarov, Park, Lee, Zeng, Choi, JCTC 2023, 19, 953 (SOC-MRSF).
Test case: C atom, ROHF/BHHLYP/cc-pVTZ, scal_rel=2, runtype=soc, soc_2e=1.
Published / experimental 3P1-3P0 spacing: 16.4 cm-1.

## Symptom

C-atom 3P SOC spacing 22.19 cm-1 (+36% vs published/exp); Si looked fine.

## Diagnosis (element-level, C/6-31G)

The 2e SOC AO matrix (wao) produced by soc2e_driver was compared element
by element against an independent SOMF implementation (PySCF,
int2e_p1vxp1, vj - 3/2 vk - 3/2 vk^T) using the identical density and a
verified AO permutation:

  * OpenQP wao == -2*vj exactly (every element, 6 decimals): the EXCHANGE
    contributions were entirely absent.
  * A Python emulation of the compute_soc2e_ao accumulation logic (same
    loops, same val2/val3/val4 lines) reproduces the full Eq.-8 result to
    machine precision -- the contraction logic itself is CORRECT.

Root cause: the Schwarz density bound in soc2e_rys_compute used only the
ket-block density max|D(K,L)|.  That is a valid bound for the Coulomb
term, but the exchange terms contract the cross blocks D(I,L), D(J,L),
D(I,K), D(J,K).  Whenever D(K,L) ~ 0 while cross blocks are finite, the
whole quartet is skipped and its exchange contributions are silently
lost.  For a spherical atom the s-p density blocks are exactly zero, so
ALL exchange contributions vanished, leaving pure Coulomb screening.
For general molecules the density blocks are dense and the loss is tiny
(H2O / CH3Br example soc_eval shifts by only ~2e-6 with the fix; the
bundled reference JSONs remain valid).

## Fix

soc2e_rys_compute now takes dabmax as the max |D| over all five density
blocks used by the contraction (KL, IL, JL, IK, JK).

## Validation after the fix

| quantity                          | before | after | reference        |
|-----------------------------------|--------|-------|------------------|
| C  3P spacing (cm-1)              | 22.19  | 16.28 | 16.4 exp/published |
| C  screening ratio (2e/1e)        | 0.734  | 0.538 | 0.535 PySCF Eq.8 |
| Si 3P spacing (cm-1)              | 76.53  | 69.22 | 77.1 exp; ~75.7 published |
| Si screening ratio                | 0.854  | 0.772 | 0.773 PySCF Eq.8 |
| H2O / CH3Br examples              | pass   | pass  | bundled JSONs    |

After the fix OpenQP agrees with the independent PySCF Eq.-8 (SOMF)
implementation for both atoms.  C also agrees with experiment and the
published value.  The remaining -10% on Si relative to the published
GAMESS-based 75.7 cm-1 is now a question of WHICH two-electron treatment
the published calculations used (full Breit-Pauli 2e vs. the mean-field
Eq. 8; GAMESS offers HSO2/HSO2FF beyond the mean-field) -- to be
confirmed with the authors.  The screening density (closed shells,
occupation 1, with the factor 2 inside the kernel) is consistent with
Eq. 8 for the closed-shell part; including the open-shell density
changes C/Si by <0.5% (PySCF check) and is not the source of any
remaining deviation.
