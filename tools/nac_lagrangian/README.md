# nac-lagrangian campaign harnesses (2026-07-31)

Clean-room derivation + gates for the fully analytic MRSF-TDDFT NAC.
Every claim here is backed by a numeric gate; see MRSF_NAC_DERIVATION.md
for the derivation and the campaign log (sections 7.x).

| script | what it gates |
|---|---|
| conv_check.py | numerical-NAC conventions: d antisym, h sym, signed gap orientation |
| freeze_ref.py | dumps frozen reference dcv/nacv vectors (npz) |
| gamma_gate.py | raw-determinant gamma: closed form == Slater-Condon == exact biorthogonal overlap response (A/B/C/D gates) |
| formula_gamma.py | EXACT replica of compute_states_overlap (2e-16) + gamma^formula extraction by generator sweep + FD gate |
| orientation_gate.py | ABSOLUTE d_ij-vs-d_ji orientation vs a code-independent exact-overlap oracle |
| assembly_gate.py | the master decomposition d_num = Xt_I.dXt_J + gamma^formula:T, T = dM/dx, vs production numerical NAC (signed) |

Conventions established (any deviation is a bug):
- dcv[i,j] = d_ij = <I|d/dR|J> antisym; nacv[i,j] = (E_j - E_i) d_ij sym.
- 2-D tagarray buffers reach numpy TRANSPOSED; transpose at every read.
- Davidson state phases are random per run; compare gauge-resolved with
  the pair-sign product rule (product over a state cycle = +1).
- tlf=0 (ov_exact, fixed) = exact overlap; tlf=2 error is ~1e-10 in
  geometric FD at dx=1e-3.
