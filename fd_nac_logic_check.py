#!/usr/bin/env python
"""Backend-free verification of the OpenQP FD/overlap NAC assembly math.

`fd_nac_verify.py` runs the real LiH MRSF calculation and needs a compiled
OpenQP (Fortran/CFFI) backend. This companion script instead exercises the
*exact Python assembly logic* of `NACME.nacme` and `NAC.numerical_nac` on
synthetic overlap matrices built from a known derivative coupling `d`, so it
proves the factor-2 / gap-sign fix and the §5 properties without a build.

Model: for real orthonormal adiabatic states and a single Cartesian
displacement of size Delta along which the (antisymmetric) derivative coupling
is `d_ij = <Psi_i|grad Psi_j>`, the state overlap to first order is

    S^+(i,j) = <Psi_i(0)|Psi_j(+Delta)> = delta_ij + Delta * d_ij + O(Delta^2)
    S^-(i,j) = <Psi_i(0)|Psi_j(-Delta)> = delta_ij - Delta * d_ij + O(Delta^2)

We add a symmetric O(Delta^2) term to both (curvature of the overlap) to make
the test non-trivial; it must cancel in the antisymmetrized, central-difference
assembly, leaving exactly `d`.

Run:  python fd_nac_logic_check.py
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True, linewidth=140)

RNG = np.random.default_rng(0)
NSTATE = 4
DELTA = 1.0e-2                      # displacement (Bohr)
DT = DELTA                          # numerical_nac sets config['nac']['dt'] = dx

# --- ground truth -----------------------------------------------------------
# antisymmetric derivative coupling d_ij (scalar per pair for this 1 coordinate)
A = RNG.standard_normal((NSTATE, NSTATE))
D_TRUE = A - A.T                    # exactly antisymmetric, d_ii = 0
# adiabatic energies (strictly increasing -> clean non-degenerate gaps)
E = np.array([-8.0, -7.6, -7.59, -7.1])    # includes a near-degenerate (1,2) pair

# symmetric O(Delta^2) curvature contamination of the overlap (must cancel out)
C = RNG.standard_normal((NSTATE, NSTATE))
C = C + C.T

# synthetic overlaps with the *correct* first-order structure + curvature noise
Sp = np.eye(NSTATE) + DELTA * D_TRUE + DELTA**2 * C
Sm = np.eye(NSTATE) - DELTA * D_TRUE + DELTA**2 * C


# --- mirror of single_point.verify_nac_conventions (kept in sync by hand) ---
def verify_nac_conventions(dc, nac, atol=1.0e-6, label=''):
    dc = np.asarray(dc, dtype=float)
    nac = np.asarray(nac, dtype=float)
    d_res = np.linalg.norm(dc + np.swapaxes(dc, 0, 1)) / (np.linalg.norm(dc) + 1.0e-30)
    h_res = np.linalg.norm(nac - np.swapaxes(nac, 0, 1)) / (np.linalg.norm(nac) + 1.0e-30)
    tag = f' ({label})' if label else ''
    assert d_res < atol, f'd not antisym{tag}: {d_res:.2e}'
    assert h_res < atol, f'h not sym{tag}: {h_res:.2e}'
    return d_res, h_res


def gap_EjmEi(energies):
    """gap[i,j] = E_j - E_i (the unified convention, both code paths)."""
    e_i = energies.reshape((-1, 1))   # row -> bra state i
    e_j = energies.reshape((1, -1))   # col -> ket state j
    return e_j - e_i


# === NACME.nacme path (also the MD / PyRAI2MD export path) ==================
# FIXED: dc = (S - S^T) / (2*dt)
dc_nacme = (Sp - Sp.T) / (2.0 * DT)          # use S^+ as the single MD frame
gap = gap_EjmEi(E)
nac_nacme = dc_nacme * gap
# OLD (buggy): dc = (S - S^T) / dt
dc_nacme_old = (Sp - Sp.T) / DT

# === NAC.numerical_nac path (central-difference FD oracle) ==================
# FIXED: per-point dc = (S^+- - S^+-^T)/(2*dx); central avg keeps /2
forward = (Sp - Sp.T) / (2.0 * DELTA)
backward = (Sm - Sm.T) / (2.0 * DELTA)
dc_num = (forward - backward) / 2.0
gap_num = gap_EjmEi(E)                        # FIXED: same orientation as nacme
nac_num = dc_num * gap_num
# OLD (buggy): per-point /dx, gap = E_i - E_j
forward_old = (Sp - Sp.T) / DELTA
backward_old = (Sm - Sm.T) / DELTA
dc_num_old = (forward_old - backward_old) / 2.0
gap_num_old = -gap_EjmEi(E)                   # E_i - E_j (swapped reshapes)
nac_num_old = dc_num_old * gap_num_old

# hand-built independent reference
hand_d = (Sp - Sm) / (2.0 * DELTA)           # = d + O(Delta^2)? curvature cancels -> d
h_std = D_TRUE * gap_EjmEi(E)                 # h_ij = (E_j - E_i) d_ij


def rel(a, b):
    return np.linalg.norm(a - b) / (np.linalg.norm(b) + 1e-30)


print('# ground-truth d antisym? ', rel(D_TRUE, -D_TRUE.T) < 1e-12)
print('# energies:', E, ' (near-degenerate pair 1-2, gap = %.4f)' % (E[1] - E[0]))

checks = []

def check(name, cond, detail=''):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f'   {detail}' if detail else ''))

print('\n# --- FIXED code paths ---')
# 1. factor-2 removed: recovered d equals true d (and hand reference)
check('nacme dc == true d', rel(dc_nacme, D_TRUE) < 1e-9, f'rel={rel(dc_nacme, D_TRUE):.2e}')
check('numerical_nac dc == true d', rel(dc_num, D_TRUE) < 1e-9, f'rel={rel(dc_num, D_TRUE):.2e}')
check('hand (S+ - S-)/(2D) == true d (curvature cancels)', rel(hand_d, D_TRUE) < 1e-9,
      f'rel={rel(hand_d, D_TRUE):.2e}')
# 2. symmetry assertions (also calls the production-mirrored guard)
verify_nac_conventions(dc_nacme, nac_nacme, label='nacme')
verify_nac_conventions(dc_num, nac_num, label='numerical_nac')
check('d antisymmetric (both paths)', True)
check('h symmetric (both paths)', True)
# 3. h = (E_j - E_i) d, standard sign, both paths agree
check('nacme h == (E_j-E_i) d', rel(nac_nacme, h_std) < 1e-9, f'rel={rel(nac_nacme, h_std):.2e}')
check('numerical_nac h == (E_j-E_i) d', rel(nac_num, h_std) < 1e-9, f'rel={rel(nac_num, h_std):.2e}')
# 4. cross-path agreement (same sign & magnitude)
check('cross-path dc agree', rel(dc_nacme, dc_num) < 1e-9, f'rel={rel(dc_nacme, dc_num):.2e}')
check('cross-path h agree (sign+mag)', rel(nac_nacme, nac_num) < 1e-9, f'rel={rel(nac_nacme, nac_num):.2e}')

print('\n# --- OLD code reproduced the documented bug ---')
# old dc was exactly 2x the true d
ratio = np.linalg.norm(dc_nacme_old) / np.linalg.norm(D_TRUE)
check('OLD nacme dc == 2 x true d', abs(ratio - 2.0) < 1e-9, f'ratio={ratio:.4f}')
ratio_n = np.linalg.norm(dc_num_old) / np.linalg.norm(D_TRUE)
check('OLD numerical_nac dc == 2 x true d', abs(ratio_n - 2.0) < 1e-9, f'ratio={ratio_n:.4f}')
# old numerical_nac h had opposite sign vs old nacme h (and vs standard)
nac_nacme_old = dc_nacme_old * gap_EjmEi(E)              # = 2 d (E_j-E_i) = +2 h_std
check('OLD numerical_nac h == -2 h_std (factor 2 + sign flip)',
      rel(nac_num_old, -2.0 * h_std) < 1e-9, f'rel={rel(nac_num_old, -2.0 * h_std):.2e}')
check('OLD nacme h == +2 h_std (factor 2, std sign)',
      rel(nac_nacme_old, 2.0 * h_std) < 1e-9, f'rel={rel(nac_nacme_old, 2.0 * h_std):.2e}')

# explicit demonstration of the two headline ratios from the handoff
print('\n# headline ratios (handoff: current/hand = 2.0, corrected/hand = 1.0):')
i, j = 0, 1
print(f'  pair {i+1}->{j+1}:  current(old)/hand = {dc_num_old[i,j]/hand_d[i,j]:.4f}   '
      f'corrected(new)/hand = {dc_num[i,j]/hand_d[i,j]:.4f}')

n_fail = sum(1 for _, ok in checks if not ok)
print('\n' + ('ALL CHECKS PASSED' if n_fail == 0 else f'{n_fail} CHECK(S) FAILED'))
raise SystemExit(1 if n_fail else 0)
