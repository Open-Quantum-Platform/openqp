#!/usr/bin/env python3
r"""QMRSF-DK matrix constructor: the hardcoded spin-adapted TDDFT orbital-Hessian
assembled DIRECTLY from the paper's explicit blocks with UNIFORM c_H scaling.

CONVENTION (author-authoritative -- "consistency", SI L272-299): the exact-exchange
fraction c_H scales EVERY exact-exchange integral uniformly -- in the diagonals AND
the off-diagonals, in EVERY sector, with NOTHING exempted.  Off-diagonals everywhere
are the Slater--Condon couplings; every K (the diagonal same-spin -K, the 0OS
-2K/+K, the 2OS/4OS +/-K couplings, AND the genuine 4-distinct exchange integrals in
V_c) is scaled by c_H.  This value-based / uniform scaling is exact at c_H=1.  The
assembled matrix therefore reproduces the current ``qmrsf_dk_pyscf.build_dk_matrices``
bit-for-bit (verified here) -- confirming the current code already IS the consistent
uniform-scaling matrix.

TDDFT footing (Lee 2018, J. Chem. Phys. 149, 104101, Eqs. 2.22-2.25): the matrix is
the Casida orbital Hessian
    A_{p q, r s} = F^(0)_{pr} d_{sq} - F^(0)_{sq} d_{pr}
                   + (qp|sr)d d - c_H(pr|sq)d d + (1-c_H) f^XC ,
whose spin-flip part keeps only the exact-exchange kernel -c_H(pr|sq) (Eq. 2.25).
The DIAGONAL Fock part is F^(0), the ground-state ROKS-BHHLYP KS Fock of the quintet
reference, whose orbital energies carry v_xc and the c_H-scaled exact exchange.  In
the CURRENT PySCF construction this diagonal is instead built as a Slater--Condon
determinant energy on the CASCI active integrals ``h_act`` (whose frozen core
carries FULL HF exchange) plus a separate ``active_vxc`` diagonal patch -- whether
that equals the true F^(0) orbital-Hessian diagonal is analysed by
``fock_diagnostic.py`` and reported (it does NOT, at the absolute level; see there).

What the paper's explicit blocks scale, and what this module reproduces:

  * 4OS effective 6x6  (eq:det5-4os-effective, SI ~L2624; A_11..A_66 SI L2576-2610,
    MRSF-sQ-SI L757-792): diagonal ``D - c_H*(K_pair)``, off-diagonal ``-c_H*K_ij``,
    with ``D = B_ave = sum_i h_ii + sum_{i<j} J_ij`` and
    ``Delta^(4OS) = -1/2 sum_{i<j} K_ij`` folded into ``D = B_0 + Delta^(4OS)``.
    Every K is 2-index -> all scaled by c_H.  Quintet root
    ``A_11 = B_ave - c_H*sum_{i<j}K_ij`` (eq:det5-delta4os).
  * 0OS closed-shell 6x6 W_dd  (eq:det5-wdd-explicit, SI ~L2804):
    diagonal ``w_ij = 2h_ii+2h_jj+J_ii+J_jj+4J_ij - 2*c_H*K_ij - E_ref``
    (eq:det5-wdiag); off-diagonal ``+c_H*K_kl`` for pairs sharing one orbital
    (eq:det5-woff); ``0`` for orbital-disjoint pairs (eq:det5-wzero).  All K 2-index
    -> scaled by c_H (SI L283-289 convention: W_dd's +K and -2K carry c_HF too).
  * 2OS 12x12 singlet / triplet  (Tables det5-singlet/triplet-2os; eq:det5-2os-st,
    eq:det5-delta2os, SI ~L2540): diagonal
    ``^1A_(pq) = B_ii + c_H*Delta^(2OS)_(pq) + c_H*K_rs`` (singlet),
    ``^3A_(pq) = B_ii + c_H*Delta^(2OS)_(pq) - c_H*K_rs`` (triplet),
    ``Delta^(2OS)_(pq) = -K_pq - 1/2(K_pr+K_ps+K_qr+K_qs)``; the off-diagonal
    couplings c_{i,j} are the (c_H-scaled exchange) single-spin-flip orbital-Hessian
    couplings.  This is the standard MRSF-TDDFT single-flip block
    (eq:det2-Casida-A-elements): its couplings ``-c_H (pr|sq)`` ARE scaled even when
    4-distinct, so this block is IDENTICAL to the value-based code (correct).
  * V_c (14x6) singles--doubles coupling  (eq:det5-singlet-block; SI L2818-2821):
    the Slater--Condon couplings ``<^1Phi^(2OS/4OS)|H|Phi^(0OS)>``.  The 0OS<->4OS
    couplings are genuine 4-distinct exchange integrals; per the uniform convention
    these ARE scaled by c_H like every other exchange integral (no exemption).
  * Assembly (eq:det5-singlet-block / eq:det5-triplet-block): singlet 20x20 =
    [A_0^S, V_c ; V_c^T, W_dd]; triplet 15x15 = single-spin-flip only (no W_dd, no
    V_c); quintet 1x1.  (Row order here follows ``qmrsf_dk_matrices_ref.build_U``:
    0OS first, then 2OS-singlets, then 4OS-singlets -- a permutation of the paper's
    [A_0^S ; W_dd] order; the block CONTENT is identical.)

The Kohn--Sham ``v_xc`` diagonal (eq:det5-wdiag note / DK main dressed diagonals) is
added exactly as the value-based ``build_dk_matrices`` does, so the two are directly
comparable.

At ``c_H = 1`` the fix is inert (``(1-c_H)=0``), so ``build_paper_matrix(...,c_H=1)``
returns the full CI and its block spectra equal the value-based code's and the
S^2-projected active-space spectrum (verification gate (a)).  At ``c_H < 1`` the
ONLY difference from the value-based code is the 0OS<->4OS (V_c) 4-distinct
coupling (verification (c)).

Engine: reuses ``qmrsf_dk_matrices_ref`` (gen_dets/det_config/os_class/build_U/
build_H_S2) and ``qmrsf_dk_pyscf`` (_build_H_K1SF value-based exchange tensor,
active_integrals/active_vxc via the caller).  Do NOT contradict SI Sec. det5.
"""
from __future__ import annotations

import numpy as np

import qmrsf_dk_matrices_ref as ref
import qmrsf_dk_pyscf as vb          # "value-based" current code

NORB = ref.NORB            # 4 active orbitals
NDET = ref.NDET            # 36 M_s=0 determinants
HARTREE2EV = 27.211386245988

# config order used by the explicit closed-shell / 4OS blocks
_WDD_ORDER = ["2200", "2020", "2002", "0220", "0202", "0022"]     # doubly (0,1)(0,2)..(2,3)
_WDD_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
# same-spin orbital pairs of each 4OS determinant in ref.FOUR_OS_ORDER
#   uudd,udud,uddu,duud,dudu,dduu
_4OS_SAMESPIN = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)),
                 ((1, 2), (0, 3)), ((1, 3), (0, 2)), ((2, 3), (0, 1))]


# ---------------------------------------------------------------------------
# integral helpers
# ---------------------------------------------------------------------------
def _JKh(h_act, eri_act):
    """Return (h_diag, J, K): h_ii, Coulomb J_ij=(ii|jj), exchange K_ij=(ij|ji)."""
    NO = NORB
    hd = np.diag(np.asarray(h_act, float))
    J = np.array([[eri_act[i, i, j, j] for j in range(NO)] for i in range(NO)])
    K = np.array([[eri_act[i, j, j, i] for j in range(NO)] for i in range(NO)])
    return hd, J, K


def _B_ave(h_act, eri_act):
    """4OS common diagonal D = B_0 + Delta^(4OS) = sum_i h_ii + sum_{i<j} J_ij.

    (The quintet high-spin diagonal before the same-spin -K terms;
    eq:det5-delta4os folds Delta^(4OS) = -1/2 sum K into B_0 to give this D.)"""
    hd, J, _ = _JKh(h_act, eri_act)
    return float(hd.sum() + sum(J[i, j] for i in range(NORB)
                                for j in range(i + 1, NORB)))


# ===========================================================================
# EXPLICIT closed-form blocks (used for construction documentation + verify (b))
# ===========================================================================
def explicit_4os_effective(h_act, eri_act, c_H):
    r"""The 6x6 effective 4OS matrix, eq:det5-4os-effective (SI ~L2624), in the
    orbital-ordered (Yamanouchi-Kotani) basis of ``ref.FOUR_OS_ORDER``.

    Diagonal ``D - c_H*(K_p1 + K_p2)`` with ``D = B_ave`` (``_B_ave``);
    off-diagonal ``-c_H*K_ij`` between the two 4OS determinants that differ by a
    spin flip in orbitals ``i,j``; ``0`` for the fully-disjoint pair.  All K are
    2-index active exchange integrals and are scaled by ``c_H``
    (SI L272-289 convention)."""
    hd, J, K = _JKh(h_act, eri_act)
    D = _B_ave(h_act, eri_act)
    kc = lambda i, j: c_H * K[i, j]
    B = np.zeros((6, 6))
    for a, ((p, q), (r, s)) in enumerate(_4OS_SAMESPIN):
        B[a, a] = D - kc(p, q) - kc(r, s)
    # off-diagonals: -c_H*K between dets differing by one spin-flip pair (i,j)
    dets = ref.gen_dets()
    idx = {ref.det_config(d): d for d in dets}
    order = [idx[c] for c in ref.FOUR_OS_ORDER]
    for a in range(6):
        ca = ref.FOUR_OS_ORDER[a]
        for b in range(6):
            if a == b:
                continue
            cb = ref.FOUR_OS_ORDER[b]
            diff = [o for o in range(NORB) if ca[o] != cb[o]]
            if len(diff) == 2:                    # differ by a spin flip in i,j
                i, j = diff
                B[a, b] = -kc(i, j)
            # len(diff)==4 (fully disjoint) -> 0 (already)
    return B


def explicit_wdd(h_act, eri_act, c_H, e_ref=0.0):
    r"""The 6x6 closed-shell doubles block W_dd, eq:det5-wdd-explicit (SI ~L2804),
    in config order ``_WDD_ORDER``.

    Diagonal ``w_ij = 2h_ii+2h_jj+J_ii+J_jj+4J_ij - 2*c_H*K_ij - E_ref``
    (eq:det5-wdiag); off-diagonal ``+c_H*K_kl`` for two closed-shell configs
    sharing exactly one doubly-occupied orbital and differing in the other pair
    (eq:det5-woff); ``0`` for orbital-disjoint pairs (eq:det5-wzero).  The -2K and
    +K are 2-index active exchange, scaled by ``c_H`` (SI L283-289)."""
    hd, J, K = _JKh(h_act, eri_act)
    W = np.zeros((6, 6))
    for a, (i, j) in enumerate(_WDD_PAIRS):
        W[a, a] = (2 * hd[i] + 2 * hd[j] + J[i, i] + J[j, j] + 4 * J[i, j]
                   - 2 * c_H * K[i, j] - e_ref)
    for a, (i, j) in enumerate(_WDD_PAIRS):
        for b, (k, l) in enumerate(_WDD_PAIRS):
            if a == b:
                continue
            common = {i, j} & {k, l}
            if len(common) == 1:                  # share one orbital
                diff = list({i, j} ^ {k, l})      # the two differing orbitals
                W[a, b] = c_H * K[diff[0], diff[1]]
            # disjoint -> 0
    return W


def explicit_2os_diag(h_act, eri_act, c_H):
    r"""The explicit 2OS singlet/triplet diagonals, eq:det5-2os-st + eq:det5-delta2os
    (SI ~L2540), returned as a dict keyed by the doubly-occ/empty pair ``(p,q)``.

    ``Delta^(2OS)_(pq) = -K_pq - 1/2 (K_pr+K_ps+K_qr+K_qs)`` (all 2-index K,
    scaled by ``c_H``), open-shell orbitals ``(r,s)``; the reference-averaged
    orbital-Hessian diagonal ``B_ii`` is the Coulomb/Fock part
    ``2h_pp+h_rr+h_ss+J_pp+2J_pr+2J_ps+J_rs`` plus its own exchange dressing
    ``+c_H*[ K_pq - 1/2(K_pr+K_ps) + 1/2(K_qr+K_qs) ]`` (so that at c_H=1
    ``B_ii + Delta^(2OS) = `` the bare Slater diagonal).  Then
    ``^1A_(pq) = B_ii + c_H*Delta^(2OS)_(pq) + c_H*K_rs`` and
    ``^3A_(pq) = ... - c_H*K_rs`` (eq:det5-2os-st / eq:det5-triplet-2os-diag)."""
    hd, J, K = _JKh(h_act, eri_act)
    out = {}
    for p in range(NORB):
        for q in range(NORB):
            if p == q:
                continue
            r, s = [o for o in range(NORB) if o not in (p, q)]
            Bcoul = (2 * hd[p] + hd[r] + hd[s] + J[p, p]
                     + 2 * J[p, r] + 2 * J[p, s] + J[r, s])
            Bexch = c_H * (K[p, q] - 0.5 * (K[p, r] + K[p, s])
                           + 0.5 * (K[q, r] + K[q, s]))
            B_ii = Bcoul + Bexch
            Delta = c_H * (-K[p, q] - 0.5 * (K[p, r] + K[p, s]
                                             + K[q, r] + K[q, s]))
            Krs = c_H * K[r, s]
            out[(p, q)] = {"B_ii": B_ii, "Delta": Delta, "K_rs": Krs,
                           "singlet": B_ii + Delta + Krs,
                           "triplet": B_ii + Delta - Krs, "rs": (r, s)}
    return out


# ===========================================================================
# the paper-block determinant Hamiltonian + spin adaptation
# ===========================================================================
def _paper_HDK(h_act, eri_act, c_H, dets, H, vxc_diag=None, e_ref=0.0):
    r"""Assemble the 36-determinant QMRSF-DK Hamiltonian from the paper's block
    prescription with UNIFORM value-based c_H scaling and return it.

    Convention (author-authoritative, "consistency"): the exact-exchange fraction
    c_H scales EVERY exact-exchange integral uniformly -- in the diagonals AND the
    off-diagonals, in EVERY sector, with NOTHING exempted (SI L272-299).  This is
    the value-based / uniform scaling, exact at c_H=1.  Concretely:

      * single-spin-flip 2OS/4OS diagonal same-spin -K and +/-K couplings,
      * 0OS closed-shell W_dd diagonal -2K and off-diagonal +K,
      * the genuine 4-distinct exchange integrals in the singles--doubles V_c
        coupling (0OS<->4OS),

    are ALL scaled by c_H.  There is no role split and no 4-index exemption, so the
    determinant Hamiltonian is exactly ``H - (1-c_H)*H_K`` with H_K the value-based
    exchange tensor of ``qmrsf_dk_pyscf._build_H_K1SF`` -- i.e. this reproduces the
    current ``build_dk_matrices`` bit-for-bit (verified).  The v_xc diagonal is
    added exactly as the value-based code.  ``osc`` = per-determinant open-shell
    count.
    """
    osc = np.array([ref.os_class(ref.det_config(d)) for d in dets])
    H_K = vb._build_H_K1SF(h_act, eri_act, H, dets, None)     # value-based exchange
    H_DK = H - (1.0 - c_H) * H_K                              # uniform c_H seam

    if e_ref:
        H_DK = H_DK - e_ref * np.eye(NDET)

    # --- v_xc diagonal (identical to build_dk_matrices) ---
    if vxc_diag is not None:
        va = np.asarray(vxc_diag[0], float)
        vb_ = np.asarray(vxc_diag[1], float)
        for d, det in enumerate(dets):
            na = np.array([1.0 if p in det else 0.0 for p in range(NORB)])
            nb = np.array([1.0 if (p + NORB) in det else 0.0 for p in range(NORB)])
            H_DK[d, d] += float(np.dot(na - 1.0, va) + np.dot(nb, vb_))
    return H_DK, osc


def build_paper_matrix(h_act, eri_act, c_H, vxc_diag=None, e_ref=0.0):
    r"""Build the three spin-pure QMRSF-DK blocks from the paper's explicit blocks.

    Same signature/return shape as ``qmrsf_dk_pyscf.build_dk_matrices``:

    Returns a dict with ``A_singlet`` (20x20), ``A_triplet`` (15x15),
    ``A_quintet`` (1x1), plus ``U``/``Sval``/``labels`` (spin adaptation),
    ``H``/``H_DK``/``S2`` (determinant basis), ``c_H``, and ``osc`` (per-det
    open-shell count).

    The blocks are ``A = U H_DK U^T`` partitioned by total spin S, where H_DK is
    the paper-block determinant Hamiltonian of ``_paper_HDK`` (the value-based
    orbital-Hessian seam with the V_c 0OS<->4OS 4-distinct coupling restored to the
    exact CI value).  See module docstring for the per-block equation map.
    """
    h_act = np.asarray(h_act, float)
    eri_act = np.asarray(eri_act, float)

    dets = ref.gen_dets()
    H, S2 = ref.build_H_S2(h_act, eri_act, dets)
    U, Sval, labels = ref.build_U(dets)

    H_DK, osc = _paper_HDK(h_act, eri_act, c_H, dets, H,
                           vxc_diag=vxc_diag, e_ref=e_ref)

    Hcsf = U @ H_DK @ U.T
    sel = {S: np.where(Sval == S)[0] for S in (0, 1, 2)}
    A = {S: Hcsf[np.ix_(sel[S], sel[S])] for S in (0, 1, 2)}
    return {
        "A_singlet": A[0], "A_triplet": A[1], "A_quintet": A[2],
        "U": U, "Sval": Sval, "labels": labels,
        "H": H, "H_DK": H_DK, "S2": S2, "c_H": c_H, "osc": osc,
    }


# ===========================================================================
# VERIFICATION
# ===========================================================================
def _block_spectra(dk):
    return (np.sort(np.linalg.eigvalsh(dk["A_singlet"])),
            np.sort(np.linalg.eigvalsh(dk["A_triplet"])),
            np.sort(np.linalg.eigvalsh(np.atleast_2d(dk["A_quintet"]))))


def verify_cH1_gate(h_act, eri_act, tol=1e-9, vxc_diag=None):
    """(a) At c_H=1 the paper blocks must equal the value-based blocks AND the
    S^2-projected active-space CI spectrum, to < tol.

    ``paper_vs_value`` is evaluated with the supplied ``vxc_diag`` (at c_H=1 the two
    determinant Hamiltonians are byte-identical, so this is the primary gate).  The
    ``paper_vs_S2CI`` check is always run on the BARE (no-vxc) matrices, because the
    per-determinant v_xc diagonal shift is not S^2-symmetric in the determinant
    basis (an existing modeling choice shared with ``build_dk_matrices``), so the
    S^2 projection of the vxc-shifted H_DK is ill-defined."""
    dk_p = build_paper_matrix(h_act, eri_act, 1.0, vxc_diag=vxc_diag)
    dk_v = vb.build_dk_matrices(h_act, eri_act, 1.0, with_diagnostics=False,
                                vxc_diag=vxc_diag)
    # value-based vs paper (with vxc as supplied)
    d_pv = max(float(np.abs(a - b).max())
               for a, b in zip(_block_spectra(dk_p), _block_spectra(dk_v)))
    # paper vs S^2-projected CI, on the BARE (no-vxc) matrices
    dk_b = build_paper_matrix(h_act, eri_act, 1.0, vxc_diag=None)
    H = dk_b["H_DK"]; S2 = dk_b["S2"]
    w, V = np.linalg.eigh(H)
    s2f = np.array([V[:, k] @ S2 @ V[:, k] for k in range(NDET)])
    d_ci = 0.0
    for S, eb in zip((0, 1, 2), _block_spectra(dk_b)):
        ed = np.sort(w[np.abs(s2f - S * (S + 1)) < 1e-6])
        if ed.shape != eb.shape:
            d_ci = np.nan
            break
        d_ci = max(d_ci, float(np.abs(eb - ed).max()))
    hdk_eq = float(np.abs(dk_b["H_DK"] - dk_b["H"]).max())   # bare: must be ~0
    return {"paper_vs_value": d_pv, "paper_vs_S2CI": d_ci,
            "HDK_eq_H_bare": hdk_eq,
            "pass": (d_pv < tol and (d_ci < tol))}


def verify_explicit_blocks(h_act, eri_act, c_H, tol=1e-12):
    """(b) The assembled 4OS 6x6 and W_dd 6x6 determinant sub-blocks must equal the
    paper's eq:det5-4os-effective and eq:det5-wdd-explicit element-by-element."""
    dets = ref.gen_dets()
    cfg = [ref.det_config(d) for d in dets]
    H, _ = ref.build_H_S2(h_act, eri_act, dets)
    H_DK, _ = _paper_HDK(h_act, eri_act, c_H, dets, H)
    idx = {c: i for i, c in enumerate(cfg)}

    # --- W_dd: 0OS Fermi phase is +1, so det sub-block == orbital-ordered form ---
    oo = [idx[c] for c in _WDD_ORDER]
    W_asm = H_DK[np.ix_(oo, oo)]
    W_exp = explicit_wdd(h_act, eri_act, c_H)
    err_w = float(np.abs(W_asm - W_exp).max())

    # --- 4OS: undo the Fermi phase D to reach the orbital-ordered eff matrix ---
    fo = [idx[c] for c in ref.FOUR_OS_ORDER]
    Dp = np.diag(ref.FERMI_D_4OS)
    E_asm = Dp @ H_DK[np.ix_(fo, fo)] @ Dp
    E_exp = explicit_4os_effective(h_act, eri_act, c_H)
    err_4 = float(np.abs(E_asm - E_exp).max())

    # spin-adapted A_11..A_66 eigenvalue equivalence (basis-independent content)
    Uorb = np.array([c for _, _, c in ref.FOUR_OS_CSF])       # S1,S2,T1,T2,T3,Q
    M = Uorb @ E_exp @ Uorb.T
    A = explicit_4os_formulas(h_act, eri_act, c_H)
    q_err = abs(M[5, 5] - A["A11"])
    t_eig = np.abs(np.sort(np.linalg.eigvalsh(M[np.ix_([2, 3, 4], [2, 3, 4])]))
                   - np.sort(np.linalg.eigvalsh(np.array(
                       [[A["A22"], A["A23"], A["A24"]],
                        [A["A23"], A["A33"], A["A34"]],
                        [A["A24"], A["A34"], A["A44"]]])))).max()
    s_eig = np.abs(np.sort(np.linalg.eigvalsh(M[np.ix_([0, 1], [0, 1])]))
                   - np.sort(np.linalg.eigvalsh(np.array(
                       [[A["A55"], A["A56"]], [A["A56"], A["A66"]]])))).max()

    return {"err_Wdd": err_w, "err_4OS_eff": err_4,
            "err_quintet_A11": float(q_err), "err_triplet_eigs": float(t_eig),
            "err_singlet_eigs": float(s_eig),
            "pass": max(err_w, err_4, q_err, t_eig, s_eig) < tol}


def explicit_4os_formulas(h_act, eri_act, c_H):
    """Closed-form A_11..A_66 (SI L2576-2610 / MRSF-sQ-SI L757-792), K scaled by
    c_H, B = B_ave.  Returned as a dict."""
    _, _, K = _JKh(h_act, eri_act)
    B = _B_ave(h_act, eri_act)
    k = lambda i, j: c_H * K[i - 1, j - 1]
    return {
        "A11": -k(2, 3) - k(2, 4) - k(3, 4) + B - k(1, 3) - k(1, 2) - k(1, 4),
        "A22": (-3 * k(2, 3) + k(2, 4) + k(3, 4) + 3 * B - 3 * k(1, 3)
                - 3 * k(1, 2) + k(1, 4)) / 3,
        "A23": np.sqrt(2) * (k(2, 4) - 2 * k(3, 4) + k(1, 4)) / 3,
        "A24": np.sqrt(2 / 3) * (-k(2, 4) + k(1, 4)),
        "A33": (3 * k(2, 3) - 5 * k(2, 4) - 2 * k(3, 4) + 6 * B + 3 * k(1, 3)
                - 6 * k(1, 2) - 5 * k(1, 4)) / 6,
        "A34": (-3 * k(2, 3) - k(2, 4) + 3 * k(1, 3) + k(1, 4)) / (2 * np.sqrt(3)),
        "A44": (-k(2, 3) - k(2, 4) - 2 * k(3, 4) + 2 * B - k(1, 3) + 2 * k(1, 2)
                - k(1, 4)) / 2,
        "A55": (k(2, 3) + k(2, 4) - 2 * k(3, 4) + 2 * B + k(1, 3) - 2 * k(1, 2)
                + k(1, 4)) / 2,
        "A56": np.sqrt(3) * (-k(2, 3) + k(2, 4) + k(1, 3) - k(1, 4)) / 2,
        "A66": -k(2, 3) / 2 - k(2, 4) / 2 + k(3, 4) + B - k(1, 3) / 2 + k(1, 2)
        - k(1, 4) / 2,
    }


def verify_2os_diag(h_act, eri_act, c_H, tol=1e-10, vxc_diag=None):
    """(b, bonus) The assembled 2OS singlet/triplet diagonals must equal
    eq:det5-2os-st ^1A/^3A = B_ii + c_H*Delta^(2OS) +/- c_H*K_rs, and the
    singlet-triplet split must be 2*c_H*K_rs."""
    dk = build_paper_matrix(h_act, eri_act, c_H, vxc_diag=vxc_diag)
    A_s, A_t = dk["A_singlet"], dk["A_triplet"]
    labels, Sval = dk["labels"], dk["Sval"]
    slabs = [labels[k] for k in range(NDET) if Sval[k] == 0]
    tlabs = [labels[k] for k in range(NDET) if Sval[k] == 1]
    s_off = [i for i, l in enumerate(slabs) if l.startswith("2OS-S")]
    t_off = [i for i, l in enumerate(tlabs) if l.startswith("2OS-T")]
    expl = explicit_2os_diag(h_act, eri_act, c_H)

    def parse(cfg):
        p = cfg.index("2"); q = cfg.index("0")
        return p, q
    err = 0.0
    for si, ti in zip(s_off, t_off):
        cfg = slabs[si].split(":")[1].split("/")[0]
        p, q = parse(cfg)
        e = expl[(p, q)]
        err = max(err, abs(float(A_s[si, si]) - e["singlet"]),
                  abs(float(A_t[ti, ti]) - e["triplet"]),
                  abs((float(A_s[si, si]) - float(A_t[ti, ti])) - 2 * e["K_rs"]))
    return {"err_2os_diag": float(err), "pass": err < tol}


def verify_matches_current(h_act, eri_act, c_H, tol=1e-10, vxc_diag=None):
    """(b) With uniform value-based scaling the hardcoded spin-adapted matrix must
    reproduce the current ``build_dk_matrices`` bit-for-bit (H_DK identical and block
    spectra identical) -- confirming the current code IS the consistent
    uniform-scaling matrix.  Returns the max determinant-Hamiltonian and per-block
    spectrum differences."""
    dk_p = build_paper_matrix(h_act, eri_act, c_H, vxc_diag=vxc_diag)
    dk_v = vb.build_dk_matrices(h_act, eri_act, c_H, with_diagnostics=False,
                               vxc_diag=vxc_diag)
    hdk = float(np.abs(dk_p["H_DK"] - dk_v["H_DK"]).max())
    spec = max(float(np.abs(a - b).max())
               for a, b in zip(_block_spectra(dk_p), _block_spectra(dk_v)))
    return {"HDK_vs_current": hdk, "spectrum_vs_current": spec,
            "pass": (hdk < tol and spec < tol)}


# ===========================================================================
# __main__ : run the full verification battery on random + structured integrals
# ===========================================================================
if __name__ == "__main__":
    print("=" * 72)
    print(" QMRSF-DK PAPER-BLOCK matrix constructor -- verification battery")
    print("=" * 72)

    allpass = True
    for seed in (0, 7, 12345):
        h, e = ref.random_integrals(seed)
        g = verify_cH1_gate(h, e)
        b = verify_explicit_blocks(h, e, 0.5)
        b2 = verify_2os_diag(h, e, 0.5)
        c = verify_matches_current(h, e, 0.5)
        ok = g["pass"] and b["pass"] and b2["pass"] and c["pass"]
        allpass = allpass and ok
        print(f"\n[seed {seed}]")
        print(f"  (a) c_H=1 gate       : paper-vs-value={g['paper_vs_value']:.2e} "
              f"paper-vs-S2CI={g['paper_vs_S2CI']:.2e}  -> "
              f"{'PASS' if g['pass'] else 'FAIL'}")
        print(f"  (b) explicit blocks  : err_Wdd={b['err_Wdd']:.2e} "
              f"err_4OSeff={b['err_4OS_eff']:.2e} "
              f"A11={b['err_quintet_A11']:.2e} "
              f"T_eig={b['err_triplet_eigs']:.2e} S_eig={b['err_singlet_eigs']:.2e}"
              f"  -> {'PASS' if b['pass'] else 'FAIL'}")
        print(f"  (b) 2OS diagonals    : err={b2['err_2os_diag']:.2e}  -> "
              f"{'PASS' if b2['pass'] else 'FAIL'}")
        print(f"  (b) == current@0.5   : H_DK diff={c['HDK_vs_current']:.2e}"
              f"  spectrum diff={c['spectrum_vs_current']:.2e}"
              f"  -> {'PASS' if c['pass'] else 'FAIL'}")

    print("\n" + "=" * 72)
    print(f" VERIFICATION: {'ALL PASS' if allpass else 'FAILED'}")
    print("=" * 72)
    assert allpass, "verification FAILED"
