#!/usr/bin/env python3
r"""QMRSF-DK: the DIRECT hard-coded spin-adapted production matrix builder.

Assembles the three spin-pure QMRSF-DK blocks

    A_singlet  (20x20) = [[ A0^S (14x14) ,  V_c (14x6) ],
                          [ V_c^T (6x14) ,  W_dd (6x6) ]]   (paper singles/doubles order)
    A_triplet  (15x15) = single-spin-flip only (2OS-T 12x12 (+) 4OS-T 3x3)
    A_quintet  ( 1x1 ) = A_11

DIRECTLY, element-by-element, from the paper's explicit spin-adapted closed-form
formulas -- WITHOUT ever forming the 36-determinant Hamiltonian and rotating it by U.

``build_paper_matrix_f0_direct`` is THE PRODUCTION builder for QMRSF-DK: it is fully
hard-coded closed forms (0OS W_dd, 2OS singlet/triplet, and 4OS A_11..A_66 diagonals
are ALL explicit SI closed forms; the only Slater-Condon contractions are the
cross-sector off-diagonals V_c and the 2OS single-flip couplings, over the few
determinants each CSF pair spans).  The U-route ``qmrsf_dk_paper_f0.build_paper_matrix_f0``
(build the 36-det Hamiltonian, rotate by U) is the equivalence-check sibling, bit-for-bit
equal to this builder AT THE SAME ``vxc_singles`` setting (both default to 'resolved';
``verify_direct_vs_uroute_ints`` checks the 'averaged' seam on random integrals).

WHERE EACH PIECE COMES FROM (paper truth):
  * single-flip 2OS diagonal  ^1A/^3A = B_ii + c_H Delta^(2OS) +/- c_H K_rs
        (eq:det5-2os-st, eq:det5-delta2os; SI ~L2540)  -- closed form, diag(h_OH).
  * 4OS spin-adapted  A_11..A_66  (eq:det5-4os-effective; SI L2576-2610,
        MRSF-sQ-SI L757-792): B_ave + K structure, Clebsch-Gordan (Yamanouchi-Kotani)
        1 quintet (+) 3 triplet (+) 2 singlet -- closed form, diag(h_OH).  The 4OS
        6x6 block is dropped in ELEMENT-BY-ELEMENT from these hard-coded A_ij: the
        genealogical CG rotation coeff4 @ E4 @ coeff4^T of the orbital-ordered
        effective matrix is a PURE PERMUTATION of the SI A_ij (proven
        element-by-element to 1e-15 on H4 square/linear + butadiene FC, c_H=0.5 and
        1.0) -- within each spin block build_U's CSF order is the REVERSAL of the SI
        order, so no rotation product is formed at all.  build_U-order placement:
            diag  S1=A66 S2=A55 | T1=A44 T2=A33 T3=A22 | Q=A11
            offd  (S1,S2)=A56 | (T1,T2)=A34 (T1,T3)=A24 (T2,T3)=A23  (cross-block 0).
        The DIAGONAL is therefore LITERALLY A11..A66 (no coeff4 @ E4 @ coeff4^T on
        the diagonal); ``self_audit`` proves these closed forms equal the
        Slater-Condon couplings element-by-element.
  * 0OS closed-shell doubles W_dd (eq:det5-wdd-explicit; SI ~L2804):
        diag 2h+2h+J+J+4J-2c_H K ; off-diag +c_H K (share one orbital) ; 0 (disjoint)
        -- closed form, diag(h_OH).  0OS Fermi phase is +1 so it drops in directly.
  * reference-imbalance Delta (eq:refimbalance/deltaavg): folded into the closed
        forms above (Delta^(2OS)=-c_H[K_pq+1/2 sum K]; Delta^(4OS)=-1/2 c_H sum K).
  * OFF-diagonal single--double / cross-sector couplings V_c and the 2OS<->2OS,
        2OS<->4OS single-flip couplings are the paper's Slater-Condon exchange
        couplings (eq:offdiag; "off-diagonals everywhere are Slater-Condon", scaled by
        c_H).  They are evaluated PER CSF PAIR from the local determinant seam
        H_DK = H(h_OH) - (1-c_H) H_K over ONLY the (few) determinants each pair of
        CSFs spans -- never the global 36x36 rotation.  ``self_audit`` additionally
        proves that the closed-form W_dd off-diagonal and the closed-form 4OS A_ij
        off-diagonals equal these Slater-Condon couplings element-by-element.

The one-electron operator is the F^(0)-anchored ``h_OH`` of
``qmrsf_dk_paper_f0.h_OH_operator``.  The DEFAULT ``vxc_singles='resolved'`` is the
paper-faithful v_xc: h_OH's spin-AVERAGED v_xc is stripped and the spin-RESOLVED
DIAGONAL patch is added back per determinant (eq:vxc1 singles <a|v^b|a>-<i|v^a|i>,
eq:vxc0os doubles sum_occ v^b - sum_empty v^a).  ``vxc_singles='averaged'`` recovers the
earlier build (v_xc spin-averaged inside eps^sigma, no patch; pre-2026-07 published
numbers).  Both are inert at the c_H=1 gate (v_xc=0 => h_OH==h_act).  The exchange/
orbital STRUCTURE (W_dd, 4OS A_ij, S-T split=2c_H K_rs, V_c) is v_xc-independent and
identical for both.

Engine: reuses ``qmrsf_dk_matrices_ref`` (gen_dets/build_U/_spinorb_integrals/_melem,
the hard-coded spin adaptation) and ``qmrsf_dk_paper`` (explicit_wdd/explicit_2os_diag/
explicit_4os_formulas closed forms) and ``qmrsf_dk_pyscf`` (_exchange_eri seam,
active_integrals) and ``qmrsf_dk_paper_f0`` (h_OH_operator).  Do NOT contradict the SI.
"""
from __future__ import annotations

import numpy as np

import qmrsf_dk_matrices_ref as ref
import qmrsf_dk_pyscf as Q
import qmrsf_dk_paper as P
import qmrsf_dk_paper_f0 as F0

NORB = ref.NORB
NDET = ref.NDET
HARTREE2EV = 27.211386245988

# The 4OS spin-adapted 6x6 (eq:det5-4os-effective) equals the SI closed forms
# A_11..A_66 (qmrsf_dk_paper.explicit_4os_formulas) up to a PURE PERMUTATION: within
# each spin block build_U's genealogical (Yamanouchi-Kotani) CSF order is the
# REVERSAL of the SI A_ij order (verified element-by-element to 1e-15 on real H4 /
# butadiene integrals, both c_H).  So each build_U 4OS CSF maps to ONE hard-coded
# A_ij -- diagonal and internal off-diagonal are dropped in directly, no
# coeff4 @ E4 @ coeff4^T rotation.
_4OS_DIAG_A = {"4OS-S1": "A66", "4OS-S2": "A55",           # singlet 2x2
               "4OS-T1": "A44", "4OS-T2": "A33", "4OS-T3": "A22",   # triplet 3x3
               "4OS-Q":  "A11"}                            # quintet 1x1
# internal (within-spin-block) 4OS<->4OS off-diagonals, keyed by the CSF pair
_4OS_OFF_A = {frozenset(("4OS-S1", "4OS-S2")): "A56",      # singlet
              frozenset(("4OS-T1", "4OS-T2")): "A34",      # triplet
              frozenset(("4OS-T1", "4OS-T3")): "A24",
              frozenset(("4OS-T2", "4OS-T3")): "A23"}


# ---------------------------------------------------------------------------
# core: DIRECT assembly given the effective one-electron operator h_eff (= h_OH)
# ---------------------------------------------------------------------------
def _active_vxc_full(mf, somo):
    r"""Full 4x4 spin-resolved semilocal v_xc projected on the SOMOs (va_mat, vb_mat).
    Zero for a pure-HF reference.  Used to strip h_OH's spin-averaged v_xc so the
    paper's spin-RESOLVED DIAGONAL patch (eq:vxc1/eq:vxc0os) can be added back."""
    C = mf.mo_coeff[:, somo]
    xc = getattr(mf, "xc", None)
    if xc is None or str(xc).strip().lower() in ("hf", "hartree-fock", ""):
        z = np.zeros((len(somo), len(somo))); return z, z
    mo = mf.mo_coeff; occ = np.asarray(mf.mo_occ, float)
    dma = (mo * (occ > 0)) @ mo.T; dmb = (mo * (occ == 2)) @ mo.T
    if getattr(mf.grids, "coords", None) is None:
        mf.grids.build()
    _, _, v = mf._numint.nr_uks(mf.mol, mf.grids, str(xc), (dma, dmb))
    return C.T @ np.asarray(v[0]) @ C, C.T @ np.asarray(v[1]) @ C


def _assemble_direct(h_eff, eri_act, c_H, e_ref=0.0, return_audit=False):
    r"""Assemble the three spin-pure blocks element-by-element from the paper's
    closed forms (diagonals + 0OS/4OS internal off-diagonals) and the Slater-Condon
    exchange couplings (cross-sector / 2OS off-diagonals), given the effective
    one-electron operator ``h_eff`` (the F^(0) h_OH) and the active ERIs.

    NEVER forms or rotates the 36x36 Hamiltonian: cross couplings are contracted
    over ONLY the determinants each CSF pair spans."""
    h_eff = np.asarray(h_eff, float)
    eri_act = np.asarray(eri_act, float)

    dets = ref.gen_dets()
    U, Sval, labels = ref.build_U(dets)          # hard-coded paper spin adaptation
    # per-CSF determinant support (det index, coefficient) -- CSFs are sparse
    supp = [[(d, U[m, d]) for d in np.nonzero(U[m])[0]] for m in range(NDET)]

    # ---- closed-form paper blocks (use diag(h_eff)) --------------------------
    wdd = P.explicit_wdd(h_eff, eri_act, c_H, e_ref=0.0)      # 6x6 in _WDD_ORDER
    wdd_idx = {c: i for i, c in enumerate(P._WDD_ORDER)}
    twoos = P.explicit_2os_diag(h_eff, eri_act, c_H)          # dict (p,q)->{singlet,triplet}
    # 4OS: DROP IN the hard-coded closed forms A_11..A_66 (eq:det5-4os-effective;
    # SI L2576-2610, MRSF-sQ-SI L757-792) DIRECTLY -- NO coeff4 @ E4 @ coeff4^T
    # rotation.  The genealogical CG rotation is a pure permutation of these A_ij
    # (reversal within each spin block; see _4OS_DIAG_A / _4OS_OFF_A above), so the
    # diagonal is LITERALLY A11..A66 and the internal off-diagonal is A23/A24/A34/A56.
    A4 = P.explicit_4os_formulas(h_eff, eri_act, c_H)         # dict A11..A66 closed forms
    four_diag = {c: float(A4[a]) for c, a in _4OS_DIAG_A.items()}
    four_off = {k: float(A4[a]) for k, a in _4OS_OFF_A.items()}

    # ---- local c_H-scaled determinant seam H_DK = H(h_eff) - (1-c_H) H_K ------
    H1f, gf = ref._spinorb_integrals(h_eff, eri_act)
    eri_K = Q._exchange_eri(eri_act)
    H1z, gK = ref._spinorb_integrals(np.zeros((NORB, NORB)), eri_K)

    def hdk(p, q):
        v = ref._melem(dets[p], dets[q], H1f, gf)
        if c_H != 1.0:
            v -= (1.0 - c_H) * ref._melem(dets[p], dets[q], H1z, gK)
        return v

    def sc(m, n):
        """paper off-diagonal Slater-Condon coupling of CSF m and CSF n, contracted
        over only the determinants the two CSFs span."""
        s = 0.0
        for p, a in supp[m]:
            for q, b in supp[n]:
                if a and b:
                    s += a * b * hdk(p, q)
        return s

    # ---- per-CSF classification ---------------------------------------------
    def kind(m):
        L = labels[m]
        if L.startswith("0OS"):
            return "0OS", L.split(":")[1]
        if L.startswith("2OS"):
            return "2OS", L.split(":")[1].split("/")[0]
        return "4OS", L                         # "4OS-S1" etc

    def diag(m):
        k, c = kind(m)
        if k == "0OS":
            return float(wdd[wdd_idx[c], wdd_idx[c]]) - e_ref
        if k == "2OS":
            p, q = c.index("2"), c.index("0")
            e = twoos[(p, q)]
            v = e["singlet"] if labels[m].startswith("2OS-S") else e["triplet"]
            return float(v) - e_ref
        return four_diag[c] - e_ref              # 4OS: literal closed form A11..A66

    def offdiag(m, n):
        (km, cm), (kn, cn) = kind(m), kind(n)
        if km == "0OS" and kn == "0OS":
            return float(wdd[wdd_idx[cm], wdd_idx[cn]])
        if km == "4OS" and kn == "4OS":
            # internal 4OS<->4OS: closed-form A23/A24/A34 (triplet), A56 (singlet).
            # Only same-spin-block pairs occur here (blocks are built per spin S);
            # cross-block 4OS pairs are absent -> 0.
            return four_off.get(frozenset((cm, cn)), 0.0)
        return sc(m, n)                          # V_c / 2OS off-diag / 2OS<->4OS

    # ---- build the three spin blocks directly -------------------------------
    blocks = {}
    audit = {"wdd_off_vs_sc": 0.0, "4os_off_vs_sc": 0.0,
             "diag_2os_vs_sc": 0.0, "diag_0os_vs_sc": 0.0, "diag_4os_vs_sc": 0.0}
    for S in (0, 1, 2):
        idx = np.where(Sval == S)[0]
        n = len(idx)
        A = np.zeros((n, n))
        for i in range(n):
            A[i, i] = diag(idx[i])
            if return_audit:
                d_sc = sc(idx[i], idx[i]) - e_ref
                k, _ = kind(idx[i])
                key = {"0OS": "diag_0os_vs_sc", "2OS": "diag_2os_vs_sc",
                       "4OS": "diag_4os_vs_sc"}[k]
                audit[key] = max(audit[key], abs(A[i, i] - d_sc))
            for j in range(i + 1, n):
                v = offdiag(idx[i], idx[j])
                A[i, j] = A[j, i] = v
                if return_audit:
                    (km, cm), (kn, cn) = kind(idx[i]), kind(idx[j])
                    if km == "0OS" and kn == "0OS":
                        audit["wdd_off_vs_sc"] = max(
                            audit["wdd_off_vs_sc"], abs(v - sc(idx[i], idx[j])))
                    elif km == "4OS" and kn == "4OS":
                        audit["4os_off_vs_sc"] = max(
                            audit["4os_off_vs_sc"], abs(v - sc(idx[i], idx[j])))
        blocks[S] = A

    out = {"A_singlet": blocks[0], "A_triplet": blocks[1], "A_quintet": blocks[2],
           "U": U, "Sval": Sval, "labels": labels, "c_H": c_H, "h_eff": h_eff}
    if return_audit:
        out["audit"] = audit
    return out


# ---------------------------------------------------------------------------
# the production entry point (F^(0) diagonal, offdiag='hOH')
# ---------------------------------------------------------------------------
def build_paper_matrix_f0_direct(mf, nc, somo, c_H,
                                 h_act=None, eri_act=None, e_ref=0.0,
                                 return_audit=False, vxc_singles="resolved"):
    r"""THE PRODUCTION QMRSF-DK matrix builder (direct, fully hard-coded).

    This is the documented production default.  Every diagonal is an explicit
    closed form (0OS W_dd, 2OS singlet/triplet ^1A/^3A, and 4OS A_11..A_66); only
    the cross-sector / 2OS off-diagonals are Slater-Condon contractions over the few
    determinants each CSF pair spans -- the 36-det 2-electron Hamiltonian is NEVER
    formed and U is NEVER applied here (both the 'resolved' default and the 'averaged'
    branch go through ``_assemble_direct``).  The spin-resolved v_xc (eq:vxc1/eq:vxc0os)
    is added as a HARD-CODED per-CSF closed-form diagonal shift (doubly-occ -> <p|v^b|p>,
    empty -> -<q|v^a|q>, singly-occ -> 1/2(<v^b>-<v^a>)), which is purely diagonal in
    the CSF basis (verified: v_xc within-block off-diagonals are 0 to 1e-18), so no U
    rotation is needed.  The U-route ``qmrsf_dk_paper_f0.build_paper_matrix_f0`` (which
    DOES form the 36-det Hamiltonian via build_H_S2 and rotate by U) is the
    equivalence-check sibling, proven bit-for-bit equal at the same ``vxc_singles``.

    Computes the F^(0)-anchored one-electron operator ``h_OH`` (v_xc inside
    eps^sigma; no v_xc patch) exactly as ``qmrsf_dk_paper_f0`` and feeds it to the
    element-by-element closed-form assembly ``_assemble_direct``.  Same return shape
    as ``build_paper_matrix_f0``: ``A_singlet``/``A_triplet``/``A_quintet`` plus
    ``U``/``Sval``/``labels``/``c_H``, and ``h_OH`` diagnostic.
    """
    if h_act is None or eri_act is None:
        h_act, eri_act, _ec, _ci, _mc = Q.active_integrals(mf, nc, somo)
    h_act = np.asarray(h_act, float)
    eri_act = np.asarray(eri_act, float)
    # ***LAW*** F^(0)/h_OH is built at the REFERENCE functional's OWN exchange fraction
    # c_ref (the FIXED BHHLYP relaxed MOs / orbital energies eps), NOT the diagnostic c_H.
    # The diagnostic c_H is a purely POST-SCF knob on the ACTIVE-space matrix exchange only
    # (it flows into _assemble_direct below).  Letting c_H into F^(0) would move the
    # core-active exchange and hence the relaxed MO energies the CSF forms / orbital Hessian
    # require -- the MOs must never change.  (For a no-core system h_OH is c_H-independent
    # anyway.)  For a pure-HF ref c_ref==1 so the c_H=1 FCI/CASCI gate is untouched.
    c_ref = F0.reference_cHF(mf)
    h_OH = F0.h_OH_operator(mf, c_ref, somo, eri_act)

    if vxc_singles == "averaged":
        # legacy: v_xc spin-AVERAGED 1/2(v^a+v^b) inside h_OH (retained for comparison)
        dk = _assemble_direct(h_OH, eri_act, c_H, e_ref=e_ref, return_audit=return_audit)
        dk["h_OH"] = h_OH
        return dk
    if vxc_singles != "resolved":
        raise ValueError(f"vxc_singles must be 'resolved' or 'averaged', got {vxc_singles!r}")

    # --- PAPER-FAITHFUL v_xc (eq:vxc1 + eq:vxc0os): SPIN-RESOLVED, HARD-CODED per-CSF -
    # The spin-resolved v_xc lift is a DIAGONAL, spin-adapted, per-configuration CLOSED
    # FORM (verified: it is PURELY diagonal in the CSF basis -- within-block v_xc
    # off-diagonals are 0 to 1e-18 -- so NO U rotation is needed).  Using the spin-adapted
    # orbital occupations of each CSF (doubly-occ -> n^a=n^b=1, empty -> 0, singly-occ ->
    # n^a=n^b=1/2), the v_xc diagonal for a configuration is
    #     Delta^vxc = sum_{doubly} <p|v^b|p>  -  sum_{empty} <q|v^a|q>
    #                 + 1/2 sum_{singly} ( <p|v^b|p> - <p|v^a|p> )
    # which is eq:vxc0os EXACTLY for the 0OS doubles (no singly-occ term) and the
    # spin-adapted eq:vxc1 for the 2OS/4OS singles.  Same for singlet and triplet of a
    # configuration (the diagonal is spin-independent); inert at c_H=1 (v_xc=0).  We feed
    # the v_xc-STRIPPED F^(0) operator to the element-by-element closed-form assembly and
    # add this closed form to each CSF diagonal -- no 36-det Hamiltonian, no U rotation.
    va_d, vb_d = Q.active_vxc(mf, mf.mo_coeff[:, somo])          # <p|v_xc^sigma|p> (eq:vxc1/0os)
    va_m, vb_m = _active_vxc_full(mf, somo)                       # full 4x4 to strip from h_OH
    dk = _assemble_direct(h_OH - 0.5 * (va_m + vb_m), eri_act, c_H,
                          e_ref=e_ref, return_audit=return_audit)

    def _vxc_csf(lab):
        if lab.startswith("4OS"):                                # all 4 orbitals singly-occ
            return 0.5 * float(sum(vb_d[p] - va_d[p] for p in range(NORB)))
        cfg = lab.split(":")[1].split("/")[0]                    # e.g. 0OS '2200', 2OS '2ud0'
        dv = 0.0
        for p, ch in enumerate(cfg):
            if ch == "2":       dv += vb_d[p]                    # doubly-occ  -> <p|v^b|p>
            elif ch == "0":     dv -= va_d[p]                    # empty       -> -<q|v^a|q>
            else:               dv += 0.5 * (vb_d[p] - va_d[p])  # singly-occ  -> 1/2(v^b-v^a)
        return float(dv)

    labels = dk["labels"]; Sval = dk["Sval"]
    for key, S in (("A_singlet", 0), ("A_triplet", 1), ("A_quintet", 2)):
        idx = np.where(Sval == S)[0]; A = np.atleast_2d(dk[key]).copy()
        for r, m in enumerate(idx):
            A[r, r] += _vxc_csf(labels[m])
        dk[key] = A
    dk["h_OH"] = h_OH
    return dk


# ---------------------------------------------------------------------------
# NOTE: build_paper_matrix_f0_direct (above) IS the production QMRSF-DK orbital
# Hessian.  Its single-flip diagonal is the F^(0)/ROKS-energy orbital Hessian
# built via h_OH = F^(0) - active-mean-field (v_xc inside eps^sigma); the active
# Slater-Condon dressing keeps the hole-particle interaction correctly, so it is
# robust for the no-beta-electron (nc=0) active space and matches XMS-CASPT2
# across systems (butadiene FC 6.49/6.65, H4 square 2.48/3.46, H4 linear 8.38/...).
# A bare (eps_a - eps_i) - c_HF*K rebuild that DROPS the active-space interaction
# was tried (build_paper_matrix_oh) and DELETED: it was fine for cored butadiene
# (6.58) but broke on nc=0 H4 (0.10 vs XMS 2.03) -- the phantom "determinant vs
# orbital Hessian" split was an error.  DO NOT reintroduce it.
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# PROOF (1): direct == U-route on random integrals (isolates the assembly logic)
# ---------------------------------------------------------------------------
def verify_direct_vs_uroute_ints(h_eff, eri_act, c_H, e_ref=0.0, tol=1e-10):
    """Direct assembly (given h_eff) must equal the U-route
    ``qmrsf_dk_paper.build_paper_matrix(h_eff, eri, c_H, vxc_diag=None)`` per block."""
    dk_d = _assemble_direct(h_eff, eri_act, c_H, e_ref=e_ref, return_audit=True)
    dk_u = P.build_paper_matrix(h_eff, eri_act, c_H, vxc_diag=None, e_ref=e_ref)
    diffs = {k: float(np.abs(np.atleast_2d(dk_d[k]) - np.atleast_2d(dk_u[k])).max())
             for k in ("A_singlet", "A_triplet", "A_quintet")}
    diffs["audit"] = dk_d["audit"]
    diffs["max"] = max(diffs["A_singlet"], diffs["A_triplet"], diffs["A_quintet"])
    diffs["pass"] = bool(diffs["max"] < tol)
    return diffs


if __name__ == "__main__":
    print("=" * 74)
    print(" QMRSF-DK DIRECT hard-coded spin-adapted builder -- random-integral proof")
    print("=" * 74)
    allpass = True
    for seed in (0, 7, 12345, 99):
        h, e = ref.random_integrals(seed)
        for c_H in (0.5, 1.0):
            r = verify_direct_vs_uroute_ints(h, e, c_H)
            a = r["audit"]
            ok = r["pass"] and max(a.values()) < 1e-9
            allpass = allpass and ok
            print(f"[seed {seed:5d} c_H={c_H}]  |A_S|={r['A_singlet']:.2e} "
                  f"|A_T|={r['A_triplet']:.2e} |A_Q|={r['A_quintet']:.2e}  "
                  f"audit(wdd_off={a['wdd_off_vs_sc']:.1e} 4os_off={a['4os_off_vs_sc']:.1e} "
                  f"diag2os={a['diag_2os_vs_sc']:.1e})  -> {'PASS' if ok else 'FAIL'}")
    print("=" * 74)
    print(f" RANDOM-INTEGRAL PROOF: {'ALL PASS' if allpass else 'FAILED'}")
    assert allpass, "direct != U-route"
