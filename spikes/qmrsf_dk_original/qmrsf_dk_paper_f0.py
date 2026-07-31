#!/usr/bin/env python3
r"""QMRSF-DK paper matrix with the TRUE F^(0) single-flip diagonal.

Production sibling of ``qmrsf_dk_paper.build_paper_matrix``.  It replaces the
value-based single-flip diagonal (a Slater--Condon determinant energy on the CASCI
``h_act`` -- whose frozen core carries FULL unscaled HF exchange and no v_xc -- plus
a per-determinant ``active_vxc`` diagonal patch) by the ground-state ROKS-BHHLYP KS
Fock F^(0) orbital-Hessian diagonal (Lee 2018, J. Chem. Phys. 149, 104101,
Eqs. 2.22/2.25; DK main text eq:dressed1sf)

    [A_0]_nn = (eps^beta_a - eps^alpha_i) - c_HF (ii|aa).

NOTE ON ``(ii|aa)``:  in this orbital-Hessian diagonal ``-c_HF (ii|aa)`` (chemist /
Mulliken notation) is the EXACT-EXCHANGE KERNEL contribution -- the i=j, a=b instance
of the Casida spin-flip exchange term ``-c_HF (ij|ab)`` (Lee 2018 Eq. 2.25).  It is
NOT a Coulomb/Hartree integral (the Hartree coupling ``(ia|jb)`` vanishes for the
spin-flip A-matrix), even though by integral VALUE ``(ii|aa)`` coincides with a
Coulomb self-integral.  It must NOT be "reconciled with" or replaced by ``(ij|ji)=K``:
that is a DIFFERENT integral -- the spectator same-spin exchange that gives the 2OS
S-T split ``+/- c_HF K_rs`` (eq:det5-2os-st).

HOW.  Build the F^(0)-anchored effective one-electron operator (exactly the
``chf_check/f0_diagonal_check.py`` machinery)

    h_OH = 1/2 (F^a + F^b)_active - (active mean field),
    F^(0),sigma = h + J[rho] - c_H K[rho^sigma] + v_xc^sigma[rho]  (ROKS quintet),
    (active mean field)^a_pq = sum_{k in act}[(pq|kk) - c_H (pk|kq)],
    (active mean field)^b_pq = sum_{k in act} (pq|kk),

and feed ``h_OH`` to the SAME Slater--Condon spin-adaptation machinery.  Its
determinant diagonal then collapses to the orbital-energy difference above; v_xc is
now inside eps^sigma, so the explicit ``active_vxc`` patch is dropped.

*** REFERENCE-CONSISTENCY (the load-bearing subtlety) ***

    delta_h = h_OH - h_act = (1 - c_H) K_core + v_xc

is DOMINATED by a large, near-orbital-INDEPENDENT constant (~ -5 eV at BHHLYP: the
frozen-core exchange that is unscaled in h_act but c_H-scaled in F^(0)); only the
~0.4 eV orbital-dependent variation (and v_xc) is physical.  That constant cancels
in EVERY excitation energy ONLY IF the F^(0) diagonal is applied CONSISTENTLY to all
36 determinants (0OS closed-shell doubles INCLUDED), because excitations are taken
relative to the quintet root and the 0OS<->singles blocks are coupled by V_c.  A
"single-flip diagonal only" surgery that leaves the 0OS ground state on the h_act
scale does NOT cancel the constant and produces a ~5 eV artifact (butadiene S1 -> 1.5
eV).  Hence this module shifts the diagonal of the WHOLE 36-det Hamiltonian (Route A
below), which is exactly what feeding ``h_OH`` to ``build_paper_matrix`` does.

Two equivalent-to-0.02-eV assemblies are provided via ``offdiag``:
  * ``offdiag='hOH'`` (default, PRODUCTION): feed ``h_OH`` globally, drop every v_xc
    patch -- literally ``build_paper_matrix(h_OH, eri, c_H, vxc_diag=None)``.  This is
    the validated ``f0_diagonal_check`` route.  Every EXCHANGE off-diagonal (the
    -c_HF(pr|sq) spin-flip couplings, the W_dd +c_HF K, the V_c 0OS<->4OS exchange)
    is a pure two-electron (nh=2) coupling and is therefore UNCHANGED; only the
    one-electron (nh=1) single-excitation couplings pick up the small delta_h
    off-diagonal (<=0.2 eV).
  * ``offdiag='hact'`` (STRICT off-diagonal reading): keep EVERY off-diagonal at the
    h_act value and shift only the (all-determinant) diagonal by delta_h.  Agrees
    with 'hOH' to ~0.02 eV; provided for the "keep off-diagonals unchanged" reading.

At c_HF = 1 with a pure-HF (ROKS->ROHF, v_xc=0) reference, ``h_OH == h_act`` to
machine precision, so BOTH routes are INERT at the c_HF=1 gate
(``build_paper_matrix_f0 == build_paper_matrix`` to <1e-12).

Engine: reuses ``qmrsf_dk_matrices_ref``, ``qmrsf_dk_pyscf`` (active_integrals/
active_vxc), ``qmrsf_dk_paper`` (build_paper_matrix / _paper_HDK).
"""
from __future__ import annotations

import numpy as np

import qmrsf_dk_matrices_ref as ref
import qmrsf_dk_pyscf as Q
import qmrsf_dk_paper as P

NORB = ref.NORB          # 4 active orbitals
NDET = ref.NDET          # 36 M_s=0 determinants
HARTREE2EV = 27.211386245988


# ---------------------------------------------------------------------------
# F^(0) ground-state KS Fock projected to the 4 SOMOs (reuses f0_diagonal_check)
# ---------------------------------------------------------------------------
def reference_cHF(mf):
    r"""The REFERENCE functional's OWN exact-exchange fraction c_ref (1.0 for HF/ROHF,
    ``dft.libxc.hybrid_coeff(mf.xc)`` for a hybrid; 0.5 for BHHLYP).

    This is the coefficient that dresses F^(0)/h_OH -- i.e. the FIXED relaxed reference
    orbitals and orbital energies eps.  It is DISTINCT from the diagnostic c_H that
    multiplies the ACTIVE-space matrix exchange K.  ***LAW***: the diagnostic c_H must
    NEVER enter F^(0).  Moving c_H into F^(0) would move the core-active exchange and
    hence the relaxed MO energies eps that the CSF forms and the orbital Hessian are
    built on -- and the MOs must not change.  (For a no-core system h_OH is already
    c_H-independent, so the discrepancy is exactly the core-active exchange in F^(0).)
    For a pure-HF reference c_ref == 1, so at the c_H=1 gate c_ref == c_H and nothing
    changes -- the HF-ref FCI/CASCI exactness check is untouched."""
    xc = getattr(mf, "xc", None)
    if xc is None or str(xc).strip().lower() in ("hf", "hartree-fock", ""):
        return 1.0
    from pyscf import dft
    try:
        return float(dft.libxc.hybrid_coeff(str(xc), mf.mol.spin))
    except TypeError:
        return float(dft.libxc.hybrid_coeff(str(xc)))


def ks_fock_active(mf, c_H, somo):
    r"""Ground-state ROKS-BHHLYP KS Fock ``F^(0),sigma = h + J[rho] - c_H K[rho^sigma]
    + v_xc^sigma[rho]`` projected onto the 4 SOMOs.  Returns ``(Fa_act, Fb_act)``.
    Robust to a pure-HF (ROHF) reference (v_xc = 0)."""
    mol = mf.mol
    mo = mf.mo_coeff
    occ = np.asarray(mf.mo_occ, float)
    dma = (mo * (occ > 0)) @ mo.T
    dmb = (mo * (occ == 2)) @ mo.T
    vj, vk = mf.get_jk(mol, np.array([dma, dmb]))
    Jtot = vj[0] + vj[1]
    xc = getattr(mf, "xc", None)
    if xc is None or str(xc).strip().lower() in ("hf", "hartree-fock", ""):
        vxc = (np.zeros_like(Jtot), np.zeros_like(Jtot))
    else:
        if mf.grids.coords is None:
            mf.grids.build()
        _, _, vxc = mf._numint.nr_uks(mol, mf.grids, xc, (dma, dmb))
    hcore = mf.get_hcore()
    Fa = hcore + Jtot - c_H * vk[0] + np.asarray(vxc[0])
    Fb = hcore + Jtot - c_H * vk[1] + np.asarray(vxc[1])
    Cs = mo[:, somo]
    return Cs.T @ Fa @ Cs, Cs.T @ Fb @ Cs


def h_OH_operator(mf, c_H, somo, eri_act):
    r"""F^(0)-anchored spin-averaged effective one-electron operator (4x4).
    ``h_OH -> h_act`` for a pure-HF reference at c_H=1."""
    Fa, Fb = ks_fock_active(mf, c_H, somo)
    J = np.einsum("pqkk->pq", eri_act)          # sum_k (pq|kk)
    Kx = np.einsum("pkkq->pq", eri_act)         # sum_k (pk|kq)
    amf_a = J - c_H * Kx
    amf_b = J
    return 0.5 * ((Fa - amf_a) + (Fb - amf_b))


# ---------------------------------------------------------------------------
# the production matrix builder
# ---------------------------------------------------------------------------
def _vxc_full_f0(mf, somo):
    r"""Full 4x4 spin-resolved semilocal v_xc on the SOMOs (va_mat, vb_mat); zero for HF.
    Used to strip h_OH's spin-AVERAGED v_xc so the paper's spin-RESOLVED DIAGONAL patch
    (eq:vxc1/eq:vxc0os) can be added back -- matches ``build_paper_matrix_f0_direct``."""
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


def build_paper_matrix_f0(mf, nc, somo, c_H,
                          h_act=None, eri_act=None, vxc_diag=None,
                          offdiag="hOH", e_ref=0.0, vxc_singles="resolved"):
    r"""Build the three spin-pure QMRSF-DK blocks with the TRUE F^(0) diagonal.

    Same return shape as ``build_paper_matrix`` (``A_singlet`` 20x20,
    ``A_triplet`` 15x15, ``A_quintet`` 1x1, plus U/Sval/labels/H/H_DK/S2/c_H/osc),
    with extra diagnostics ``h_OH`` and ``delta_h_diag``.

    Parameters
    ----------
    mf   : converged quintet ROKS (or ROHF) mean-field.
    nc   : number of doubly-occupied (core) MOs.
    somo : the 4 SOMO indices.
    c_H  : exact-exchange fraction (0.5 BHHLYP, 1.0 HF gate).
    h_act, eri_act, vxc_diag : optional precomputed active integrals / v_xc
        (recomputed from ``mf`` if omitted; supply them to keep a gauge scan
        apples-to-apples with the value/paper builders).
    offdiag : 'hOH' (default/production, feed h_OH globally) or 'hact' (strict: keep
        every off-diagonal at h_act, shift only the diagonal).
    vxc_singles : 'resolved' (DEFAULT, paper-faithful): the single-flip/0OS v_xc is the
        spin-RESOLVED DIAGONAL patch (eq:vxc1 <a|v^b|a>-<i|v^a|i>, eq:vxc0os
        sum_occ v^b - sum_empty v^a); h_OH's spin-averaged v_xc is stripped first.
        'averaged': legacy spin-AVERAGED 1/2(v^a+v^b) inside h_OH, no patch (the
        pre-2026-07 published build).  Both are inert at the c_H=1 gate (v_xc=0).
        Only offdiag='hOH' honours vxc_singles.
    """
    if h_act is None or eri_act is None:
        h_act, eri_act, _ec, _ci, _mc = Q.active_integrals(mf, nc, somo)
    h_act = np.asarray(h_act, float)
    eri_act = np.asarray(eri_act, float)
    # F^(0)/h_OH is built at the REFERENCE functional's exchange c_ref (fixed MOs), NOT
    # the diagnostic c_H -- see reference_cHF().  The diagnostic c_H below enters ONLY the
    # active-space matrix exchange (P.build_paper_matrix / _paper_HDK).
    h_OH = h_OH_operator(mf, reference_cHF(mf), somo, eri_act)
    ddiag = np.diag(np.asarray(h_OH, float) - h_act)     # length-4

    if offdiag == "hOH":
        if vxc_singles == "resolved":
            # PAPER-FAITHFUL (eq:vxc1 singles / eq:vxc0os doubles): strip h_OH's
            # spin-AVERAGED v_xc and add the spin-RESOLVED DIAGONAL patch (na-1)v^a+nb v^b.
            # Bit-for-bit identical to build_paper_matrix_f0_direct(vxc_singles='resolved').
            va_d, vb_d = Q.active_vxc(mf, mf.mo_coeff[:, somo])
            va_m, vb_m = _vxc_full_f0(mf, somo)
            dk = P.build_paper_matrix(h_OH - 0.5 * (va_m + vb_m), eri_act, c_H,
                                      vxc_diag=(va_d, vb_d), e_ref=e_ref)
        elif vxc_singles == "averaged":
            # legacy: feed F^(0) operator globally, v_xc spin-AVERAGED inside eps, no patch.
            dk = P.build_paper_matrix(h_OH, eri_act, c_H, vxc_diag=None, e_ref=e_ref)
        else:
            raise ValueError(f"vxc_singles must be 'resolved' or 'averaged', got {vxc_singles!r}")
    elif offdiag == "hact":
        # STRICT off-diagonal reading: off-diagonals stay at h_act, only the
        # (all-determinant) diagonal is shifted by delta_h; no v_xc patch.
        dets = ref.gen_dets()
        H, S2 = ref.build_H_S2(h_act, eri_act, dets)
        U, Sval, labels = ref.build_U(dets)
        H_DK, osc = P._paper_HDK(h_act, eri_act, c_H, dets, H,
                                 vxc_diag=None, e_ref=e_ref)
        for d, det in enumerate(dets):
            ntot = np.array([(1.0 if p in det else 0.0)
                             + (1.0 if (p + NORB) in det else 0.0)
                             for p in range(NORB)])
            H_DK[d, d] += float(np.dot(ntot, ddiag))
        Hcsf = U @ H_DK @ U.T
        sel = {S: np.where(Sval == S)[0] for S in (0, 1, 2)}
        A = {S: Hcsf[np.ix_(sel[S], sel[S])] for S in (0, 1, 2)}
        dk = {"A_singlet": A[0], "A_triplet": A[1], "A_quintet": A[2],
              "U": U, "Sval": Sval, "labels": labels,
              "H": H, "H_DK": H_DK, "S2": S2, "c_H": c_H, "osc": osc}
    else:
        raise ValueError(f"offdiag must be 'hOH' or 'hact', got {offdiag!r}")

    dk["h_OH"] = h_OH
    dk["delta_h_diag"] = ddiag
    return dk


# ---------------------------------------------------------------------------
# c_H = 1 inert gate (must pass; if it breaks the F^(0) construction is wrong)
# ---------------------------------------------------------------------------
def verify_cH1_gate(mf, nc, somo, tol=1e-12, offdiag="hOH"):
    """(a) With a pure-HF (ROHF) reference at c_H=1, ``build_paper_matrix_f0`` must
    equal ``build_paper_matrix`` (and the active-space CASCI spectrum) to < tol."""
    h_act, eri_act, ecore, casci, _mc = Q.active_integrals(mf, nc, somo)
    vxc = Q.active_vxc(mf, mf.mo_coeff[:, somo])          # zeros for HF
    dk_p = P.build_paper_matrix(h_act, eri_act, 1.0, vxc_diag=vxc)
    dk_f = build_paper_matrix_f0(mf, nc, somo, 1.0, h_act=h_act, eri_act=eri_act,
                                 offdiag=offdiag)
    blk = max(float(np.abs(dk_f[k] - dk_p[k]).max())
              for k in ("A_singlet", "A_triplet", "A_quintet"))
    hoh = float(np.abs(dk_f["h_OH"] - h_act).max())
    ci_err = np.nan
    if casci is not None:
        spec = np.sort(np.concatenate([
            np.linalg.eigvalsh(dk_f["A_singlet"]),
            np.linalg.eigvalsh(dk_f["A_triplet"]),
            np.atleast_1d(np.linalg.eigvalsh(np.atleast_2d(dk_f["A_quintet"]))),
        ]) + ecore)
        ci = np.sort(np.asarray(casci))
        if spec.shape == ci.shape:
            ci_err = float(np.abs(spec - ci).max())
    return {"block_f0_vs_paper": blk, "hOH_vs_hact": hoh,
            "spectrum_vs_CASCI": ci_err,
            "pass": bool(blk < tol and hoh < tol)}


if __name__ == "__main__":
    atom = ("H 0.7 0.7 0\nH -0.7 0.7 0\nH -0.7 -0.7 0\nH 0.7 -0.7 0")
    mf, cH = Q.build_reference(atom, "cc-pvdz", "hf")
    somo, nc = Q.active_space(mf)
    print("c_H=1 INERT GATE (ROHF H4):")
    for od in ("hOH", "hact"):
        g = verify_cH1_gate(mf, nc, somo, offdiag=od)
        print(f"  [offdiag={od:4s}] |A_f0-A_paper|={g['block_f0_vs_paper']:.2e} "
              f"|h_OH-h_act|={g['hOH_vs_hact']:.2e} "
              f"vsCASCI={g['spectrum_vs_CASCI']:.2e} "
              f"-> {'PASS' if g['pass'] else 'FAIL'}")
        assert g["pass"], f"c_H=1 gate FAILED for offdiag={od}"
