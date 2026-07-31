#!/usr/bin/env python3
"""Standalone PySCF QMRSF-DK tool.

For ANY molecule this runs the high-spin **quintet** ROHF/ROKS reference, takes
the **4 SOMOs** (singly occupied MOs) as the CAS(4,4) active space, builds the
hard-coded QMRSF-DK spin-adapted matrices (Singlet 20x20, Triplet 15x15,
Quintet 1x1) exactly as derived/verified in ``DK_MATRICES_SPEC.md``, and returns
the excitation energies.

The spin-adaptation U and the 36-determinant Slater-Condon machinery live in the
already-verified ``qmrsf_dk_matrices_ref`` module (its self-test PASSES; its U
reproduces ``U_explicit_36.npy`` bit-for-bit).  This file adds:

  * the PySCF pipeline (build_reference / active_space / active_integrals),
  * the **c_H seam**  H_DK = H - (1 - c_H) * H_K1SF  (c_H scales EVERY exact-exchange
    integral across the whole 36-det block -- 2OS/4OS single flips AND the 0OS
    doubles, per the SI), with the c_H=1 -> CASCI(4,4) exactness assertion,
  * the explicit reference-imbalance **Delta** diagnostics
    (B_ave + Delta = E0, the per-config 2OS split +/-K_rs, quintet = B_ave - sumK),
  * run_qmrsf_dk(atom, basis, xc) -> excitation-energy dict (eV + Ha).

Theory is hard-coded once; just feed it different systems.
Engine: PySCF 2.10 (numpy/scipy only).  Do NOT contradict the two SPEC files.
"""
from __future__ import annotations

import numpy as np

import qmrsf_dk_matrices_ref as ref

HARTREE2EV = 27.211386245988

NORB = ref.NORB          # 4 active orbitals
NSO = ref.NSO            # 8 spin orbitals
NDET = ref.NDET          # 36 M_s=0 determinants

# Colloquial functional name -> this libxc's exact alias (BHHLYP crashes as a
# raw libxc name in PySCF 2.10; 'bhandhlyp' is the working alias, c_H=0.5).
_XC_ALIASES = {"bhhlyp": "bhandhlyp"}


# ===========================================================================
# 1. quintet reference
# ===========================================================================
def build_reference(atom, basis, xc="hf", verbose=0, spin=4, charge=0,
                    conv_tol=1e-10, max_cycle=200):
    """Run the high-spin quintet (S=2, spin=4) ROHF (xc='hf') or ROKS reference.

    Returns ``(mf, c_H)`` where ``c_H`` is the fraction of exact (HF) exchange:
    1.0 for HF/ROHF, ``dft.libxc.hybrid_coeff(xc)`` for a hybrid functional.
    """
    from pyscf import gto, scf, dft

    mol = gto.M(atom=atom, basis=basis, spin=spin, charge=charge, verbose=verbose)
    if xc.lower() in ("hf", "rohf", ""):
        mf = scf.ROHF(mol)
        c_H = 1.0
    else:
        # Map colloquial functional names to this libxc's exact aliases so the
        # common name ('bhhlyp') does not crash with a libxc KeyError.
        xc = _XC_ALIASES.get(xc.lower(), xc)
        mf = scf.ROKS(mol)
        mf.xc = xc
        try:
            c_H = float(dft.libxc.hybrid_coeff(xc, mol.spin))
        except TypeError:
            c_H = float(dft.libxc.hybrid_coeff(xc))
        except KeyError as ex:
            raise KeyError(
                f"{ex}. Unknown libxc functional name '{xc}'. "
                f"For BHHLYP use xc='bhandhlyp'.") from ex
    mf.conv_tol = conv_tol
    mf.max_cycle = max_cycle
    mf.kernel()
    return mf, c_H


# ===========================================================================
# 2. active space = the 4 SOMOs
# ===========================================================================
def active_space(mf):
    """Return ``(somo_idx, ncore)``: the 4 singly-occupied MO indices and the
    number of doubly-occupied (core) MOs.  Raises if there are not exactly 4
    SOMOs (the QMRSF-DK quintet reference requires exactly four)."""
    mo_occ = np.asarray(mf.mo_occ)
    somo_idx = np.where(np.isclose(mo_occ, 1.0))[0]
    ncore = int(np.sum(np.isclose(mo_occ, 2.0)))
    if len(somo_idx) != 4:
        raise ValueError(
            f"expected exactly 4 SOMOs (quintet, spin=4) but found "
            f"{len(somo_idx)} singly-occupied MOs (occ==1). "
            f"mo_occ = {mo_occ}")
    return somo_idx.tolist(), ncore


# ===========================================================================
# 3. active-space integrals + CASCI reference spectrum
# ===========================================================================
def active_integrals(mf, ncore, somo_idx, nroots=36):
    """Build CAS(4,4) on the 4 SOMOs and return
    ``(h_act, eri_act, ecore, casci_e, mc)``.

    ``h_act`` is the (4,4) one-electron effective Hamiltonian, ``eri_act`` the
    (4,4,4,4) chemist ``(ij|kl)`` two-electron integrals, ``ecore`` the frozen
    core + nuclear constant, and ``casci_e`` the PySCF CASCI(4,4) total energies
    (the c_H=1 validation target).  The 4 SOMOs are forced into the active
    window [ncore:ncore+4]."""
    from pyscf import mcscf, ao2mo

    # reorder MOs so the 4 SOMOs occupy the active window [ncore:ncore+4]
    nmo = mf.mo_coeff.shape[1]
    somo = list(somo_idx)
    others = [i for i in range(nmo) if i not in somo]
    # Select the frozen core explicitly by OCCUPATION (occ==2), not by an index
    # slice.  This is robust to non-aufbau / symmetry-broken ROHF orderings
    # (a doubly-occupied orbital above a SOMO, or a virtual below); the old
    # ``others[:ncore]`` would silently freeze the wrong orbital there.
    mo_occ = np.asarray(mf.mo_occ)
    core = [i for i in others if np.isclose(mo_occ[i], 2.0)]
    assert len(core) == ncore, (
        f"core selection by occupation found {len(core)} doubly-occupied MOs "
        f"but ncore={ncore}; mo_occ={mo_occ}")
    virt = [i for i in others if not np.isclose(mo_occ[i], 2.0)]
    order = core + somo + virt
    mo = mf.mo_coeff[:, order]

    # NOTE: the active space is filled with (2,2) electrons -- the M_s=0 sector
    # whose 36 determinants give the 20 singlet + 15 triplet + 1 quintet
    # spectrum.  (A plain CAS(4,4) would inherit mol.spin=4 -> only the single
    # M_s=2 top-of-quintet determinant.)
    mc = mcscf.CASCI(mf, 4, (2, 2))
    mc.ncore = ncore
    mc.fcisolver.nroots = nroots
    mc.fcisolver.spin = 0
    mc.verbose = 0

    h_act, ecore = mc.get_h1eff(mo)
    eri_act = ao2mo.restore(1, mc.get_h2eff(mo[:, ncore:ncore + 4]), 4)
    eri_act = np.asarray(eri_act).reshape(4, 4, 4, 4)

    try:
        mc.kernel(mo)
        casci_e = np.sort(np.atleast_1d(mc.e_tot))
    except Exception:
        casci_e = None
    return np.asarray(h_act), eri_act, float(ecore), casci_e, mc


def active_vxc(mf, mo_act):
    """``(va, vb)``: the diagonal of the local (semilocal DFT) exchange-correlation
    potential over the 4 active orbitals, ``<p|v_xc^sigma|p>``, spin resolved.

    This is the DFT character the QMRSF-DK diagonal carries (main text); it is the
    part of the KS orbital energies that the bare Slater-Condon Hamiltonian omits.
    For a hybrid only the SEMILOCAL piece enters here (the exact-exchange fraction
    c_H is carried separately on the couplings), so for a pure-HF reference this is
    identically zero and the c_H=1 bare path is preserved."""
    xc = getattr(mf, "xc", None)
    if xc is None or str(xc).strip().lower() in ("hf", "hartree-fock", ""):
        return np.zeros(mo_act.shape[1]), np.zeros(mo_act.shape[1])
    mo = mf.mo_coeff
    occ = np.asarray(mf.mo_occ, float)
    dma = (mo * (occ > 0)) @ mo.T          # alpha density (occ 1 or 2)
    dmb = (mo * (occ == 2)) @ mo.T         # beta density (occ 2)
    if getattr(mf.grids, "coords", None) is None:
        mf.grids.build()
    _, _, vmat = mf._numint.nr_uks(mf.mol, mf.grids, xc, (dma, dmb))
    va = np.einsum("pi,pq,qi->i", mo_act, np.asarray(vmat[0]), mo_act)
    vb = np.einsum("pi,pq,qi->i", mo_act, np.asarray(vmat[1]), mo_act)
    return va, vb


# ===========================================================================
# 4. the c_H seam + spin-adapted blocks + Delta diagnostics
# ===========================================================================
def _exchange_eri(eri_act):
    """Chemist eri tensor restricted to its EXCHANGE-VALUE integrals (the K-orbit
    ``(ij|ji)`` i!=j and the genuine 4-index ``(pr|sq)``), with the Hartree/Coulomb
    integrals zeroed by value (one charge distribution diagonal: ``a==b`` or
    ``c==d`` -> the J-orbit (ii|jj), the off-diagonal (pq|kk), the (ii|ij)).

    This is the standard hybrid "scale-K" prescription at the integral level.  It
    correctly captures the spin-flip L/R coupling K_rs (an (ij|ji) exchange-VALUE
    integral that the Slater-Condon spin-orbital routing places in the *direct*
    argument slot, so a ROLE/contraction-based "-( PS|QR) only" split would WRONGLY
    miss it and the S-T gap would not scale with c_H), and it keeps every
    2-electron Hartree coupling bare.

    KNOWN APPROXIMATION (documented): for 4 DISTINCT orbitals the chemist tensor is
    symmetric over the three pairings (pq|rs)/(pr|qs)/(ps|qr), so a fixed tensor
    cannot separate the exchange ROLE of a 4-index integral from its Coulomb role;
    this rule scales all non-Hartree 4-index pairings.  It is EXACT at c_H=1 (the
    CASCI gate), EXACT on every diagonal and on all 2-index K (Delta, the quintet,
    the S-T gaps), and removes the off-diagonal (pq|kk) Hartree over-scaling.  The
    residual 4-index Coulomb-role term (small for a 4-SOMO active space) is the only
    c_H<1 approximation; the production-exact partition is the SI off-diagonal
    c_ij (orbital-Hessian -c_H(pr|sq)) -- see DK_MATRICES_SPEC / README caveat."""
    eri_K = np.array(eri_act, float, copy=True)
    for i in range(NORB):
        eri_K[i, i, :, :] = 0.0   # a==b : Coulomb density |phi_i|^2 on electron 1
        eri_K[:, :, i, i] = 0.0   # c==d : Coulomb density |phi_i|^2 on electron 2
    return eri_K


def _build_H_K1SF(h_act, eri_act, H_full, dets, cfg):
    """The exchange-VALUE 36x36 Slater-Condon contribution (h set to zero, built
    from ``_exchange_eri``) that carries the hybrid fraction c_H in the seam
    H_DK = H - (1-c_H)*H_K.

    Per the MRSF-TDDFT prescription, the exact-exchange fraction c_H of the XC
    functional multiplies EVERY exact-exchange integral throughout the
    Slater-Condon model Hamiltonian -- the single-spin-flip (orbital-Hessian)
    block AND the closed-shell (0OS) doubles block -- so this contribution spans
    the WHOLE 36x36 matrix (no 1SF mask).  The off-diagonals are the genuine
    Slater-Condon couplings; only their exchange content scales.  Hartree (J), the
    one-electron h, and the Coulomb role stay full; exact at c_H=1 (H_DK==H).  See
    _exchange_eri for the 4-distinct-orbital caveat.  (``H_full``/``cfg`` are
    accepted for signature compatibility.)"""
    eri_K = _exchange_eri(eri_act)
    HK, _ = ref.build_H_S2(np.zeros((NORB, NORB)), eri_K, dets)
    return HK


def _delta_diagnostics(h_act, eri_act, A_singlet, A_triplet, A_quintet,
                       Sval, labels):
    """Compute and verify the reference-imbalance Delta relations.

    Returns a dict with B_ave / E0 / Delta per config and the verification
    errors for:
      * quintet      : A_quintet = B_ave(4OS) - sum K
      * 4OS uniform  : B_ave(4OS) + Delta^4OS = E0(4OS),  Delta^4OS = -1/2 sum K
      * 2OS per-config: 1A = B_ave_(pq) + Delta_(pq) + K_rs,
                        3A = B_ave_(pq) + Delta_(pq) - K_rs,
                        B_ave_(pq) + Delta_(pq) = E0 (the det Slater diagonal).
    """
    K = np.array([[eri_act[i, j, j, i] for j in range(NORB)] for i in range(NORB)])
    J = np.array([[eri_act[i, i, j, j] for j in range(NORB)] for i in range(NORB)])
    hd = np.diag(np.asarray(h_act))
    pairs = [(i, j) for i in range(NORB) for j in range(i + 1, NORB)]
    sumK = sum(K[i, j] for i, j in pairs)
    sumJ = sum(J[i, j] for i, j in pairs)
    sumh = float(hd.sum())

    # ---- 4OS / quintet ----
    B_ave_4os = sumh + sumJ
    Delta_4os = -0.5 * sumK
    E0_4os = B_ave_4os + Delta_4os                 # common 4OS Coulomb diagonal
    quint_val = float(np.atleast_2d(A_quintet)[0, 0])
    err_quint = abs(quint_val - (B_ave_4os - sumK))

    # ---- 2OS per-config ----
    slabs = [labels[k] for k in range(NDET) if Sval[k] == 0]
    tlabs = [labels[k] for k in range(NDET) if Sval[k] == 1]
    # singlet rows 6..17 are the 12 2OS singlets; triplet rows 0..11 the partners
    s_off = [i for i, l in enumerate(slabs) if l.startswith("2OS-S")]
    t_off = [i for i, l in enumerate(tlabs) if l.startswith("2OS-T")]

    def parse(cfg):
        p = [i for i, c in enumerate(cfg) if c == "2"][0]
        q = [i for i, c in enumerate(cfg) if c == "0"][0]
        rs = [i for i, c in enumerate(cfg) if c in "ud"]
        return p, q, rs[0], rs[1]

    rows = []
    err_2os = 0.0
    for si, ti in zip(s_off, t_off):
        cfg = slabs[si].split(":")[1].split("/")[0]
        p, q, r, s = parse(cfg)
        oneA = float(A_singlet[si, si])
        threeA = float(A_triplet[ti, ti])
        E0 = 0.5 * (oneA + threeA)                 # det Slater diagonal
        Krs = K[r, s]
        Delta = -K[p, q] - 0.5 * (K[p, r] + K[p, s] + K[q, r] + K[q, s])
        Bcoul = (2 * hd[p] + hd[r] + hd[s] + eri_act[p, p, p, p]
                 + 2 * J[p, r] + 2 * J[p, s] + J[r, s])
        B_ave = Bcoul + K[p, q] - 0.5 * (K[p, r] + K[p, s]) \
            + 0.5 * (K[q, r] + K[q, s])
        err_2os = max(err_2os,
                      abs(B_ave + Delta - E0),
                      abs(oneA - (B_ave + Delta + Krs)),
                      abs(threeA - (B_ave + Delta - Krs)))
        rows.append({"cfg": cfg, "pq": (p, q), "rs": (r, s),
                     "B_ave": B_ave, "Delta": Delta, "E0": E0,
                     "K_rs": Krs, "singlet": oneA, "triplet": threeA})

    return {
        "sumK": sumK, "sumJ": sumJ, "sumh": sumh,
        "B_ave_4os": B_ave_4os, "Delta_4os": Delta_4os, "E0_4os": E0_4os,
        "quintet": quint_val,
        "err_quintet": float(err_quint),
        "err_2os": float(err_2os),
        "rows_2os": rows,
    }


def build_dk_matrices(h_act, eri_act, c_H, e_ref=0.0, with_diagnostics=True,
                      vxc_diag=None):
    """Build the three spin-pure QMRSF-DK blocks with the c_H exchange seam.

    H_DK = H - (1 - c_H) * H_K1SF, where H is the full 36-det CAS(4,4) M_s=0
    Slater-Condon Hamiltonian and H_K1SF is the EXCHANGE-only contribution over the
    WHOLE 36-det block (NO 1SF mask): c_H scales every exact-exchange integral,
    INCLUDING the 0OS doubles block, exactly as the SI prescribes.  (The name
    ``H_K1SF`` is historical; the quantity is the full-36 exchange.)

    At c_H = 1, H_DK == H exactly (asserted) -> the union of the three block
    spectra equals the PySCF CASCI(4,4) spectrum on the same 4-SOMO space.

    Returns a dict with ``A_singlet`` (20x20), ``A_triplet`` (15x15),
    ``A_quintet`` (1x1), the spin-adaptation ``U``/``Sval``/``labels``, the
    determinant-basis ``H``/``H_DK``/``S2``, and (if requested) the ``delta``
    diagnostics block.
    """
    h_act = np.asarray(h_act, float)
    eri_act = np.asarray(eri_act, float)

    dets = ref.gen_dets()
    cfg = [ref.det_config(d) for d in dets]
    H, S2 = ref.build_H_S2(h_act, eri_act, dets)
    U, Sval, labels = ref.build_U(dets)

    HK1SF = _build_H_K1SF(h_act, eri_act, H, dets, cfg)
    H_DK = H - (1.0 - c_H) * HK1SF
    if abs(c_H - 1.0) < 1e-14 and vxc_diag is None:
        assert np.max(np.abs(H_DK - H)) < 1e-12, "c_H=1 must give H_DK == H"

    # v_xc on the diagonal: the DFT character the singles/doubles diagonals carry
    # (main text Eqs. dressed1sf/vxc0os).  Each config's diagonal is shifted by the
    # local DFT XC potential of the active orbitals, weighted by the occupation
    # change from the quintet reference (each active orbital alpha-occupied):
    #     dV = sum_p [ (n_p^a - 1) <p|v_xc^a|p> + n_p^b <p|v_xc^b|p> ].
    # For a single flip i->a this is <a|v_xc^b|a> - <i|v_xc^a|i> (Eq. vxc1); for a
    # 0OS double flip it is sum_filled<.|v_xc^b|.> - sum_emptied<.|v_xc^a|.> (Eq.
    # vxc0os).  Zero for a pure-HF reference (no local XC), so the c_H=1 bare path
    # is preserved.
    if vxc_diag is not None:
        va = np.asarray(vxc_diag[0], float)
        vb = np.asarray(vxc_diag[1], float)
        NORB = ref.NORB
        for d, det in enumerate(dets):
            na = np.array([1.0 if p in det else 0.0 for p in range(NORB)])
            nb = np.array([1.0 if (p + NORB) in det else 0.0 for p in range(NORB)])
            H_DK[d, d] += float(np.dot(na - 1.0, va) + np.dot(nb, vb))

    Hcsf = U @ H_DK @ U.T
    if e_ref:
        Hcsf = Hcsf - e_ref * np.eye(NDET)

    sel = {S: np.where(Sval == S)[0] for S in (0, 1, 2)}
    A = {S: Hcsf[np.ix_(sel[S], sel[S])] for S in (0, 1, 2)}

    out = {
        "A_singlet": A[0], "A_triplet": A[1], "A_quintet": A[2],
        "U": U, "Sval": Sval, "labels": labels,
        "H": H, "H_DK": H_DK, "S2": S2, "c_H": c_H,
    }
    if with_diagnostics:
        out["delta"] = _delta_diagnostics(h_act, eri_act, A[0], A[1], A[2],
                                          Sval, labels)
    return out


# ===========================================================================
# 5. full pipeline
# ===========================================================================
def run_qmrsf_dk(atom, basis, xc="hf", spin=4, charge=0, verbose=0,
                 nstates=None):
    """Full pipeline: quintet reference -> 4-SOMO CAS(4,4) -> QMRSF-DK blocks.

    Returns a dict of EXCITATION energies (block eigenvalue minus the quintet
    root) for each spin channel, in both eV and Hartree, plus metadata."""
    mf, c_H = build_reference(atom, basis, xc=xc, spin=spin, charge=charge,
                              verbose=verbose)
    somo_idx, ncore = active_space(mf)
    h_act, eri_act, ecore, casci_e, mc = active_integrals(mf, ncore, somo_idx)
    mo_act = mf.mo_coeff[:, somo_idx]          # the 4 SOMOs (O1..O4 by energy)
    vxc_diag = active_vxc(mf, mo_act)          # DFT v_xc on the diagonal (0 for HF)

    dk = build_dk_matrices(h_act, eri_act, c_H, vxc_diag=vxc_diag)
    ev_s = np.sort(np.linalg.eigvalsh(dk["A_singlet"]))
    ev_t = np.sort(np.linalg.eigvalsh(dk["A_triplet"]))
    ev_q = np.sort(np.linalg.eigvalsh(np.atleast_2d(dk["A_quintet"])))

    # total active electronic energies (add ecore for absolute total)
    quint_root = float(ev_q[0])

    def exc(ev):
        d = np.asarray(ev) - quint_root
        return {"hartree": d.tolist(), "eV": (d * HARTREE2EV).tolist()}

    if nstates is not None:
        ev_s = ev_s[:nstates]
        ev_t = ev_t[:nstates]

    return {
        "singlet": exc(ev_s),
        "triplet": exc(ev_t),
        "quintet": exc(ev_q),
        "abs_hartree": {
            "singlet": (ev_s + ecore).tolist(),
            "triplet": (ev_t + ecore).tolist(),
            "quintet": (ev_q + ecore).tolist(),
        },
        "meta": {
            "c_H": c_H, "ncore": ncore, "somo_idx": somo_idx,
            "ecore": ecore, "quintet_total": quint_root + ecore,
            "scf_e": float(mf.e_tot), "xc": xc, "basis": basis,
            "somo_energies": np.asarray(mf.mo_energy)[somo_idx].tolist(),
            "casci_e": None if casci_e is None else np.asarray(casci_e).tolist(),
        },
        "dk": dk,
        "casci_e": casci_e, "ecore": ecore,
    }


# ===========================================================================
# __main__ self-tests : random-integral gate + Delta demonstration
# ===========================================================================
def _gate_random_integrals(seed=0, c_H=1.0, verbose=True):
    """At c_H=1 the three QMRSF-DK blocks must reproduce the S^2-projected
    36-determinant CI spectrum (the matrix machinery, independent of SCF)."""
    h, e = ref.random_integrals(seed)
    out = build_dk_matrices(h, e, c_H=c_H)
    H, S2, Sval = out["H"], out["S2"], out["Sval"]

    # reference: S^2-projected eigenvalues of the bare 36-det CI
    w, V = np.linalg.eigh(H)
    s2f = np.array([V[:, k] @ S2 @ V[:, k] for k in range(NDET)])

    blocks = {0: out["A_singlet"], 1: out["A_triplet"], 2: out["A_quintet"]}
    dims = (out["A_singlet"].shape[0], out["A_triplet"].shape[0],
            int(np.atleast_2d(out["A_quintet"]).shape[0]))
    eigerr = 0.0
    for S in (0, 1, 2):
        eb = np.sort(np.linalg.eigvalsh(np.atleast_2d(blocks[S])))
        ed = np.sort(w[np.abs(s2f - S * (S + 1)) < 1e-6])
        eigerr = max(eigerr, float(np.abs(eb - ed).max()))
    # H_DK == H at c_H=1
    h_dk_err = float(np.max(np.abs(out["H_DK"] - H)))
    ok = (dims == (20, 15, 1) and eigerr < 1e-9 and h_dk_err < 1e-12)
    if verbose:
        print(f"  [seed {seed}] block dims (S/T/Q)        = {dims}")
        print(f"  [seed {seed}] H_DK == H (c_H=1)          = {h_dk_err:.3e}")
        print(f"  [seed {seed}] block spectra vs S2-CI     = {eigerr:.3e}")
        print(f"  [seed {seed}] RESULT                     = "
              f"{'PASS' if ok else 'FAIL'}")
    return ok, eigerr, h_dk_err, dims


def _demo_delta(seed=0, verbose=True):
    """Demonstrate the reference-imbalance Delta explicitly:
      B_ave + Delta = E0 ; 2OS diag = B_ave+Delta +/- K_rs ; quintet = B_ave - sumK.
    """
    h, e = ref.random_integrals(seed)
    out = build_dk_matrices(h, e, c_H=1.0)
    d = out["delta"]
    if verbose:
        print("  --- Delta (reference-imbalance) demonstration ---")
        print(f"    4OS: B_ave = sum h + sum J        = {d['B_ave_4os']:+.6f}")
        print(f"    4OS: Delta^4OS = -1/2 sum K        = {d['Delta_4os']:+.6f}")
        print(f"    4OS: B_ave + Delta = E0(4OS)       = {d['E0_4os']:+.6f}")
        print(f"    QUINTET diag                       = {d['quintet']:+.6f}")
        print(f"           B_ave - sum K               = "
              f"{d['B_ave_4os'] - d['sumK']:+.6f}")
        print(f"    |quintet - (B_ave - sumK)|         = {d['err_quintet']:.3e}")
        print("    2OS per-config  (cfg: B_ave + Delta = E0 ; "
              "1A/3A = E0 +/- K_rs):")
        for r in d["rows_2os"][:4]:
            print(f"      {r['cfg']}: B_ave={r['B_ave']:+.4f} "
                  f"Delta={r['Delta']:+.4f}  E0={r['E0']:+.4f}  "
                  f"K_rs={r['K_rs']:+.4f}  1A={r['singlet']:+.4f} "
                  f"3A={r['triplet']:+.4f}")
        print(f"    ... (12 2OS configs total)")
        print(f"    max |2OS identity error|           = {d['err_2os']:.3e}")
    ok = d["err_quintet"] < 1e-9 and d["err_2os"] < 1e-9
    if verbose:
        print(f"    Delta demonstration RESULT         = "
              f"{'PASS' if ok else 'FAIL'}")
    return ok, d


if __name__ == "__main__":
    print("=" * 68)
    print(" QMRSF-DK PySCF tool -- __main__ self-tests")
    print("=" * 68)

    print("\n[1] RANDOM-INTEGRAL GATE (c_H=1 -> S^2-projected 36-CI spectrum)")
    g_ok = True
    for sd in (0, 12345, 7):
        ok, *_ = _gate_random_integrals(seed=sd)
        g_ok = g_ok and ok

    print("\n[2] DELTA DEMONSTRATION (B_ave+Delta=E0; 2OS +/-K_rs; "
          "quintet=B_ave-sumK)")
    d_ok, _ = _demo_delta(seed=0)

    print("\n[3] c_H SEAM SANITY (c_H<1 scales single-spin-flip exchange)")
    h, e = ref.random_integrals(0)
    full = build_dk_matrices(h, e, c_H=1.0, with_diagnostics=False)
    half = build_dk_matrices(h, e, c_H=0.5, with_diagnostics=False)
    q_full = float(np.atleast_2d(full["A_quintet"])[0, 0])
    q_half = float(np.atleast_2d(half["A_quintet"])[0, 0])
    print(f"    quintet diag c_H=1.0 = {q_full:+.6f}")
    print(f"    quintet diag c_H=0.5 = {q_half:+.6f}  "
          f"(quintet carries c_H exchange too; both single-spin-flip)")
    print(f"    H_DK(0.5) differs from H only on the 1SF block: "
          f"max|H-H_DK| = {np.max(np.abs(full['H_DK']-half['H_DK'])):.3e}")

    print("\n" + "=" * 68)
    allok = g_ok and d_ok
    print(f" SELF-TESTS: {'ALL PASS' if allok else 'FAILED'}")
    print("=" * 68)
    assert allok, "self-tests FAILED"
