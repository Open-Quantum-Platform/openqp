#!/usr/bin/env python3
"""Hard-coded spin-adapted QMRSF-DK reference matrices (NumPy).

Given the CAS(4,4) active-space integrals (``h_act`` 4x4 and ``eri_act`` 4x4x4x4,
chemist ``(ij|kl)``), this module builds the M_s=0 36-determinant Slater-Condon
Hamiltonian and applies the HARD-CODED analytic spin-adaptation U to return the
three spin-pure QMRSF-DK blocks::

    A_singlet  (20x20)
    A_triplet  (15x15)
    A_quintet  ( 1x1 )      = (U^T H U) partitioned by total spin S.

The spin-adaptation U is hard-coded analytically -- NO numerical S^2
diagonalization is used (cf. ``engine/pspace_spinadapt.py``, the icPT2 sibling):

  * 0OS  (closed shell)   : identity -- the six closed-shell determinants are the
                            six singlets.
  * 2OS  (two singly-occ) : the 12 left/right pairs (partner = u<->d swap) are
                            Hadamard-coupled -- singlet (phi_L - phi_R)/sqrt2,
                            triplet (phi_L + phi_R)/sqrt2 (orbital-ordered).
  * 4OS  (four singly-occ): the six genealogical (Yamanouchi-Kotani) CSFs Q, S1,
                            S2, T1, T2, T3 over the ordered determinant basis
                            [uudd, udud, uddu, duud, dudu, dduu].

Every orbital-ordered CSF is then dressed by the diagonal Fermionic-ordering
phase D = diag(+1,-1,+1,+1,-1,+1) on the 4OS block (more generally the per-
determinant ``(-1)^nu`` reordering sign), mapping the genealogical basis to the
engine's sorted-spin-orbital basis so that U^T H U block-diagonalizes the engine
Hamiltonian directly.

Self-test (run as ``__main__``): on RANDOM active-space integrals it asserts
that U is orthonormal, that U^T H U is block diagonal in spin (cross-blocks
< 1e-10), and that each spin block's spectrum equals the S^2-projected H
spectrum.  All matrices are verified against the from-scratch CAS(4,4) CI in
the spec (scratchpad/DK_MATRICES_SPEC.md); the hard-coded U reproduces the
saved reference U_explicit_36.npy bit-for-bit per spin block.
"""
from __future__ import annotations

import numpy as np

NORB = 4          # active orbitals
NELA = NELB = 2   # alpha / beta electrons (M_s = 0)
NSO = 2 * NORB    # spin orbitals (alpha 0..3, beta 4..7)
NDET = 36         # M_s=0 determinants of CAS(4,4)

SQ2, SQ3, SQ6 = np.sqrt(2.0), np.sqrt(3.0), np.sqrt(6.0)

# 4OS genealogical determinant order (orbital-ordered basis).
FOUR_OS_ORDER = ["uudd", "udud", "uddu", "duud", "dudu", "dduu"]

# Diagonal Fermi phase on the 4OS block, genealogical -> engine(sorted) basis.
# Equals (-1)^nu for each of FOUR_OS_ORDER; quoted explicitly in the spec.
FERMI_D_4OS = np.array([+1.0, -1.0, +1.0, +1.0, -1.0, +1.0])

# Hard-coded genealogical 4OS CSF coefficients (q=1/2, a3=1/sqrt3,
# b3=1/(2 sqrt3), a6=1/sqrt6) over FOUR_OS_ORDER.  S=0 singlets, S=1 triplets,
# S=2 quintet.  These are exactly the Yamanouchi-Kotani CSFs of the spec.
_Q, _A3, _B3, _A6 = 0.5, 1.0 / SQ3, 1.0 / (2.0 * SQ3), 1.0 / SQ6
FOUR_OS_CSF = [
    ("4OS-S1", 0, np.array([0.0,  _Q,  -_Q,  -_Q,   _Q,  0.0])),
    ("4OS-S2", 0, np.array([_A3, -_B3, -_B3, -_B3, -_B3,  _A3])),
    ("4OS-T1", 1, np.array([0.0,  _Q,   _Q,  -_Q,  -_Q,  0.0])),
    ("4OS-T2", 1, np.array([_A3, -_B3,  _B3, -_B3,  _B3, -_A3])),
    ("4OS-T3", 1, np.array([_A6,  _A6, -_A6,  _A6, -_A6, -_A6])),
    ("4OS-Q",  2, np.array([_A6,  _A6,  _A6,  _A6,  _A6,  _A6])),
]


# ---------------------------------------------------------------------------
# determinant bookkeeping
# ---------------------------------------------------------------------------
def gen_dets():
    """The 36 M_s=0 determinants of CAS(4,4) as sorted spin-orbital tuples.

    Spin orbital index = orbital (alpha) or orbital + NORB (beta).  Ordering is
    alpha-pair outer, beta-pair inner, both ascending -- identical to the engine
    and to the Fortran ``gen_dets``.
    """
    dets = []
    for a1 in range(NORB):
        for a2 in range(a1 + 1, NORB):
            for b1 in range(NORB):
                for b2 in range(b1 + 1, NORB):
                    dets.append(tuple(sorted([a1, a2, b1 + NORB, b2 + NORB])))
    return dets


def det_config(det):
    """Spatial occupation string ('0','u','d','2' per orbital) of a determinant."""
    s = []
    for o in range(NORB):
        a, b = o in det, (o + NORB) in det
        s.append("2" if (a and b) else "u" if a else "d" if b else "0")
    return "".join(s)


def os_class(cfg):
    """Number of singly-occupied orbitals: 0 (0OS), 2 (2OS), or 4 (4OS)."""
    return sum(c in "ud" for c in cfg)


def fermi_phase(cfg):
    """(-1)^nu: parity of the permutation taking the orbital-ordered creation
    string (orbital by orbital, alpha before beta) into ascending spin-orbital
    order.  This is the per-determinant element of the diagonal phase D."""
    seq = []
    for o, ch in enumerate(cfg):
        if ch == "2":
            seq += [o, o + NORB]
        elif ch == "u":
            seq += [o]
        elif ch == "d":
            seq += [o + NORB]
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


# ---------------------------------------------------------------------------
# Slater-Condon CAS(4,4) M_s=0 Hamiltonian and S^2
# ---------------------------------------------------------------------------
def _spinorb_integrals(h_act, eri_act, coulomb_only=False):
    """Spin-orbital one-electron H1 and antisymmetrized two-electron g.

    ``eri_act`` is chemist ``(ij|kl)``; the physicist antisymmetrized element is
    g(P,Q,R,S) = <PQ|RS> - <PQ|SR> = (PR|QS) - (PS|QR) with spin selection.

    ``coulomb_only=True`` drops the exchange term ``-(PS|QR)`` of every coupling,
    leaving only the direct/Hartree part ``(PR|QS)``.  This is a ROLE-based split
    (it removes exchange by its role in each antisymmetrized element, not by
    integral value), so it is unambiguous for 4-distinct-orbital integrals."""
    spat = [p % NORB for p in range(NSO)]
    spin = [p // NORB for p in range(NSO)]
    H1 = np.zeros((NSO, NSO))
    for P in range(NSO):
        for Q in range(NSO):
            if spin[P] == spin[Q]:
                H1[P, Q] = h_act[spat[P], spat[Q]]
    g = np.zeros((NSO, NSO, NSO, NSO))
    for P in range(NSO):
        for Q in range(NSO):
            for R in range(NSO):
                for S in range(NSO):
                    a = eri_act[spat[P], spat[R], spat[Q], spat[S]] \
                        if (spin[P] == spin[R] and spin[Q] == spin[S]) else 0.0
                    b = eri_act[spat[P], spat[S], spat[Q], spat[R]] \
                        if (spin[P] == spin[S] and spin[Q] == spin[R]) else 0.0
                    g[P, Q, R, S] = a if coulomb_only else (a - b)
    return H1, g


def _melem(D1, D2, H1, g):
    """Slater-Condon matrix element <D1|H|D2> for sorted spin-orbital dets."""
    D1, D2 = list(D1), list(D2)
    holes = [x for x in D2 if x not in D1]
    parts = [x for x in D1 if x not in D2]
    common = [x for x in D1 if x in D2]
    nh = len(holes)
    if nh > 2:
        return 0.0
    occ, sgn = list(D2), 1.0
    for hole in holes:
        idx = occ.index(hole)
        if idx % 2 == 1:
            sgn = -sgn
        occ.pop(idx)
    for p in reversed(parts):
        idx = sum(1 for x in occ if x < p)
        if idx % 2 == 1:
            sgn = -sgn
        occ.insert(idx, p)
    if nh == 0:
        e = sum(H1[i, i] for i in D1)
        for i in range(4):
            for k in range(i + 1, 4):
                e += g[D1[i], D1[k], D1[i], D1[k]]
        return e
    if nh == 1:
        Pp, Hh = parts[0], holes[0]
        val = H1[Pp, Hh] + sum(g[Pp, Qc, Hh, Qc] for Qc in common)
        return sgn * val
    return sgn * g[parts[0], parts[1], holes[0], holes[1]]


def _apply_ex(det, cre, ann):
    """a^dag_cre a_ann on a sorted spin-orbital det -> (new det, phase) or (None,0)."""
    if ann not in det:
        return None, 0.0
    if cre != ann and cre in det:
        return None, 0.0
    occ, ph = list(det), 1.0
    idx = occ.index(ann)
    if idx % 2 == 1:
        ph = -ph
    occ.pop(idx)
    idx = sum(1 for x in occ if x < cre)
    if idx % 2 == 1:
        ph = -ph
    occ.insert(idx, cre)
    return tuple(occ), ph


def build_H_S2(h_act, eri_act, dets=None):
    """Return (H, S2): the 36x36 Slater-Condon CAS(4,4) M_s=0 Hamiltonian and the
    determinant-basis S^2 = S_- S_+ (M_s=0 sector)."""
    if dets is None:
        dets = gen_dets()
    H1, g = _spinorb_integrals(h_act, eri_act)
    H = np.array([[_melem(dets[i], dets[j], H1, g) for j in range(NDET)]
                  for i in range(NDET)])
    index = {d: i for i, d in enumerate(dets)}
    S2 = np.zeros((NDET, NDET))
    for j in range(NDET):
        for p in range(NORB):
            Di, ph1 = _apply_ex(dets[j], p, p + NORB)       # S_+
            if ph1 == 0.0:
                continue
            for q in range(NORB):
                Dk, ph2 = _apply_ex(Di, q + NORB, q)        # S_-
                if ph2 == 0.0:
                    continue
                if Dk in index:
                    S2[index[Dk], j] += ph1 * ph2
    return H, S2


def build_H_coulomb(h_act, eri_act, dets=None):
    """The Coulomb-only (direct, no exchange) 36-det Slater-Condon H: the
    antisymmetrized exchange term ``-(PS|QR)`` is dropped from every coupling.

    Used for the ROLE-based c_H exchange seam: the exchange part of the full
    Hamiltonian is ``H_exchange = build_H_S2(...)[0] - build_H_coulomb(...)`` (the
    one-electron H1 and the Hartree ``(PR|QS)`` parts cancel, leaving exactly the
    ``-(PS|QR)`` exchange).  Because the split is by ROLE (not integral value) it
    is unambiguous for the 4-distinct-orbital ``(pq|rs)`` integrals."""
    if dets is None:
        dets = gen_dets()
    H1, g = _spinorb_integrals(h_act, eri_act, coulomb_only=True)
    return np.array([[_melem(dets[i], dets[j], H1, g) for j in range(NDET)]
                     for i in range(NDET)])


# ---------------------------------------------------------------------------
# the HARD-CODED spin-adaptation U (rows = CSFs, columns = engine dets)
# ---------------------------------------------------------------------------
def build_U(dets=None):
    """Return (U, Sval, labels): the 36x36 hard-coded spin-adaptation.

    ``U`` rows are the CSFs expressed in the engine sorted-determinant basis,
    ordered 20 singlet, then 15 triplet, then 1 quintet.  ``Sval[k]`` is the
    total spin S of CSF k (0/1/2); ``labels[k]`` a human-readable tag.

    U = U_orbital * D, where U_orbital are the orbital-ordered CSFs (0OS
    identity, 2OS Hadamard, 4OS genealogical) and D is the diagonal Fermi phase
    -- NO numerical S^2 diagonalization is performed.
    """
    if dets is None:
        dets = gen_dets()
    cfg = [det_config(d) for d in dets]
    osc = [os_class(c) for c in cfg]
    cfg_index = {c: i for i, c in enumerate(cfg)}
    Dphase = np.array([fermi_phase(c) for c in cfg], float)

    sing, trip, quint = [], [], []   # each entry: (vector-over-dets, label)

    def vec(coeffs):
        v = np.zeros(NDET)
        for d, c in coeffs:
            v[d] = c
        return v

    # 0OS : the six closed-shell determinants are singlets (identity).
    for i in range(NDET):
        if osc[i] == 0:
            sing.append((vec([(i, 1.0)]), f"0OS:{cfg[i]}"))

    # 2OS : Hadamard per left/right pair (partner = u<->d swap), orbital-ordered.
    swap = str.maketrans("ud", "du")
    seen = set()
    for i in range(NDET):
        if osc[i] != 2 or i in seen:
            continue
        j = cfg_index[cfg[i].translate(swap)]
        seen |= {i, j}
        sing.append((vec([(i, 1 / SQ2), (j, -1 / SQ2)]), f"2OS-S:{cfg[i]}/{cfg[j]}"))
        trip.append((vec([(i, 1 / SQ2), (j,  1 / SQ2)]), f"2OS-T:{cfg[i]}/{cfg[j]}"))

    # 4OS : hard-coded genealogical CSFs over FOUR_OS_ORDER.
    fo = [cfg_index[c] for c in FOUR_OS_ORDER]
    for name, S, coeff in FOUR_OS_CSF:
        v = vec(list(zip(fo, coeff)))
        (sing if S == 0 else trip if S == 1 else quint).append((v, name))

    blocks = [(sing, 0), (trip, 1), (quint, 2)]
    rows, Sval, labels = [], [], []
    for entries, S in blocks:
        for v, lab in entries:
            rows.append(v)
            Sval.append(S)
            labels.append(lab)

    U_orb = np.array(rows)                 # orbital-ordered CSFs (rows)
    U = U_orb * Dphase[None, :]            # dress with Fermi phase -> engine basis
    return U, np.array(Sval), labels


# ---------------------------------------------------------------------------
# assembled QMRSF-DK blocks
# ---------------------------------------------------------------------------
def build_dk_matrices(h_act, eri_act, e_ref=0.0):
    """Build the three spin-pure QMRSF-DK matrices for given CAS(4,4) integrals.

    Parameters
    ----------
    h_act   : (4,4) one-electron active integrals.
    eri_act : (4,4,4,4) chemist (ij|kl) active electron-repulsion integrals.
    e_ref   : optional scalar subtracted from every diagonal (reference shift).

    Returns a dict with the spin-pure blocks ``A_singlet`` (20x20),
    ``A_triplet`` (15x15), ``A_quintet`` (1x1), plus the spin-adaptation ``U``,
    the determinant-basis ``H`` and ``S2``, the per-CSF spin ``Sval`` and the
    ``labels``.  The blocks are A = U H U^T partitioned by total spin.
    """
    h_act = np.asarray(h_act, float)
    eri_act = np.asarray(eri_act, float)
    dets = gen_dets()
    H, S2 = build_H_S2(h_act, eri_act, dets)
    U, Sval, labels = build_U(dets)

    Hcsf = U @ H @ U.T
    if e_ref:
        Hcsf = Hcsf - e_ref * np.eye(NDET)

    sel = {S: np.where(Sval == S)[0] for S in (0, 1, 2)}
    A = {S: Hcsf[np.ix_(sel[S], sel[S])] for S in (0, 1, 2)}
    return {
        "A_singlet": A[0],            # 20x20
        "A_triplet": A[1],            # 15x15
        "A_quintet": A[2],            #  1x1
        "U": U, "H": H, "S2": S2,
        "Sval": Sval, "labels": labels,
    }


# ---------------------------------------------------------------------------
# random integrals + self-test
# ---------------------------------------------------------------------------
def random_integrals(seed=0):
    """Random symmetric h_act and 8-fold-symmetric chemist eri_act (for tests)."""
    rng = np.random.default_rng(seed)
    h = rng.standard_normal((NORB, NORB))
    h = 0.5 * (h + h.T)
    e = rng.standard_normal((NORB, NORB, NORB, NORB))
    # impose chemist 8-fold permutational symmetry (ij|kl)
    e = 0.5 * (e + e.transpose(1, 0, 2, 3))
    e = 0.5 * (e + e.transpose(0, 1, 3, 2))
    e = 0.5 * (e + e.transpose(2, 3, 0, 1))
    return h, e


def self_test(seed=0, verbose=True):
    """Verify orthonormality, spin block-diagonality, and CI-spectrum match on
    random integrals.  Returns a dict of error metrics."""
    h, e = random_integrals(seed)
    out = build_dk_matrices(h, e)
    U, H, S2, Sval = out["U"], out["H"], out["S2"], out["Sval"]

    res = {}
    res["orthonormal"] = float(np.abs(U @ U.T - np.eye(NDET)).max())

    # <S^2> of each CSF == S(S+1)
    s2diag = np.diag(U @ S2 @ U.T)
    res["s2_err"] = float(max(np.abs(s2diag[Sval == S] - S * (S + 1)).max()
                              for S in (0, 1, 2)))

    # cross-spin-block coupling of U H U^T must vanish
    Hcsf = U @ H @ U.T
    cross = 0.0
    for i in range(NDET):
        for j in range(NDET):
            if Sval[i] != Sval[j]:
                cross = max(cross, abs(Hcsf[i, j]))
    res["cross_block"] = float(cross)

    # block dims
    res["dims"] = (int((Sval == 0).sum()), int((Sval == 1).sum()),
                   int((Sval == 2).sum()))

    # each block spectrum == S^2-projected H spectrum
    w, V = np.linalg.eigh(H)
    s2f = np.array([V[:, k] @ S2 @ V[:, k] for k in range(NDET)])
    blocks = {0: out["A_singlet"], 1: out["A_triplet"], 2: out["A_quintet"]}
    eigerr = 0.0
    for S in (0, 1, 2):
        eb = np.sort(np.linalg.eigvalsh(np.atleast_2d(blocks[S])))
        ed = np.sort(w[np.abs(s2f - S * (S + 1)) < 1e-6])
        eigerr = max(eigerr, float(np.abs(eb - ed).max()))
    res["block_spectrum_vs_CI"] = eigerr

    ok = (res["orthonormal"] < 1e-10 and res["s2_err"] < 1e-9 and
          res["cross_block"] < 1e-10 and res["dims"] == (20, 15, 1) and
          res["block_spectrum_vs_CI"] < 1e-9)
    res["PASS"] = bool(ok)

    if verbose:
        print("==== QMRSF-DK hard-coded spin-adapted matrices (NumPy) ====")
        print(f"  U orthonormal err            = {res['orthonormal']:.3e}")
        print(f"  <S^2> err (0/2/6)            = {res['s2_err']:.3e}")
        print(f"  block dims (S/T/Q)           = {res['dims']}")
        print(f"  max cross-spin-block |UHU^T| = {res['cross_block']:.3e}")
        print(f"  block spectra vs S^2-proj CI = {res['block_spectrum_vs_CI']:.3e}")
        print(f"  RESULT: {'PASS' if ok else 'FAIL'}")
    return res


def crosscheck_reference(npy_path, labels_pkl=None):
    """Optional: confirm the hard-coded U reproduces the saved verified reference
    U_explicit_36.npy (column-CSF convention) bit-for-bit, per spin block, up to
    per-CSF sign.  Returns the worst per-block row-set mismatch."""
    import pickle
    Uref = np.load(npy_path)
    U, Sval, _ = build_U()
    ref_rows = Uref.T  # reference stores CSFs as columns
    worst = 0.0
    start = 0
    for S, n in ((0, 20), (1, 15), (2, 1)):
        mine = U[Sval == S]
        pool = list(ref_rows[start:start + n])
        for r in mine:
            ds = [min(np.abs(rr - r).max(), np.abs(rr + r).max()) for rr in pool]
            k = int(np.argmin(ds))
            worst = max(worst, ds[k])
            pool.pop(k)
        start += n
    return float(worst)


if __name__ == "__main__":
    res = self_test(seed=0)
    # a second random seed for robustness
    res2 = self_test(seed=12345, verbose=False)
    assert res["PASS"] and res2["PASS"], "self-test FAILED"
    print(f"  second random seed PASS       = {res2['PASS']}")

    # optional cross-check against the saved verified reference, if present
    import os
    ref = ("/tmp/claude-1935200050/-bighome-alireza/"
           "e30f4aec-794a-4be6-abdb-c9ad8b6a5ba2/scratchpad/U_explicit_36.npy")
    if os.path.exists(ref):
        w = crosscheck_reference(ref)
        print(f"  hard-coded U vs U_explicit_36 = {w:.3e}  (per-block, up to sign)")
        assert w < 1e-12, "U does not match the verified reference"
    print("ALL CHECKS PASSED")
