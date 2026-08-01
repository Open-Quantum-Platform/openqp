"""Direct (matrix-free) determinant QDPT2 engine for the GAMESS-convention
family (MRMP2 / MCQDPT2 / XMCQDPT2).

The dense engine in :mod:`oqp.library.caspt2_dyall` assembles the full
Hamiltonian and H0 over EVERY determinant of the correlated orbital space and
diagonalizes the external block -- O(N_det^2) memory and O(N_det^3) work, which
caps it at toy systems (``[pt2] max_det``).  The QDPT family's zeroth order is
``H0 = sum_p eps_p n_p``: DIAGONAL in the determinant basis.  That removes any
need for the external eigenbasis, and the first-order interacting space is only
the singles/doubles generated from the N_CAS reference determinants:

    V_I(Phi)  = <Phi|H|Psi_I> = sum_D c_I(D) <Phi|H|D>,   Phi in Q
    E0(Phi)   = sum_{p occ(Phi)} eps_p                      (diagonal resolvent)
    Heff[I,J] = <I|H|J> + 1/2 sum_Phi [ V_I R_J V_J + V_J R_I V_I ](Phi)

with R_J(Phi) the (level/imaginary/ISA-regularized) 1/(E0_J - E0(Phi)).  This
module therefore streams excitation CLASSES (alpha/beta singles; alpha-alpha,
beta-beta, alpha-beta doubles) from each reference determinant with fully
NumPy-vectorized integral gathers, phases (``np.bitwise_count`` prefix masks)
and diagonal denominators, merges duplicates with a lexicographic unique over
(alpha-word, beta-word) keys, and never materializes any N x N matrix.

Scaling: work ~ N_CAS * (n_occ^2 n_virt^2) streamed in per-reference batches;
memory ~ the merged external space (guarded by ``[pt2] max_terms``).  Numerical
conventions (integrals, phases, shifts) are pinned to the dense engine by the
equivalence tests in ``tests/test_qdpt2_direct.py`` (agreement <= 1e-10 Eh).

Scope: the diagonal-Fock QDPT family only.  The CASPT2 Fock-MATRIX H0 and the
Dyall H0 are not diagonal in the determinant basis and keep the dense path.
Alpha/beta strings are encoded in uint64 words (norb <= 63 per spin).
"""
from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------- bit utilities
def _prefix_phase_mask(lo, hi):
    """uint64 mask of bits strictly between positions lo < hi (exclusive)."""
    return ((np.uint64(1) << hi) - np.uint64(1)) ^ ((np.uint64(1) << (lo + np.uint64(1))) - np.uint64(1))


def _single_phases(word, occ_idx, virt_idx):
    """Phases for all (i -> a) single excitations within one spin word.

    Returns ``phase[n_occ, n_virt]`` with ``(-1)**n_between`` where n_between is
    the number of set bits strictly between i and a in the SOURCE word -- the
    net crossing of an annihilation at i followed by a creation at a."""
    occ = occ_idx.astype(np.uint64)[:, None]
    vrt = virt_idx.astype(np.uint64)[None, :]
    lo = np.minimum(occ, vrt)
    hi = np.maximum(occ, vrt)
    mask = _prefix_phase_mask(lo, hi)
    nbet = np.bitwise_count(np.uint64(word) & mask)
    return np.where(nbet % 2 == 0, 1.0, -1.0)


def _seq_phase_ann(word, pos):
    """Phase of a_pos acting on word (bits below pos), and the new word."""
    below = np.bitwise_count(np.uint64(word) & ((np.uint64(1) << pos) - np.uint64(1)))
    return np.where(below % 2 == 0, 1.0, -1.0), np.uint64(word) ^ (np.uint64(1) << pos)


def _seq_phase_cre(word, pos):
    below = np.bitwise_count(np.uint64(word) & ((np.uint64(1) << pos) - np.uint64(1)))
    return np.where(below % 2 == 0, 1.0, -1.0), np.uint64(word) | (np.uint64(1) << pos)


def _occ_virt(word, norb):
    bits = (int(word) >> np.arange(norb)) & 1
    occ = np.nonzero(bits)[0]
    virt = np.nonzero(1 - bits)[0]
    return occ.astype(np.int64), virt.astype(np.int64)


# --------------------------------------------------------------------------- per-reference streaming
class _Stream:
    """Accumulates (alpha_key, beta_key, element, e0) batches across classes."""

    def __init__(self):
        self.ka, self.kb, self.val, self.e0, self.src = [], [], [], [], []

    def push(self, ka, kb, val, e0, src_index):
        n = val.size
        if n == 0:
            return
        self.ka.append(np.broadcast_to(ka, (n,)).astype(np.uint64).ravel())
        self.kb.append(np.broadcast_to(kb, (n,)).astype(np.uint64).ravel())
        self.val.append(val.ravel().astype(float))
        self.e0.append(e0.ravel().astype(float))
        self.src.append(np.full(n, src_index, dtype=np.int32))

    def arrays(self):
        if not self.val:
            z = np.zeros(0)
            return z.astype(np.uint64), z.astype(np.uint64), z, z, z.astype(np.int32)
        return (np.concatenate(self.ka), np.concatenate(self.kb),
                np.concatenate(self.val), np.concatenate(self.e0),
                np.concatenate(self.src))


def _mean_field_single(h1e, eri, occ_same, occ_other):
    """Effective one-electron matrix for a single excitation within one spin:
    F_eff[a, i] = h[a, i] + sum_{k in same-spin occ} [(ai|kk) - (ak|ki)]
                          + sum_{k in other-spin occ} (ai|kk)
    The k == i term belongs to the same-spin sum and is removed by the caller
    columnwise (i is occupied, so it appears in occ_same)."""
    F = h1e.copy()
    if occ_same.size:
        F = F + eri[:, :, occ_same, occ_same].sum(axis=2)   # + (pq|kk)
        F = F - eri[:, occ_same, occ_same, :].sum(axis=1)   # - (pk|kq)
    if occ_other.size:
        F = F + eri[:, :, occ_other, occ_other].sum(axis=2)
    return F


def _stream_reference(stream, src_index, wa, wb, h1e, eri, eps, norb,
                      internal_test):
    """Push every external single/double excitation of one reference
    determinant (alpha word ``wa``, beta word ``wb``) onto the stream."""
    occ_a, virt_a = _occ_virt(wa, norb)
    occ_b, virt_b = _occ_virt(wb, norb)
    e0_ref = float(eps[occ_a].sum() + eps[occ_b].sum())
    one_a = np.uint64(1)

    # ---- singles (alpha), with the full mean-field single-excitation element
    if occ_a.size and virt_a.size:
        Fa = _mean_field_single(h1e, eri, occ_a, occ_b)
        # remove the k == i self term columnwise: -( (ai|ii) - (ai|ii) ) = 0 for
        # the same-spin pair, so no correction is needed beyond using occ sets
        # of the SOURCE determinant (i in occ includes k = i whose direct and
        # exchange contributions cancel exactly).
        elem = Fa[np.ix_(virt_a, occ_a)].T                      # [i, a]
        phase = _single_phases(wa, occ_a, virt_a)               # [i, a]
        ka = (np.uint64(wa) ^ (one_a << occ_a.astype(np.uint64))[:, None]) \
            | (one_a << virt_a.astype(np.uint64))[None, :]
        e0 = e0_ref - eps[occ_a][:, None] + eps[virt_a][None, :]
        keep = ~internal_test(ka, np.uint64(wb))
        stream.push(ka[keep], np.uint64(wb), (elem * phase)[keep], e0[keep], src_index)

    # ---- singles (beta)
    if occ_b.size and virt_b.size:
        Fb = _mean_field_single(h1e, eri, occ_b, occ_a)
        elem = Fb[np.ix_(virt_b, occ_b)].T
        phase = _single_phases(wb, occ_b, virt_b)
        kb = (np.uint64(wb) ^ (one_a << occ_b.astype(np.uint64))[:, None]) \
            | (one_a << virt_b.astype(np.uint64))[None, :]
        e0 = e0_ref - eps[occ_b][:, None] + eps[virt_b][None, :]
        keep = ~internal_test(np.uint64(wa), kb)
        stream.push(np.uint64(wa), kb[keep], (elem * phase)[keep], e0[keep], src_index)

    # ---- same-spin doubles
    for word, occ, virt, is_alpha in ((wa, occ_a, virt_a, True),
                                      (wb, occ_b, virt_b, False)):
        if occ.size < 2 or virt.size < 2:
            continue
        oi, oj = np.triu_indices(occ.size, k=1)
        ai, bj = np.triu_indices(virt.size, k=1)
        i = occ[oi]; j = occ[oj]                                # i < j
        a = virt[ai]; b = virt[bj]                              # a < b
        # element: <ij||ab> in chemist integrals = (ai|bj) - (aj|bi)
        AI = eri[a[:, None], i[None, :], b[:, None], j[None, :]]   # [vp, op]
        AJ = eri[a[:, None], j[None, :], b[:, None], i[None, :]]
        elem = (AI - AJ).T                                          # [op, vp]
        # sequential phases: ann i, ann j, cre b, cre a on the same word
        w0 = np.uint64(word)
        p1, w1 = _seq_phase_ann(w0, i.astype(np.uint64))
        p2, w2 = _seq_phase_ann(w1, j.astype(np.uint64))
        p3, w3 = _seq_phase_cre(w2[:, None], b.astype(np.uint64)[None, :])
        p4, w4 = _seq_phase_cre(w3, a.astype(np.uint64)[None, :])
        phase = (p1 * p2)[:, None] * p3 * p4                        # [op, vp]
        e0 = (e0_ref - eps[i] - eps[j])[:, None] + (eps[a] + eps[b])[None, :]
        if is_alpha:
            keep = ~internal_test(w4, np.uint64(wb))
            stream.push(w4[keep], np.uint64(wb), (elem * phase)[keep], e0[keep], src_index)
        else:
            keep = ~internal_test(np.uint64(wa), w4)
            stream.push(np.uint64(wa), w4[keep], (elem * phase)[keep], e0[keep], src_index)

    # ---- alpha-beta doubles
    if occ_a.size and virt_a.size and occ_b.size and virt_b.size:
        # element (ai|bj); phases separate per spin word
        elem = eri[np.ix_(virt_a, occ_a, virt_b, occ_b)]        # [aA, iA, bB, jB]
        pa = _single_phases(wa, occ_a, virt_a).T                # [aA, iA]
        pb = _single_phases(wb, occ_b, virt_b).T                # [bB, jB]
        ka = (np.uint64(wa) ^ (one_a << occ_a.astype(np.uint64))[None, :]) \
            | (one_a << virt_a.astype(np.uint64))[:, None]      # [aA, iA]
        kb = (np.uint64(wb) ^ (one_a << occ_b.astype(np.uint64))[None, :]) \
            | (one_a << virt_b.astype(np.uint64))[:, None]      # [bB, jB]
        e0 = (e0_ref
              + (eps[virt_a][:, None] - eps[occ_a][None, :])[:, :, None, None]
              + (eps[virt_b][:, None] - eps[occ_b][None, :])[None, None, :, :])
        val = elem * pa[:, :, None, None] * pb[None, None, :, :]
        KA = np.broadcast_to(ka[:, :, None, None], val.shape)
        KB = np.broadcast_to(kb[None, None, :, :], val.shape)
        keep = ~internal_test(KA, KB)
        stream.push(KA[keep], KB[keep], val[keep], e0[keep], src_index)


# --------------------------------------------------------------------------- parallel streaming
def _make_internal_test(ncore, nact, norb):
    core_mask = np.uint64((1 << ncore) - 1)
    virt_mask = np.uint64(((1 << (norb - ncore - nact)) - 1) << (ncore + nact))

    def internal_test(ka, kb):
        ka = np.asarray(ka, dtype=np.uint64)
        kb = np.asarray(kb, dtype=np.uint64)
        return (((ka & core_mask) == core_mask) & ((kb & core_mask) == core_mask)
                & ((ka & virt_mask) == np.uint64(0))
                & ((kb & virt_mask) == np.uint64(0)))

    return internal_test


def _stream_chunk(payload):
    """Worker entry: stream one chunk of reference determinants.

    Top-level (spawn-picklable).  Each worker re-imports numpy only; the
    integral arrays ride in via the payload (one-time cost per chunk), so
    ``[pt2] nproc > 1`` pays off when N_CAS is large, not for toy references."""
    (h1e, eri, eps, norb, ncore, nact, sup_a, sup_b, base) = payload
    internal_test = _make_internal_test(ncore, nact, norb)
    stream = _Stream()
    for k, (wa, wb) in enumerate(zip(sup_a, sup_b)):
        _stream_reference(stream, base + k, wa, wb, h1e, eri, eps, norb,
                          internal_test)
    return stream.arrays()


def _stream_all(h1e, eri, eps, norb, ncore, nact, sup_a, sup_b, nproc):
    nsup = sup_a.size
    if nproc <= 1 or nsup < 4 * nproc:
        return _stream_chunk((h1e, eri, eps, norb, ncore, nact,
                              sup_a, sup_b, 0))
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    bounds = np.linspace(0, nsup, nproc + 1).astype(int)
    payloads = [(h1e, eri, eps, norb, ncore, nact,
                 sup_a[b0:b1], sup_b[b0:b1], int(b0))
                for b0, b1 in zip(bounds[:-1], bounds[1:]) if b1 > b0]
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=nproc, mp_context=ctx) as pool:
        parts = list(pool.map(_stream_chunk, payloads))
    return tuple(np.concatenate([p[i] for p in parts]) for i in range(5))


# --------------------------------------------------------------------------- Fortran kernel bridge
def _terms_per_reference(norb, na_occ, nb_occ):
    """Exact streamed-term count per reference (identical for every reference,
    since all share the electron counts)."""
    nva, nvb = norb - na_occ, norb - nb_occ
    singles = na_occ * nva + nb_occ * nvb
    dss = ((na_occ * (na_occ - 1)) // 2) * ((nva * (nva - 1)) // 2) \
        + ((nb_occ * (nb_occ - 1)) // 2) * ((nvb * (nvb - 1)) // 2)
    dab = na_occ * nva * nb_occ * nvb
    return singles + dss + dab


def _stream_fortran(h1e, eri, eps, norb, ncore, nact, sup_a, sup_b, C,
                    max_terms, nthreads):
    """Run the liboqp OpenMP streaming kernel; returns (ka, kb, e0, V) with the
    external space already merged/unique, or None when the kernel is absent."""
    try:
        from oqp import lib, ffi
    except Exception:
        return None
    if not hasattr(lib, "qdpt2_stream_kernel"):
        return None

    nsup, nstate = int(sup_a.size), int(C.shape[1])
    na_occ = int(np.bitwise_count(sup_a[0]))
    nb_occ = int(np.bitwise_count(sup_b[0]))
    total = nsup * _terms_per_reference(norb, na_occ, nb_occ)
    if total > max_terms:
        raise ValueError(
            f"direct QDPT2 stream would produce {total} terms > [pt2] "
            f"max_terms={max_terms}; raise the guard or shrink the space")
    # sort/reduce kernel: uniques <= streamed terms, so the exact term count
    # is a tight output bound (no power-of-two hash headroom needed)
    cap = max(64, int(total))

    h1e_c = np.ascontiguousarray(h1e, dtype=np.float64)
    eri_c = np.ascontiguousarray(eri, dtype=np.float64)
    eps_c = np.ascontiguousarray(eps, dtype=np.float64)
    cvec = np.ascontiguousarray(C, dtype=np.float64)          # [nsup, nstate]
    ka_in = np.ascontiguousarray(sup_a).view(np.int64)
    kb_in = np.ascontiguousarray(sup_b).view(np.int64)
    out_ka = np.zeros(cap, dtype=np.int64)
    out_kb = np.zeros(cap, dtype=np.int64)
    out_e0 = np.zeros(cap, dtype=np.float64)
    out_v = np.zeros((cap, nstate), dtype=np.float64)         # slot-major

    def ptr(arr, ctype):
        return ffi.cast(ctype, arr.ctypes.data)

    n = int(lib.qdpt2_stream_kernel(
        int(norb), int(ncore), int(nact), nsup, nstate, int(max(1, nthreads)),
        cap,
        ptr(ka_in, "int64_t *"), ptr(kb_in, "int64_t *"),
        ptr(cvec, "double *"), ptr(h1e_c, "double *"),
        ptr(eri_c, "double *"), ptr(eps_c, "double *"),
        ptr(out_ka, "int64_t *"), ptr(out_kb, "int64_t *"),
        ptr(out_e0, "double *"), ptr(out_v, "double *")))
    if n < 0:
        raise RuntimeError(
            "QDPT2 fortran kernel hash table overflowed; raise [pt2] max_terms "
            "(capacity follows it) or fall back to [pt2] engine=direct")
    return (out_ka[:n].view(np.uint64), out_kb[:n].view(np.uint64),
            out_e0[:n].copy(), np.ascontiguousarray(out_v[:n].T))


# --------------------------------------------------------------------------- assembly
def _regularized_inverse(d, options):
    """1/d with the same regularizations as the dense engine's _first_order."""
    d = d + options.level_shift
    if options.edshft:
        with np.errstate(divide="ignore", invalid="ignore"):
            d = d + options.edshft / d
    if options.imaginary_shift:
        return d / (d * d + options.imaginary_shift ** 2)
    return np.divide(1.0, d, out=np.zeros_like(d),
                     where=(np.abs(d) > 1.0e-300) & np.isfinite(d))


def direct_qdpt2(h1e, eri, coeffs, energies, dets, eps, D_sa, ncore, nact,
                 active_nelec, norb, enuc, roots, options, rotate):
    """Matrix-free diagonal-H0 QDPT2 over the reference roots.

    Inputs mirror :func:`oqp.library.caspt2_dyall._multistate` (semicanonical,
    frozen-core-reduced quantities).  Returns the same result dictionary."""
    if norb > 63:
        raise ValueError("direct QDPT2 supports up to 63 orbitals per spin "
                         f"(got norb={norb}); reduce the basis or freeze more")
    nstate = len(roots)
    eps = np.asarray(eps, dtype=float)

    # ---- reference determinants (support) as spin words in the reduced space
    core_mask = np.uint64((1 << ncore) - 1)
    act_mask = (1 << nact) - 1

    sup_a = np.array([int(core_mask) | ((d & act_mask) << ncore)
                      for d in dets], dtype=np.uint64)
    sup_b = np.array([int(core_mask) | ((int(d) >> nact) << ncore)
                      for d in dets], dtype=np.uint64)

    # ---- reference CI matrix over support (columns = requested roots)
    C = np.asarray(coeffs, dtype=float)[:, list(roots)].copy()
    E_ref_in = np.asarray([float(energies[r]) for r in roots])

    # ---- XMS rotation: diagonalize the SA Fock in the model space.
    # F is a one-electron operator; its model matrix only needs the support:
    # <D'|F|D> = delta * sum_occ F_pp + single-excitation elements.
    rotated = bool(rotate)
    if rotated:
        from oqp.library.caspt2_dyall import _effective_fock
        F = _effective_fock(h1e, eri, D_sa)
        fmodel = _support_one_electron_model(F, sup_a, sup_b, norb, C)
        _w, R = np.linalg.eigh(0.5 * (fmodel + fmodel.T))
        C = C @ R
        hmodel = R.T @ np.diag(E_ref_in) @ R
    else:
        hmodel = np.diag(E_ref_in)

    # ---- zeroth-order energies of the (possibly rotated) references
    e0_sup = np.array([float(eps[_occ_virt(a, norb)[0]].sum()
                             + eps[_occ_virt(b, norb)[0]].sum())
                       for a, b in zip(sup_a, sup_b)])
    e0_state = np.array([float((C[:, i] ** 2) @ e0_sup) for i in range(nstate)])

    # ---- stream all classes from every support determinant.
    # Engine history (114.5M-term CAS(8,8)-in-24 benchmark): the v1 hash
    # kernel was DRAM-bound (73 s); the v2 sort/reduce kernel measures 18 s
    # serial and 7.4 s at 8 threads once explicit thread requests stopped
    # being clamped by the core's ambient OMP setting -- versus ~9 s for the
    # NumPy lexsort path.  'auto' therefore prefers the kernel with a
    # multi-thread default ([pt2] nproc overrides); 'direct' forces the NumPy
    # path; both are equivalence-pinned.  The serial T-way merge still caps
    # thread scaling at ~2.4x -- a key-range-parallel merge is the follow-up.
    nproc = int(getattr(options, "nproc", 1))
    max_terms = int(getattr(options, "max_terms", 30_000_000))
    engine = str(getattr(options, "engine", "auto"))
    merged = None
    if engine in {"auto", "fortran"}:
        import os
        kernel_threads = nproc if nproc > 1 else max(1, min(8, os.cpu_count() or 1))
        merged = _stream_fortran(h1e, eri, eps, norb, ncore, nact,
                                 sup_a, sup_b, C, max_terms, kernel_threads)
        if merged is None and engine == "fortran":
            raise RuntimeError(
                "[pt2] engine=fortran requested but liboqp lacks "
                "qdpt2_stream_kernel; rebuild the core or use engine=auto")

    if merged is not None:
        _ka_u, _kb_u, e0_ext, V = merged
        nuniq = int(e0_ext.size)
    else:
        ka, kb, val, e0, src = _stream_all(h1e, eri, eps, norb, ncore, nact,
                                           sup_a, sup_b, nproc)
        if val.size > max_terms:
            raise ValueError(
                f"direct QDPT2 stream produced {val.size} terms > [pt2] "
                f"max_terms={max_terms}; raise the guard or shrink the space")

        # merge duplicate external determinants
        order = np.lexsort((kb, ka))
        ka, kb, val, e0, src = ka[order], kb[order], val[order], e0[order], src[order]
        new = np.empty(val.size, dtype=bool)
        if val.size:
            new[0] = True
            new[1:] = (ka[1:] != ka[:-1]) | (kb[1:] != kb[:-1])
        uidx = np.cumsum(new) - 1
        nuniq = int(uidx[-1]) + 1 if val.size else 0
        e0_ext = np.full(nuniq, np.inf)
        if val.size:
            np.minimum.at(e0_ext, uidx, e0)   # duplicates carry identical E0

        V = np.zeros((nstate, nuniq))
        for i in range(nstate):
            V[i] = np.bincount(uidx, weights=val * C[src, i], minlength=nuniq)

    # ---- diagonal resolvent, E2, effective Hamiltonian
    e2 = np.zeros(nstate)
    min_denoms = np.full(nstate, float("inf"))
    RV = np.zeros_like(V)
    for i in range(nstate):
        d = e0_ext - e0_state[i]
        if d.size:
            min_denoms[i] = float(np.min(np.abs(d)))
        RV[i] = -V[i] * _regularized_inverse(d, options)
        e2[i] = float(V[i] @ RV[i])

    heff = hmodel.copy()
    for i in range(nstate):
        for j in range(nstate):
            heff[i, j] += 0.5 * (V[i] @ RV[j] + V[j] @ RV[i])
    heff = 0.5 * (heff + heff.T)
    ms_energies, mixing = np.linalg.eigh(heff)
    ss_energies = np.array([hmodel[i, i] + e2[i] for i in range(nstate)])

    return {
        "ms_energies": ms_energies,
        "mixing": mixing,
        "heff": heff,
        "ss_energies": ss_energies,
        "ref_energies": np.array([hmodel[i, i] for i in range(nstate)]),
        "e2": e2,
        "min_denoms": min_denoms,
        "n_external": nuniq,
        "rotated": rotated,
        "ref_drift": 0.0,
        "engine": "direct",
    }


def _support_one_electron_model(F, sup_a, sup_b, norb, C):
    """Model-space matrix of a one-electron operator over the support:
    diagonal = sum_occ F_pp; off-diagonal = single-excitation elements with
    the same phase convention as the streaming classes."""
    nsup = sup_a.size
    Fdet = np.zeros((nsup, nsup))
    occ_cache = [( _occ_virt(a, norb)[0], _occ_virt(b, norb)[0])
                 for a, b in zip(sup_a, sup_b)]
    for m in range(nsup):
        oa, ob = occ_cache[m]
        Fdet[m, m] = F[oa, oa].sum() + F[ob, ob].sum()
        for n in range(m + 1, nsup):
            da = int(sup_a[m] ^ sup_a[n])
            db = int(sup_b[m] ^ sup_b[n])
            ca, cb = bin(da).count("1"), bin(db).count("1")
            if ca + cb != 2:
                continue
            if ca == 2:
                word, other = int(sup_a[m]), int(sup_a[n])
                diff = da
            else:
                word, other = int(sup_b[m]), int(sup_b[n])
                diff = db
            p1 = (diff & -diff).bit_length() - 1
            p2 = diff.bit_length() - 1
            i, a = (p1, p2) if (word >> p1) & 1 else (p2, p1)
            nbet = bin(word & ((1 << max(i, a)) - 1) & ~((1 << (min(i, a) + 1)) - 1)).count("1")
            Fdet[m, n] = Fdet[n, m] = F[a, i] * (1.0 if nbet % 2 == 0 else -1.0)
    return C.T @ Fdet @ C
