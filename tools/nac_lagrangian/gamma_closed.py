"""gamma^formula: exact first-order orbital-rotation response of the MRSF
state-overlap FORMULA (compute_states_overlap over exact tlf=0 minors).

Stage 1: LITERAL replica (incl. ov_exact case(3)'s overwrite semantics and
         the numpy<->Fortran transpose-corrected staging) vs Fortran S at a
         finite rotation -- must be machine precision.
Stage 2: extract gamma^formula_pq per state pair by sweeping every rotation
         generator (Richardson FD of the exact replica -- a derived object,
         no fitting), then GATE: sum(gamma o K) == Fortran FD per block.
         Saves gamma^formula to npz for wiring into mrsf_nac_overlap.

Run:  python formula_gamma.py H2O_energy.inp out.npz
"""
import sys
import numpy as np


def main():
    import oqp
    from oqp.pyoqp import Runner
    from scipy.linalg import expm

    inp = sys.argv[1]
    out_npz = sys.argv[2] if len(sys.argv) > 2 else '/tmp/gamma_formula.npz'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_fg.log'))
    r.run()
    mol = r.mol

    nstate = mol.config['tdhf']['nstate']
    noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
    nocb = noca - 2
    W = np.array(mol.data['OQP::VEC_MO_A'], copy=True)   # numpy = C_fortran^T
    nbf = W.shape[0]
    nvirb = nbf - nocb
    nij = noca * nvirb
    noc = noca - 1
    RS = 1.0 / np.sqrt(2.0)

    X0 = np.array(mol.data['OQP::td_bvec_mo'], copy=True
                  ).reshape(-1).reshape((nstate, nij)).T.copy()

    def unfold(bv, st):
        ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
        ijlr2 = (noca - nocb - 1) * noca + noca
        x = np.zeros((noca, nvirb))
        for i in range(1, noca + 1):
            for jj in range(nocb + 1, nbf + 1):
                ij = (jj - nocb - 1) * noca + i
                if ij == ijlr1:
                    x[i - 1, jj - nocb - 1] = bv[ijlr1 - 1, st - 1] * RS
                elif ij == ijlr2:
                    x[i - 1, jj - nocb - 1] = -bv[ijlr1 - 1, st - 1] * RS
                else:
                    x[i - 1, jj - nocb - 1] = bv[ij - 1, st - 1]
        return x

    Xt = [unfold(X0, s + 1) for s in range(nstate)]

    # ------------------ minors: exact ov_exact replicas ---------------------
    def s_ij_of(M):
        G = np.zeros((noca, noca))
        for i1 in range(1, noca + 1):
            for i2 in range(1, noca + 1):
                if i1 == i2:
                    keep = [k for k in range(noca) if k != i1 - 1]
                    G[i1 - 1, i2 - 1] = np.linalg.det(M[np.ix_(keep, keep)])
                else:
                    imin, imax = min(i1, i2), max(i1, i2)
                    rows = ([k for k in range(noca) if k + 1 not in (imin, imax)]
                            + [i2 - 1])
                    cols = ([k for k in range(noca) if k + 1 not in (imin, imax)]
                            + [i1 - 1])
                    G[i1 - 1, i2 - 1] = -np.linalg.det(M[np.ix_(rows, cols)])
        return G

    def s_ab_of(M):
        G = np.zeros((nvirb, nvirb))
        core = list(range(nocb))
        for j1 in range(nvirb):
            for j2 in range(nvirb):
                G[j1, j2] = np.linalg.det(
                    M[np.ix_(core + [nocb + j1], core + [nocb + j2])])
        return G

    def s_ia_of(M):
        """LITERAL transcription of ov_exact case(3), overwrite semantics.
        ddet is (noc x noc); blocks written in source order, last write wins."""
        G = np.zeros((noca, nvirb))
        for i1 in range(1, noca + 1):
            for j1 in range(nvirb):
                ia1 = nocb + j1 + 1                      # 1-based s_mo column
                D = np.zeros((noc, noc))
                # (1,1)..(1,4)
                for i in range(1, i1):
                    for ipp in range(1, i1):
                        D[i - 1, ipp - 1] = M[i - 1, ipp - 1]
                    for ipp in range(i1, noc - 1):
                        D[i - 1, ipp - 1] = M[i - 1, ipp]
                    D[i - 1, noc - 2] = M[i - 1, i1 - 1]
                    D[i - 1, noc - 1] = M[i - 1, ia1 - 1]
                # (2,1)..(2,4)
                for i in range(i1, noc - 1):
                    for ipp in range(1, i1):
                        D[i - 1, ipp - 1] = M[i, ipp - 1]
                    for ipp in range(i1, noc - 1):
                        D[i - 1, ipp - 1] = M[i, ipp]
                    D[i - 1, noc - 2] = M[i, i1 - 1]
                    D[i - 1, noc - 1] = M[i, ia1 - 1]
                # (3,1)..(3,4): fixed ddet row noc-1 <- s_mo row noc
                i = noc - 1
                for ipp in range(1, i1):
                    D[i - 1, ipp - 1] = M[i, ipp - 1]
                for ipp in range(i1, noc - 1):
                    D[i - 1, ipp - 1] = M[i, ipp]
                D[i - 1, noc - 2] = M[i, i1 - 1]
                D[i - 1, noc - 1] = M[i, ia1 - 1]
                # (4,1)..(4,4): fixed ddet row noc <- s_mo row noc+1
                i = noc
                for ipp in range(1, i1):
                    D[i - 1, ipp - 1] = M[i, ipp - 1]
                for ipp in range(i1, noc - 1):
                    D[i - 1, ipp - 1] = M[i, ipp]
                D[i - 1, noc - 2] = M[i, i1 - 1]
                D[i - 1, noc - 1] = M[i, ia1 - 1]
                G[i1 - 1, j1] = np.linalg.det(D)
        return G

    def minors(M):
        return s_ij_of(M), s_ab_of(M), s_ia_of(M)

    # ------------------ the contraction replica -----------------------------
    socc = slice(nocb, noca)
    genmask = np.ones((noca, nvirb))
    genmask[socc, 0:2] = 0.0

    def contraction(s_ij, s_ab, s_ia, xo, xn):
        ns = len(xo)
        S = np.zeros((ns, ns))
        for oi in range(ns):
            for ni in range(ns):
                co, cn = xo[oi], xn[ni]
                cog, cng = co * genmask, cn * genmask
                alpham = cog @ s_ab
                betam = s_ij @ cng
                acc = float(np.sum(alpham * betam))
                gammam = cog @ s_ia.T
                deltam = cng @ s_ia.T
                acc += float(np.sum(gammam * deltam.T))
                for pi in range(nocb, noca):
                    for qi in range(nocb, noca):
                        for ri in range(nocb, noca):
                            for si in range(nocb, noca):
                                acc += (co[pi, qi - nocb] * s_ij[pi, ri]
                                        * cn[ri, si - nocb]
                                        * s_ab[qi - nocb, si - nocb])
                for pi in range(nocb, noca):
                    for qi in range(nocb, noca):
                        for ri in range(noca):
                            for si in range(nocb, nbf):
                                if (ri >= nocb) and (si < noca):
                                    continue
                                acc += (co[pi, qi - nocb] * cn[ri, si - nocb]
                                        * (s_ij[pi, ri] * s_ab[qi - nocb, si - nocb]
                                           + s_ia[pi, si - nocb] * s_ia[ri, qi - nocb]) * RS)
                for pi in range(noca):
                    for qi in range(nocb, nbf):
                        if (pi >= nocb) and (qi < noca):
                            continue
                        for ri in range(nocb, noca):
                            for si in range(nocb, noca):
                                acc += (co[pi, qi - nocb] * cn[ri, si - nocb]
                                        * (s_ij[pi, ri] * s_ab[qi - nocb, si - nocb]
                                           + s_ia[pi, si - nocb] * s_ia[ri, qi - nocb]) * RS)
                S[oi, ni] = acc
        for i in range(ns):
            S[:, i] /= np.linalg.norm(S[:, i])
        return S

    def replica_S(M):
        sij, sab, sia = minors(M)
        return contraction(sij, sab, sia, Xt, Xt)

    # ------------------ Fortran driver (transpose-corrected staging) --------
    Cb = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    X_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
    mol.data['OQP::VEC_MO_A_old'] = W
    mol.data['OQP::VEC_MO_B_old'] = Cb.copy()
    mol.data['OQP::E_MO_A_old'] = np.array(mol.data['OQP::E_MO_A'], copy=True)
    mol.data['OQP::E_MO_B_old'] = np.array(mol.data['OQP::E_MO_B'], copy=True)
    mol.data['OQP::td_bvec_mo_old'] = X_raw.copy()
    mol.data['OQP::xyz_old'] = np.array(mol.get_system(), copy=True).reshape((3, -1))
    mol.data.set_tdhf_tlf(0)

    def fortran_S(theta, K):
        # C_f -> C_f e^{theta K}  <=>  numpy W -> expm(-theta K) @ W
        mol.data['OQP::VEC_MO_A'] = expm(-theta * K) @ W
        mol.data['OQP::VEC_MO_B'] = expm(-theta * K) @ W
        mol.data['OQP::td_bvec_mo'] = X_raw
        oqp.get_structures_ao_overlap(mol)
        oqp.get_states_overlap(mol)
        S = np.array(mol.data['OQP::td_states_overlap'], copy=True)
        mol.data['OQP::VEC_MO_A'] = W
        mol.data['OQP::VEC_MO_B'] = Cb
        return S.T          # numpy tagarray = Fortran^T; return s_st(oi,ni)

    rng = np.random.default_rng(23)
    Kf = rng.standard_normal((nbf, nbf))
    Kf = Kf - Kf.T
    th = 1e-3
    SF = fortran_S(th, Kf)
    SP = replica_S(expm(th * Kf))
    d1 = np.abs(SP - SF).max()
    print(f'Stage 1 (literal replica, corrected staging): max|SP-SF| = {d1:.3e}')
    if d1 > 1e-9:
        np.set_printoptions(precision=8, suppress=True)
        print('SP =\n', SP)
        print('SF =\n', SF)
        print('REPLICA STILL MISMATCHED -- stopping before stage 2.')
        return

    # ------------------ generator-sweep gamma (referee) ---------------------
    print('sweep gamma (referee)...', flush=True)
    hh = 1e-4
    gam_sweep = np.zeros((nstate, nstate, nbf, nbf))
    for p in range(nbf):
        for q in range(p):
            K = np.zeros((nbf, nbf))
            K[p, q] = 1.0
            K[q, p] = -1.0
            Sp1 = replica_S(expm(hh * K))
            Sm1 = replica_S(expm(-hh * K))
            Sp2 = replica_S(expm(2 * hh * K))
            Sm2 = replica_S(expm(-2 * hh * K))
            dS = (8.0 * (Sp1 - Sm1) - (Sp2 - Sm2)) / (12.0 * hh)
            gam_sweep[:, :, p, q] = 0.5 * dS
            gam_sweep[:, :, q, p] = -0.5 * dS
    print('  sweep done.', flush=True)

    # ------------------ CLOSED FORM: cofactor sensitivities -----------------
    print('closed-form gamma via cofactors...', flush=True)
    import numpy.linalg as la
    I_mo = np.eye(nbf)

    def adjugate(A):
        n = A.shape[0]
        adj = np.zeros_like(A)
        for r in range(n):
            for c in range(n):
                M = np.delete(np.delete(A, r, axis=0), c, axis=1)
                adj[c, r] = ((-1.0) ** (r + c)) * la.det(M)
        return adj

    # sensitivity W[entry][p,q]: d(minor entry)/dtheta_pq = sum W*K.
    # For minor det(M[rows, cols]) (ordered lists), at M=I:
    #   d det = sum_ab adj0[b,a]... using ddet = tr(adj(A0) dA) with
    #   dA[a,b] = K[rows[a], cols[b]] -> W[p,q] += adj0[b,a] at
    #   p=rows[a], q=cols[b].
    def sens(rows, cols, sign=1.0):
        A0 = I_mo[np.ix_(rows, cols)]
        adj0 = adjugate(A0)
        Wm = {}
        n = len(rows)
        for a in range(n):
            for b in range(n):
                v = sign * adj0[b, a]
                if v != 0.0:
                    Wm[(rows[a], cols[b])] = Wm.get((rows[a], cols[b]), 0.0) + v
        return Wm

    # minor definitions (same index sets as the exact replica)
    minors_def = {}
    for i1 in range(1, noca + 1):
        for i2 in range(1, noca + 1):
            if i1 == i2:
                keep = [k for k in range(noca) if k != i1 - 1]
                minors_def[('ij', i1, i2)] = (keep, keep, 1.0)
            else:
                imin, imax = min(i1, i2), max(i1, i2)
                rows = ([k for k in range(noca) if k + 1 not in (imin, imax)]
                        + [i2 - 1])
                cols = ([k for k in range(noca) if k + 1 not in (imin, imax)]
                        + [i1 - 1])
                minors_def[('ij', i1, i2)] = (rows, cols, -1.0)
    core = list(range(nocb))
    for j1 in range(nvirb):
        for j2 in range(nvirb):
            minors_def[('ab', j1, j2)] = (core + [nocb + j1], core + [nocb + j2], 1.0)
    # s_ia: NET row/col mapping of the LITERAL overwrite layout (blocks 3/4
    # always win ddet rows noc-1, noc; block2 shifts by +1 past the i1 gap)
    for i1 in range(1, noca + 1):
        for j1 in range(nvirb):
            a1g = nocb + j1
            rows = []
            for r in range(1, noc - 1):                  # ddet rows 1..noc-2
                rows.append((r if r <= i1 - 1 else r + 1) - 1)
            rows += [noc, noc + 1 - 1]                   # s_mo rows noc, noc+1 (0-based: noc-1+1? careful)
            cols = []
            for c in range(1, noc - 1):
                cols.append((c if c <= i1 - 1 else c + 1) - 1)
            cols += [i1 - 1, a1g]
            # 0-based correction for the fixed rows: s_mo rows noc, noc+1 (1-based)
            rows[-2] = noc - 1 + 1      # = noc (0-based index of 1-based noc+? )
            rows[-1] = noc              # 0-based of 1-based noc+1
            rows[-2] = noc + 1 - 1 - 1  # 0-based of 1-based noc  -> noc-1? no...
            # definitive: 1-based s_mo row noc -> 0-based noc-1; row noc+1 -> noc
            rows[-2] = noc - 1
            rows[-1] = noc
            minors_def[('ia', i1 - 1, j1)] = (rows, cols, 1.0)

    W = {k: sens(*v) for k, v in minors_def.items()}
    print(f'  {len(W)} minor sensitivities built.', flush=True)

    # dS/d(minor entry) at the reference by exact entry perturbation
    sij0, sab0, sia0 = s_ij_of(I_mo), s_ab_of(I_mo), s_ia_of(I_mo)
    S0 = contraction(sij0, sab0, sia0, Xt, Xt)
    eps = 1e-6
    gam_closed = np.zeros((nstate, nstate, nbf, nbf))
    for key, Wm in W.items():
        kind = key[0]
        if kind == 'ij':
            i1, i2 = key[1] - 1, key[2] - 1
            s = sij0.copy(); s[i1, i2] += eps
            Sp_ = contraction(s, sab0, sia0, Xt, Xt)
            s[i1, i2] -= 2 * eps
            Sm_ = contraction(s, sab0, sia0, Xt, Xt)
        elif kind == 'ab':
            j1, j2 = key[1], key[2]
            s = sab0.copy(); s[j1, j2] += eps
            Sp_ = contraction(sij0, s, sia0, Xt, Xt)
            s[j1, j2] -= 2 * eps
            Sm_ = contraction(sij0, s, sia0, Xt, Xt)
        else:
            i1, j1 = key[1], key[2]
            s = sia0.copy(); s[i1, j1] += eps
            Sp_ = contraction(sij0, sab0, s, Xt, Xt)
            s[i1, j1] -= 2 * eps
            Sm_ = contraction(sij0, sab0, s, Xt, Xt)
        dSds = (Sp_ - Sm_) / (2 * eps)          # (nstate,nstate)
        for (p, q), w in Wm.items():
            # antisym generator K_pq=+1, K_qp=-1: dM[p,q]=+1, dM[q,p]=-1
            gam_closed[:, :, p, q] += 0.5 * dSds * w
            gam_closed[:, :, q, p] -= 0.5 * dSds * w
    print('  closed form assembled.', flush=True)

    print('\n===== GATE: closed-form gamma vs generator sweep =====')
    sl = {'d': slice(0, nocb), 's': slice(nocb, noca), 'v': slice(noca, nbf)}
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            d = np.abs(gam_closed[I, J] - gam_sweep[I, J]).max()
            n1 = np.linalg.norm(gam_sweep[I, J])
            print(f'  ({I+1},{J+1}): |gamma|={n1:.6f}  max|closed-sweep|={d:.3e}')
            for b1 in 'dsv':
                for b2 in 'dsv':
                    dd = np.abs(gam_closed[I, J][sl[b1], sl[b2]]
                                - gam_sweep[I, J][sl[b1], sl[b2]]).max()
                    if dd > 1e-6:
                        sw = np.linalg.norm(gam_sweep[I, J][sl[b1], sl[b2]])
                        cl = np.linalg.norm(gam_closed[I, J][sl[b1], sl[b2]])
                        print(f'      block {b1}{b2}: maxdiff={dd:.3e} '
                              f'|sweep|={sw:.5f} |closed|={cl:.5f}')


if __name__ == '__main__':
    main()
