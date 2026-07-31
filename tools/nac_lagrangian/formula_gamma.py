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

    # ------------------ Stage 2: gamma^formula by generator sweep -----------
    print('\nStage 2: extracting gamma^formula (all generators, Richardson FD)')
    hh = 1e-4
    gam = np.zeros((nstate, nstate, nbf, nbf))
    for p in range(nbf):
        for q in range(p):
            K = np.zeros((nbf, nbf))
            K[p, q] = 1.0
            K[q, p] = -1.0
            # Richardson: (8(S(h)-S(-h)) - (S(2h)-S(-2h)))/(12h), exact to O(h^4)
            Sp1 = replica_S(expm(hh * K))
            Sm1 = replica_S(expm(-hh * K))
            Sp2 = replica_S(expm(2 * hh * K))
            Sm2 = replica_S(expm(-2 * hh * K))
            dS = (8.0 * (Sp1 - Sm1) - (Sp2 - Sm2)) / (12.0 * hh)
            # split between the (p,q)/(q,p) slots so that the FULL contraction
            # sum_pq gam[I,J,p,q] K[p,q] over an antisymmetric K reproduces
            # dS_IJ/dtheta exactly (no double counting)
            gam[:, :, p, q] = 0.5 * dS
            gam[:, :, q, p] = -0.5 * dS
    print('  done.')

    # gate vs Fortran FD on random block generators
    blocks = {'ds': (slice(0, nocb), slice(nocb, noca)),
              'dv': (slice(0, nocb), slice(noca, nbf)),
              'sv': (slice(nocb, noca), slice(noca, nbf)),
              'ss': (slice(nocb, noca), slice(nocb, noca)),
              'full': (slice(0, nbf), slice(0, nbf))}
    thd = 1e-5
    print('\ngate: gamma^formula . K vs Fortran FD')
    ok = True
    for bname, (lo, hi) in blocks.items():
        if bname == 'full':
            K = Kf
        else:
            K = np.zeros((nbf, nbf))
            blk = rng.standard_normal((hi.stop - hi.start, lo.stop - lo.start))
            K[hi, lo] = blk
            if (lo.start, lo.stop) != (hi.start, hi.stop):
                K[lo, hi] = -blk.T
            else:
                K[hi, lo] = blk - blk.T
        dF = (fortran_S(thd, K) - fortran_S(-thd, K)) / (2 * thd)
        for I in range(nstate):
            for J in range(nstate):
                if I >= J:
                    continue
                an = float(np.sum(gam[I, J] * K))
                diff = abs(dF[I, J] - an)
                flag = '' if diff < 5e-6 else '  <-- MISMATCH'
                if diff >= 5e-6:
                    ok = False
                print(f'  {bname:>4} ({I+1},{J+1}): FD={dF[I, J]:+.8f}  '
                      f'gamma.K={an:+.8f}  diff={diff:.1e}{flag}')

    np.savez(out_npz, gamma_formula=gam,
             energies=np.array(mol.energies), probes='generator-sweep',
             note='gamma[I,J,p,q] = dS_IJ/dtheta_pq of the exact tlf=0 formula')
    print(f'\nsaved {out_npz}   gate {"PASSED" if ok else "FAILED"}')


if __name__ == '__main__':
    main()
