"""ASSEMBLY GATE: prove the master decomposition of the MRSF formula-NAC

    d_num[i,j] = Xt_I . dXt_J/dx  +  sum_pq gamma^formula[I,J,p,q] T_pq(x)
    T = dM/dx,  M = C0^T S_AO(x0,x) C(x)   (= U^x + Sk in one object)

end to end against the production numerical NAC, all in ONE process (shared
R0 state phases -> SIGNED comparison, no gauge freedom).

Steps: R0 run -> in-process gamma^formula (generator sweep of the exact
replica) -> +-h displacements along all 3N coords (aligned M and amplitudes)
-> d_pred -> production numerical_nac -> signed per-pair comparison.

Run:  python assembly_gate.py H2O_energy.inp
"""
import sys
import numpy as np


def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint, NAC
    from scipy.linalg import expm

    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_ag.log'))
    r.run()
    mol = r.mol

    nstate = mol.config['tdhf']['nstate']
    noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
    nocb = noca - 2
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    nbf = W0.shape[0]
    nvirb = nbf - nocb
    nij = noca * nvirb
    noc = noca - 1
    RS = 1.0 / np.sqrt(2.0)
    natom = mol.data['natom']
    ncoord = 3 * natom

    X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)

    # ---- d_num FIRST, on the pristine state (harness-order hygiene) --------
    nac = NAC(mol)
    nacv, dcv, flags = nac.numerical_nac()
    print('numerical flags:', set(flags), flush=True)
    # restore the R0 anchor data untouched
    mol.data['OQP::td_bvec_mo'] = X0_raw

    def unfold_m(bv, st):
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

    def unfold_all(raw):
        Xm = raw.reshape(-1).reshape((nstate, nij)).T.copy()
        return [unfold_m(Xm, s + 1) for s in range(nstate)]

    Xt0 = unfold_all(X0_raw)

    # ---------------- exact formula machinery (validated 2e-16) -------------
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
        G = np.zeros((noca, nvirb))
        for i1 in range(1, noca + 1):
            for j1 in range(nvirb):
                ia1 = nocb + j1 + 1
                D = np.zeros((noc, noc))
                for i in range(1, i1):
                    for ipp in range(1, i1):
                        D[i - 1, ipp - 1] = M[i - 1, ipp - 1]
                    for ipp in range(i1, noc - 1):
                        D[i - 1, ipp - 1] = M[i - 1, ipp]
                    D[i - 1, noc - 2] = M[i - 1, i1 - 1]
                    D[i - 1, noc - 1] = M[i - 1, ia1 - 1]
                for i in range(i1, noc - 1):
                    for ipp in range(1, i1):
                        D[i - 1, ipp - 1] = M[i, ipp - 1]
                    for ipp in range(i1, noc - 1):
                        D[i - 1, ipp - 1] = M[i, ipp]
                    D[i - 1, noc - 2] = M[i, i1 - 1]
                    D[i - 1, noc - 1] = M[i, ia1 - 1]
                for i in (noc - 1, noc):
                    for ipp in range(1, i1):
                        D[i - 1, ipp - 1] = M[i, ipp - 1]
                    for ipp in range(i1, noc - 1):
                        D[i - 1, ipp - 1] = M[i, ipp]
                    D[i - 1, noc - 2] = M[i, i1 - 1]
                    D[i - 1, noc - 1] = M[i, ia1 - 1]
                G[i1 - 1, j1] = np.linalg.det(D)
        return G

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
                acc = float(np.sum((cog @ s_ab) * (s_ij @ cng)))
                acc += float(np.sum((cog @ s_ia.T) * (cng @ s_ia.T).T))
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
        return contraction(s_ij_of(M), s_ab_of(M), s_ia_of(M), Xt0, Xt0)

    # ---------------- gamma^formula in-process ------------------------------
    print('extracting gamma^formula (generator sweep)...', flush=True)
    hh = 1e-4
    gam = np.zeros((nstate, nstate, nbf, nbf))
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
            gam[:, :, p, q] = 0.5 * dS
            gam[:, :, q, p] = -0.5 * dS
    print('  gamma done.', flush=True)

    # ---------------- displacement loop -------------------------------------
    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False
    h = 1e-3

    def displaced(coord):
        """Return (M_f aligned, Xt aligned) at the displaced geometry."""
        mol.update_system(coord)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
        mol.data['OQP::VEC_MO_A_old'] = W0
        mol.data['OQP::VEC_MO_B_old'] = Wb0
        mol.data['OQP::E_MO_A_old'] = e0a
        mol.data['OQP::E_MO_B_old'] = e0b
        oqp.get_structures_ao_overlap(mol)
        M_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        M_f = M_np.reshape(-1).reshape((nbf, nbf)).T
        # orbital sign alignment to R0 (diag > 0), tracked so the amplitudes
        # can be transformed into the SAME orbital-sign gauge
        sg = np.sign(np.diag(M_f))
        sg[sg == 0] = 1.0
        M_f = M_f * sg[None, :]
        Xd = unfold_all(np.array(mol.data['OQP::td_bvec_mo'], copy=True))
        # amplitude transform under orbital sign flips: X(i,a) -> sg_i sg_a X
        sgo = sg[:noca]
        sgv = sg[nocb:]
        Xd = [sgo[:, None] * x * sgv[None, :] for x in Xd]
        # per-state phase alignment to R0 amplitudes
        for s in range(nstate):
            if np.sum(Xt0[s] * Xd[s]) < 0:
                Xd[s] = -Xd[s]
        return M_f, Xd

    # reference minors once (amplitude directional derivatives live on them)
    I_mo = np.eye(nbf)
    sij0, sab0, sia0 = s_ij_of(I_mo), s_ab_of(I_mo), s_ia_of(I_mo)

    def amp_directional(J, dXJ, eps=1e-5):
        """Exact formula amplitude term: directional derivative of the
        normalized contraction w.r.t. the ket amplitudes of state J along
        dXJ, at the reference minors. NOT the plain dot product X_I . dX_J:
        the sqrt2 terms and the socc baseline of s_ia make the formula's
        amplitude metric nontrivial."""
        Xp_ = [x.copy() for x in Xt0]
        Xm_ = [x.copy() for x in Xt0]
        Xp_[J] = Xt0[J] + eps * dXJ
        Xm_[J] = Xt0[J] - eps * dXJ
        Sp_ = contraction(sij0, sab0, sia0, Xt0, Xp_)
        Sm_ = contraction(sij0, sab0, sia0, Xt0, Xm_)
        return (Sp_[:, J] - Sm_[:, J]) / (2 * eps)

    d_amp = np.zeros((ncoord, nstate, nstate))
    d_orb = np.zeros((ncoord, nstate, nstate))
    d_repl = np.zeros((ncoord, nstate, nstate))
    for k in range(ncoord):
        cp = xyz0.copy()
        cp[k] += h
        Mp, Xp = displaced(cp)
        cm = xyz0.copy()
        cm[k] -= h
        Mm, Xm = displaced(cm)
        T = (Mp - Mm) / (2 * h)
        # TOTAL formula derivative from the exact replica (bypasses production)
        Srp = contraction(s_ij_of(Mp), s_ab_of(Mp), s_ia_of(Mp), Xt0, Xp)
        Srm = contraction(s_ij_of(Mm), s_ab_of(Mm), s_ia_of(Mm), Xt0, Xm)
        d_repl[k] = (Srp - Srm) / (2 * h)
        for J in range(nstate):
            dXJ = (Xp[J] - Xm[J]) / (2 * h)
            damp_col = amp_directional(J, dXJ)     # d_amp[:, J] for all bra I
            for I in range(nstate):
                if I == J:
                    continue
                d_amp[k, I, J] = damp_col[I]
                d_orb[k, I, J] = float(np.sum(gam[I, J] * T))
        print(f'  coord {k+1}/{ncoord} done', flush=True)

    d_pred = d_amp + d_orb

    # PRODUCTION COMPARABILITY: the pipeline antisymmetrizes the overlap
    # numerator ((S - S^T)/2dt), i.e. d_num[i,j] = (dS_ij - dS_ji)/2. The
    # formula's dS itself carries a SYMMETRIC part (not an exact-eigenstate
    # overlap), so all three legs must be compared on the antisymmetrized
    # combination -- which also makes d antisymmetry EMERGENT, not stamped.
    d_repl = 0.5 * (d_repl - d_repl.transpose(0, 2, 1))
    d_amp = 0.5 * (d_amp - d_amp.transpose(0, 2, 1))
    d_orb = 0.5 * (d_orb - d_orb.transpose(0, 2, 1))
    d_pred = d_amp + d_orb

    print('\n================= ASSEMBLY GATE (three-way, antisymmetrized) =====')
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            dn = dcv[I, J].reshape(-1)
            dr = d_repl[:, I, J]
            dp = d_pred[:, I, J]
            da = d_amp[:, I, J]
            do = d_orb[:, I, J]

            def cs(u, v):
                return float(np.dot(u, v)) / (np.linalg.norm(u)
                                              * np.linalg.norm(v) + 1e-300)
            print(f'pair ({I+1},{J+1}):')
            print(f'  |d_num|={np.linalg.norm(dn):.6f}  |d_repl|={np.linalg.norm(dr):.6f}  '
                  f'|d_chain|={np.linalg.norm(dp):.6f}  '
                  f'(|amp|={np.linalg.norm(da):.6f} |orb|={np.linalg.norm(do):.6f})')
            print(f'  LEG1 replica==chain : cos={cs(dr, dp):+.8f} max|diff|={np.abs(dr-dp).max():.2e}')
            print(f'  LEG2 replica==num   : cos={cs(dr, dn):+.8f} max|diff|={np.abs(dr-dn).max():.2e}')
            print(f'  LEG3 chain==num     : cos={cs(dp, dn):+.8f} max|diff|={np.abs(dp-dn).max():.2e}')
            print(f'  num : {np.round(dn, 5)}')
            print(f'  repl: {np.round(dr, 5)}')
            print(f'  chn : {np.round(dp, 5)}')


if __name__ == '__main__':
    main()
