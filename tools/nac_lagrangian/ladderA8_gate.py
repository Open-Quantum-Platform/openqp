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
    E2 = list(mol.energies)
    Om = [E2[k + 1] - E2[0] for k in range(nstate)]

    # ---------------- Fortran ANALYTIC skeleton pieces ----------------------
    print('computing Fortran skeleton pieces (nac_amp + esum + wsx)...', flush=True)
    oqp.mrsf_nac_amp(mol)
    raw_amp = np.array(mol.data['OQP::nac_amp'], copy=True).reshape(-1)
    amp2e = raw_amp.reshape((nstate, nstate, natom, 3))
    esum_v = {}
    wsx_v = {}
    for i in range(1, nstate + 1):
        for j in range(1, nstate + 1):
            if i == j:
                continue
            oqp.mrsf_nac_esum(mol, i, j)
            esum_v[(i, j)] = np.array(mol.data['OQP::nac_esum'], copy=True).reshape(-1)
            wsx_v[(i, j)] = np.array(mol.data['OQP::nac_wsx'], copy=True).reshape(-1)
    print('  skeleton pieces done.', flush=True)

    # ---------------- interstate L by RHS polarization ----------------------
    import os as _os
    _os.environ['NAC_DUMP_RHS'] = '1'
    nocca = noca
    noccb = nocb

    def rhs_of(target, amp_vec):
        """RHS of the production z-vector for given target-column amplitudes."""
        rr = X0_raw.copy().reshape(-1)
        rr[(target - 1) * nij:target * nij] = amp_vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        mol.data.set_tdhf_target(target)
        oqp.tdhf_mrsf_z_vector(mol)
        out = np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).reshape(-1)
        mol.data['OQP::td_bvec_mo'] = X0_raw
        return out

    print('building interstate L by RHS polarization...', flush=True)
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()   # folded columns
    L_IJ = {}
    rhs_single = {}
    for s in range(nstate):
        rhs_single[s] = rhs_of(s + 1, Xf[:, s])
    for I in range(nstate):
        for J in range(I + 1, nstate):
            rpm = rhs_of(I + 1, Xf[:, I] + Xf[:, J])
            L_IJ[(I, J)] = 0.5 * (rpm - rhs_single[I] - rhs_of(I + 1, Xf[:, J]))
    print('  L done. rhs dim =', len(rhs_single[0]), flush=True)

    def pack_rot(U):
        """Pack a spatial U_pq (nbf x nbf) into the z-vector rotation order:
        doc-socc (socc x doc), doc-virt (virt x doc), socc-virt (virt x socc),
        entry = U[hi, lo]."""
        v = []
        for ii in range(noccb, nocca):          # socc
            for jj in range(noccb):             # doc
                v.append(U[ii, jj])
        for kk in range(nocca, nbf):            # virt
            for jj in range(noccb):
                v.append(U[kk, jj])
        for kk in range(nocca, nbf):
            for ii in range(noccb, nocca):
                v.append(U[kk, ii])
        return np.array(v)


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
        # Sk: ket-half AO-derivative part, C fixed at reference
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        oqp.get_structures_ao_overlap(mol)
        Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
        return M_f, Xd, Skf

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
    Ux_pack = {}
    Ux_full = {}
    for k in range(ncoord):
        cp = xyz0.copy()
        cp[k] += h
        Mp, Xp, Skp = displaced(cp)
        cm = xyz0.copy()
        cm[k] -= h
        Mm, Xm, Skm = displaced(cm)
        T = (Mp - Mm) / (2 * h)
        Sk = (Skp - Skm) / (2 * h)
        Ux = T - Sk
        Ux_pack[k] = pack_rot(Ux)
        Ux_full[k] = Ux
        for J in range(nstate):
            dXJ = (Xp[J] - Xm[J]) / (2 * h)
            damp_col = amp_directional(J, dXJ)     # d_amp[:, J] for all bra I
            for I in range(nstate):
                if I == J:
                    continue
                d_amp[k, I, J] = damp_col[I]
        print(f'  coord {k+1}/{ncoord} done', flush=True)

    print('\n============ LADDER A2: scaffold amp vs analytic skeleton ============')
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            gap = Om[J] - Om[I]
            sc = d_amp[:, I, J]
            g1 = amp2e[J, I].reshape(-1)
            g2 = amp2e[I, J].reshape(-1)
            es = esum_v[(I + 1, J + 1)]
            ws = wsx_v[(I + 1, J + 1)]

            def rep(name, v):
                c = float(np.dot(sc, v)) / (np.linalg.norm(sc) * np.linalg.norm(v) + 1e-300)
                print(f'    {name:30s} |v|={np.linalg.norm(v):.6f} cos={c:+.6f} '
                      f'resid={np.linalg.norm(sc - v):.6f}')

            print(f'pair ({I+1},{J+1}): |scaffold amp| = {np.linalg.norm(sc):.6f}  gap={gap:.6f}')
            rep('amp2e[J,I]/gap', g1 / gap)
            rep('amp2e[I,J]/gap', g2 / gap)
            rep('esum/gap', es / gap)
            rep('(amp2e[J,I]+esum)/gap', (g1 + es) / gap)
            rep('(amp2e[J,I]+esum+wsx)/gap', (g1 + es + ws) / gap)
            rep('(amp2e[J,I]+esum-wsx)/gap', (g1 + es - ws) / gap)
            resid = sc - (g1 + es) / gap
            cw = float(np.dot(resid, ws)) / (np.linalg.norm(resid) * np.linalg.norm(ws) + 1e-300)
            print(f'    residual-vs-wsx: |resid|={np.linalg.norm(resid):.6f} '
                  f'|wsx/gap|={np.linalg.norm(ws / gap):.6f} cos={cw:+.6f}')


    # ---------------- z-vector interchange via the production seam ----------
    # Push L^{IJ} (polarized production RHS, packed 86-dim) as an nbf x nbf
    # matrix whose (hi,lo)-(lo,hi) packing reproduces it, then let the
    # existing hook solve sfrolhs and contract via the gZ-gS gradient seam:
    # the response contribution z.B^x per coordinate, Fock response in LHS.
    def unpack_rot_to_mat(v):
        Lm = np.zeros((nbf, nbf))
        idx = 0
        for ii in range(noccb, nocca):
            for jj in range(noccb):
                Lm[ii, jj] = v[idx]; idx += 1
        for kk in range(nocca, nbf):
            for jj in range(noccb):
                Lm[kk, jj] = v[idx]; idx += 1
        for kk in range(nocca, nbf):
            for ii in range(noccb, nocca):
                Lm[kk, ii] = v[idx]; idx += 1
        return Lm

    def grad_now():
        oqp.tdhf_mrsf_gradient(mol)
        return np.array(mol.get_grad(), copy=True).reshape(-1)

    print('\nsolving interstate z-vectors via the seam...', flush=True)
    resp_z = {}
    for I in range(nstate):
        for J in range(I + 1, nstate):
            Lm = unpack_rot_to_mat(L_IJ[(I, J)])
            mol.data['OQP::nac_orbgrad_L'] = Lm.ravel(order='F').copy()
            mol.data['OQP::td_bvec_mo'] = X0_raw
            mol.data.set_tdhf_target(I + 1)
            oqp.set_mrsf_nac_cphf(mol, I + 1, J + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            conv = bool(mol.mol_energy.Z_Vector_converged)
            if not conv:
                oqp.set_mrsf_nac_cphf(mol, 0, 0)
                print(f'  ({I+1},{J+1}): z-vector NOT converged', flush=True)
                resp_z[(I, J)] = None
                continue
            gZ = grad_now()
            mol.data['OQP::td_p'] = np.zeros_like(
                np.array(mol.data['OQP::td_p'], copy=True))
            mol.data['OQP::WAO'] = np.zeros_like(
                np.array(mol.data['OQP::WAO'], copy=True))
            gS = grad_now()
            oqp.set_mrsf_nac_cphf(mol, 0, 0)
            resp_z[(I, J)] = gZ - gS
            print(f'  ({I+1},{J+1}): |z.B^x| = {np.linalg.norm(gZ-gS):.6f}', flush=True)
    try:
        del mol.data['OQP::nac_orbgrad_L']
    except Exception:
        pass
    print('  skeleton pieces done.', flush=True)
    # ---------------- interstate L by matvec rotation-FD (ALL blocks) -------
    # Bilinear-form FD with FIXED reference vectors: no Rayleigh contamination.
    # Frozen-Fock convention (mrsf_matvec_apply uses the stored AO Fock).
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass
    from scipy.linalg import expm as _expm

    def set_trial(col):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = col
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)

    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()

    def bilinear_S(Wrot):
        """X_I^T A(C_rot) X_J for all pairs; C_rot staged via numpy W."""
        mol.data['OQP::VEC_MO_A'] = Wrot
        mol.data['OQP::VEC_MO_B'] = Wrot
        out = np.zeros((nstate, nstate))
        for s in range(nstate):
            set_trial(Xf[:, s])
            oqp.mrsf_matvec_apply(mol)
            Ax = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
            for I in range(nstate):
                out[I, s] = float(np.dot(Xf[:, I], Ax))
        return out

    print('sweeping L generators (antisym + sym, all blocks)...', flush=True)
    tt = 1e-4
    La = np.zeros((nstate, nstate, nbf, nbf))   # antisym generators
    Ls = np.zeros((nstate, nstate, nbf, nbf))   # symmetric generators
    for p in range(nbf):
        for q in range(p + 1):
            K = np.zeros((nbf, nbf))
            if p == q:
                K[p, p] = 1.0
                Sp_ = bilinear_S(_expm(tt * K) @ W0)
                Sm_ = bilinear_S(_expm(-tt * K) @ W0)
                Ls[:, :, p, p] = (Sp_ - Sm_) / (2 * tt)
                continue
            # antisym
            K[:] = 0.0
            K[p, q] = 1.0
            K[q, p] = -1.0
            Sp_ = bilinear_S(_expm(-tt * K) @ W0)
            Sm_ = bilinear_S(_expm(tt * K) @ W0)
            dS = (Sp_ - Sm_) / (2 * tt)
            La[:, :, p, q] = 0.5 * dS
            La[:, :, q, p] = -0.5 * dS
            # sym: (C e^{tS})^T = e^{+tS} W  (S symmetric!)
            K[:] = 0.0
            K[p, q] = 1.0
            K[q, p] = 1.0
            Sp_ = bilinear_S(_expm(tt * K) @ W0)
            Sm_ = bilinear_S(_expm(-tt * K) @ W0)
            dS = (Sp_ - Sm_) / (2 * tt)
            Ls[:, :, p, q] = 0.5 * dS
            Ls[:, :, q, p] = 0.5 * dS
        print(f'  row {p+1}/{nbf}', flush=True)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    print('  L sweep done.', flush=True)
    np.savez('/bighome/cheolho.choi/nac_audit/ladderA8_data.npz',
             La=La, Ls=Ls,
             Ux=np.array([Ux_full[k] for k in range(ncoord)]),
             d_amp=d_amp, amp2e=amp2e,
             esum=np.array([[esum_v.get((i + 1, j + 1), np.zeros(ncoord))
                             for j in range(nstate)] for i in range(nstate)]),
             zB=np.array([[resp_z.get((min(i, j), max(i, j)), np.zeros(ncoord))
                           if i != j else np.zeros(ncoord)
                           for j in range(nstate)] for i in range(nstate)]),
             Om=np.array(Om))
    print('\n========= LADDER A8: PHYSICAL (no-scan) CLOSURE =========')
    inter = np.zeros((nbf, nbf))
    inter[0:nocb, nocb:] = 1.0
    inter[nocb:, 0:nocb] = 1.0
    inter[nocb:noca, noca:] = 1.0
    inter[noca:, nocb:noca] = 1.0
    same = 1.0 - inter
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            gap = Om[J] - Om[I]
            sc = d_amp[:, I, J]
            g1 = amp2e[J, I].reshape(-1)
            es = esum_v[(I + 1, J + 1)]
            resid = sc - (g1 + es) / gap
            zb = resp_z.get((I, J))
            if zb is None:
                continue
            t1 = zb / gap
            t2a = np.zeros(ncoord)
            t2s = np.zeros(ncoord)
            t3 = np.zeros(ncoord)
            for k in range(ncoord):
                U = Ux_full[k]
                Ua = 0.5 * (U - U.T)
                Us = 0.5 * (U + U.T)
                t2a[k] = float(np.sum((La[I, J] * Ua) * same)) / gap
                t2s[k] = float(np.sum((Ls[I, J] * Us) * same)) / gap
                t3[k] = float(np.sum((Ls[I, J] * Us) * inter)) / gap
            print(f'pair ({I+1},{J+1}): |resid|={np.linalg.norm(resid):.6f}')
            print(f'    |t1 zB|={np.linalg.norm(t1):.6f} |t2a sameA|={np.linalg.norm(t2a):.6f} '
                  f'|t2s sameS|={np.linalg.norm(t2s):.6f} |t3 interS|={np.linalg.norm(t3):.6f}')
            v = t1 + t2a + t2s + t3          # single physical combination
            c = float(np.dot(resid, v)) / (np.linalg.norm(resid)
                                           * np.linalg.norm(v) + 1e-300)
            cl = np.linalg.norm(resid - v)
            print(f'    PHYSICAL t1+t2a+t2s+t3: |v|={np.linalg.norm(v):.6f} '
                  f'cos={c:+.6f} closure={cl:.6f} '
                  f'({100*(1-cl/np.linalg.norm(resid)):.1f}% closed)')
            best = None
            for ss in [(a, b, cc, d) for a in (1, -1) for b in (1, -1)
                       for cc in (1, -1) for d in (1, -1)]:
                vv = ss[0] * t1 + ss[1] * t2a + ss[2] * t2s + ss[3] * t3
                cll = np.linalg.norm(resid - vv)
                if best is None or cll < best[0]:
                    best = (cll, ss, vv)
            cll, ss, vv = best
            print(f'    diagnostic best {ss}: closure={cll:.6f} '
                  f'({100*(1-cll/np.linalg.norm(resid)):.1f}% closed)')


if __name__ == '__main__':
    main()
