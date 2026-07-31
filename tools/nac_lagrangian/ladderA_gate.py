"""LADDER A GATE: the amplitude term by matvec perturbation theory.

Identity under test (raw gauge, folded space, fold V constant in x so
entrywise matvec differencing carries no output-fold artifact):

  amp_scaffold[k,I,J] = dS/dX . dX_J/dx      (proven scaffold, v4)
  amp_PT[k,I,J] = [X_I^T (A(x+h)-A(x-h))/(2h) X_J] / (Om_J - Om_I)

If they match, the analytic amplitude term reduces to the ANALYTIC pieces
of dA (esum + bilinear-2e, already FD-validated) plus L:U^x -- all gated.

Run:  python ladderA_gate.py H2O_energy.inp
"""
import sys
import numpy as np


def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint

    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_lA.log'))
    r.run()
    mol = r.mol

    nstate = mol.config['tdhf']['nstate']
    noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
    nocb = noca - 2
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    nbf = W0.shape[0]
    nvirb = nbf - nocb
    nij = noca * nvirb
    RS = 1.0 / np.sqrt(2.0)
    natom = mol.data['natom']
    ncoord = 3 * natom

    X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
    Xshape = X0_raw.shape
    X0 = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()   # folded (nij, nst)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    E = list(mol.energies)
    Om = [E[k + 1] - E[0] for k in range(nstate)]

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

    # the matvec must be LINEAR: kill int2 screening (trial-density-dependent
    # screening makes apply_A nonlinear -- branch Phase-11 artifact)
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception as e:
        print('WARNING: could not tighten int2e_cutoff:', e)

    def set_trial(col):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = col
        mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)

    def Acols():
        """Apply the CURRENT-state matvec to all nstate reference vectors."""
        out = np.zeros((nij, nstate))
        for s in range(nstate):
            set_trial(X0[:, s])
            oqp.mrsf_matvec_apply(mol)
            out[:, s] = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
        return out

    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False
    h = 1e-3

    genmask = np.ones((noca, nvirb))
    genmask[nocb:noca, 0:2] = 0.0

    # reference minors for the metric-directional amplitude term
    noc = noca - 1

    def s_ia0():
        # at M=I the ityp-3 minors: build literally (same helper as v4, M=I)
        M = np.eye(nbf)
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

    sia0 = s_ia0()
    sij0 = np.eye(noca)
    sab0 = np.eye(nvirb)

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

    def amp_directional(J, dXJ, eps=1e-5):
        Xp_ = [x.copy() for x in Xt0]
        Xm_ = [x.copy() for x in Xt0]
        Xp_[J] = Xt0[J] + eps * dXJ
        Xm_[J] = Xt0[J] - eps * dXJ
        Sp_ = contraction(sij0, sab0, sia0, Xt0, Xp_)
        Sm_ = contraction(sij0, sab0, sia0, Xt0, Xm_)
        return (Sp_[:, J] - Sm_[:, J]) / (2 * eps)

    def displaced(coord):
        mol.update_system(coord)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        # capture the displaced amplitudes BEFORE the matvec trials clobber
        # td_bvec_mo (Acols writes trial vectors into it)
        Xd_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
        # matvec columns of A(x) applied to the fixed reference vectors
        A = Acols()
        mol.data['OQP::td_bvec_mo'] = Xd_raw          # restore
        # aligned displaced amplitudes (orbital signs + state phases)
        mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
        mol.data['OQP::VEC_MO_A_old'] = W0
        mol.data['OQP::VEC_MO_B_old'] = Wb0
        mol.data['OQP::E_MO_A_old'] = e0a
        mol.data['OQP::E_MO_B_old'] = e0b
        oqp.get_structures_ao_overlap(mol)
        M_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        M_f = M_np.reshape(-1).reshape((nbf, nbf)).T
        sg = np.sign(np.diag(M_f))
        sg[sg == 0] = 1.0
        Xd = unfold_all(Xd_raw)
        sgo = sg[:noca]
        sgv = sg[nocb:]
        Xd = [sgo[:, None] * x * sgv[None, :] for x in Xd]
        for s in range(nstate):
            if np.sum(Xt0[s] * Xd[s]) < 0:
                Xd[s] = -Xd[s]
        return A, Xd

    amp_sc = np.zeros((ncoord, nstate, nstate))
    amp_pt = np.zeros((ncoord, nstate, nstate))
    for k in range(ncoord):
        cp = xyz0.copy()
        cp[k] += h
        Ap, Xp = displaced(cp)
        cm = xyz0.copy()
        cm[k] -= h
        Am, Xm = displaced(cm)
        dA_cols = (Ap - Am) / (2 * h)           # (nij, nstate): dA . X0_J
        for J in range(nstate):
            dXJ = (Xp[J] - Xm[J]) / (2 * h)
            col = amp_directional(J, dXJ)
            for I in range(nstate):
                if I == J:
                    continue
                amp_sc[k, I, J] = col[I]
                amp_pt[k, I, J] = float(np.dot(X0[:, I], dA_cols[:, J])) / (Om[J] - Om[I])
        print(f'  coord {k+1}/{ncoord} done', flush=True)

    print('\n=========== LADDER A: scaffold vs matvec-PT (per pair) ===========')
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            u = amp_sc[:, I, J]
            v = amp_pt[:, I, J]
            cos = float(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-300)
            print(f'pair ({I+1},{J+1}): |scaffold|={np.linalg.norm(u):.6f}  '
                  f'|PT|={np.linalg.norm(v):.6f}  cos={cos:+.8f}  '
                  f'max|diff|={np.abs(u-v).max():.2e}')
    print('\nantisymmetrized combination (what enters d):')
    a_sc = 0.5 * (amp_sc - amp_sc.transpose(0, 2, 1))
    a_pt = 0.5 * (amp_pt - amp_pt.transpose(0, 2, 1))
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            u = a_sc[:, I, J]
            v = a_pt[:, I, J]
            cos = float(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-300)
            print(f'pair ({I+1},{J+1}): |scaffold|={np.linalg.norm(u):.6f}  '
                  f'|PT|={np.linalg.norm(v):.6f}  cos={cos:+.8f}  '
                  f'max|diff|={np.abs(u-v).max():.2e}')


if __name__ == '__main__':
    main()
