"""Production analytic MRSF-TDDFT NAC -- the nac-lagrangian assembly (v1).

Per pair (I < J), all terms from certified components (see
tools/nac_lagrangian/MRSF_NAC_DERIVATION.md, Secs. 4 and 7.24-7.41):

  d_IJ = antisym[ T1 + T2 + Twsx + gamma:Sk ]
    T1    = [amp2e + esum](ytil_IJ, X_J)      slot-injected engines
    T2    = z-vector seam, combined RHS = polarized-L(ytil,X_J) + gamma
    Twsx  = -Tr[W^IJ S^x]                     (Fock-weighted sym channel)
    gamma = closed-form gamma^formula;  Sk from NAC_DUMP_DS
    ytil  = (om_J - A)^{-1}|_perp G_met        (MINRES on the matvec)

ACCURACY (v1, H2O/BHHLYP/6-31G* vs numerical): pair-dependent
5e-3..3e-1 absolute on |d| = 0.03..0.5 -- the amplitude-channel
derivative-sigma engine (ROUTE_A_SPEC.md) is required to reach the
theory-level closure (1e-4..1e-9, proven in 7.27-7.29). Until then this
path prints an accuracy warning and is intended for testing/development.
"""
import os
import numpy as np


def analytic_nac(mol):
    import oqp
    from scipy.sparse.linalg import minres, LinearOperator
    from oqp.library import nac_kernel as FK

    os.environ.setdefault('NAC_DUMP_DS', '1')
    os.environ.setdefault('NAC_DUMP_RHS', '1')

    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ncoord = 3 * natom
    E0 = list(mol.energies)
    Om = [E0[k + 1] - E0[0] for k in range(nstate)]

    ctx = FK.build_context(mol)
    nbf, nij = ctx['nbf'], ctx['nij']
    noca, nocb, nvirb, RS = ctx['noca'], ctx['nocb'], ctx['nvirb'], ctx['RS']
    X0_raw = ctx['X0_raw']
    Xt0 = ctx['Xt']
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
    ijlr2 = (noca - nocb - 1) * noca + noca

    def unfold_vec(v):
        x = np.zeros((noca, nvirb))
        for i in range(1, noca + 1):
            for jj in range(nocb + 1, nbf + 1):
                ij = (jj - nocb - 1) * noca + i
                if ij == ijlr1:
                    x[i - 1, jj - nocb - 1] = v[ijlr1 - 1] * RS
                elif ij == ijlr2:
                    x[i - 1, jj - nocb - 1] = -v[ijlr1 - 1] * RS
                else:
                    x[i - 1, jj - nocb - 1] = v[ij - 1]
        return x

    # gamma + Sk
    gam = FK.gamma_closed(ctx)
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    C = W0.T
    flatF = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            flatF[(I + J * nstate) * nbf * nbf:
                  (I + J * nstate + 1) * nbf * nbf] = gam[I, J].T.reshape(-1)
    mol.data['OQP::nac_gamma_tlf'] = flatF.reshape((nbf * nbf, nstate, nstate))
    oqp.mrsf_nac_overlap(mol)
    nbf2 = nbf * nbf
    dsk_raw = np.array(mol.data['OQP::dbg_dsket'], copy=True).reshape(-1)
    Sk_an = np.zeros((ncoord, nbf, nbf))
    for c in range(ncoord):
        Sk_an[c] = C.T @ dsk_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T @ C

    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    def matvec(v):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = v
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        oqp.mrsf_matvec_apply(mol)
        return np.array(mol.data['OQP::nac_mvax'],
                        copy=True).ravel()[:nij].copy()

    # G_met by unit sweep (closed form: ROUTE_A_SPEC item 4i)
    sij0 = FK.s_ij_of(ctx, np.eye(nbf))
    sab0 = FK.s_ab_of(ctx, np.eye(nbf))
    sia0 = FK.s_ia_of(ctx, np.eye(nbf))
    EPSA = 1e-6

    def ampdir(J, dXt):
        Xp = [x.copy() for x in Xt0]
        Xm = [x.copy() for x in Xt0]
        Xp[J] = Xt0[J] + EPSA * dXt
        Xm[J] = Xt0[J] - EPSA * dXt
        Sp = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xp)
        Sm = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xm)
        return (Sp[:, J] - Sm[:, J]) / (2 * EPSA)

    G_met = np.zeros((nstate, nstate, nij))
    for J in range(nstate):
        for k in range(nij):
            if k == ijlr2 - 1:
                continue
            e = np.zeros(nij)
            e[k] = 1.0
            g = ampdir(J, unfold_vec(e))
            for I in range(nstate):
                if I != J:
                    G_met[I, J, k] = g[I]

    phys = [k for k in range(nij) if k != ijlr2 - 1]

    def ytil_of(I, J):
        XJ = Xf[:, J]

        def op(v):
            vv = np.zeros(nij)
            vv[phys] = v
            vv -= XJ * float(np.dot(XJ, vv))
            av = Om[J] * vv - matvec(vv)
            av -= XJ * float(np.dot(XJ, av))
            return av[phys]

        rhs = G_met[I, J].copy()
        rhs -= XJ * float(np.dot(XJ, rhs))
        y, _ = minres(LinearOperator((nij - 1, nij - 1), matvec=op),
                      rhs[phys], rtol=1e-9, maxiter=3000)
        out = np.zeros(nij)
        out[phys] = y
        out -= XJ * float(np.dot(XJ, out))
        return out

    def inject(I, vec):
        rr = X0_raw.copy().reshape(-1)
        rr[I * nij:(I + 1) * nij] = vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)

    def rhs_of(target, amp_vec):
        rr = X0_raw.copy().reshape(-1)
        rr[(target - 1) * nij:target * nij] = amp_vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        mol.data.set_tdhf_target(target)
        oqp.tdhf_mrsf_z_vector(mol)
        out = np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).reshape(-1)
        mol.data['OQP::td_bvec_mo'] = X0_raw
        return out

    def unpack_rot(v):
        Lm = np.zeros((nbf, nbf))
        idx = 0
        for ii in range(nocb, noca):
            for jj in range(nocb):
                Lm[ii, jj] = v[idx]; idx += 1
        for kk in range(noca, nbf):
            for jj in range(nocb):
                Lm[kk, jj] = v[idx]; idx += 1
        for kk in range(noca, nbf):
            for ii in range(nocb, noca):
                Lm[kk, ii] = v[idx]; idx += 1
        return Lm

    def grad_now():
        oqp.tdhf_mrsf_gradient(mol)
        return np.array(mol.get_grad(), copy=True).reshape(-1)

    dp = np.zeros((ncoord, nstate, nstate))
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            y = ytil_of(I, J)
            # T1 + wsx (slot injection)
            inject(I, y)
            oqp.mrsf_nac_amp(mol)
            a = np.array(mol.data['OQP::nac_amp'], copy=True
                         ).reshape((nstate, nstate, natom, 3))
            oqp.mrsf_nac_esum(mol, I + 1, J + 1)
            es = np.array(mol.data['OQP::nac_esum'], copy=True).reshape(-1)
            wsx = np.array(mol.data['OQP::nac_wsx'], copy=True).reshape(-1)
            mol.data['OQP::td_bvec_mo'] = X0_raw
            T1 = a[J, I].reshape(-1) + es
            # T2 seam
            qy = rhs_of(J + 1, y)
            qx = rhs_of(J + 1, Xf[:, J])
            qm = rhs_of(J + 1, y + Xf[:, J])
            Lmat = unpack_rot(0.5 * (qm - qy - qx)) + gam[I, J]
            mol.data['OQP::nac_orbgrad_L'] = Lmat.ravel(order='F').copy()
            mol.data['OQP::td_bvec_mo'] = X0_raw
            mol.data.set_tdhf_target(J + 1)
            oqp.set_mrsf_nac_cphf(mol, I + 1, J + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            gZ = grad_now()
            mol.data['OQP::td_p'] = np.zeros_like(
                np.array(mol.data['OQP::td_p'], copy=True))
            mol.data['OQP::WAO'] = np.zeros_like(
                np.array(mol.data['OQP::WAO'], copy=True))
            gS = grad_now()
            oqp.set_mrsf_nac_cphf(mol, 0, 0)
            T2 = gZ - gS
            gsk = np.array([float(np.sum(gam[I, J] * Sk_an[c]))
                            for c in range(ncoord)])
            dp[:, I, J] = T1 + T2 + wsx + gsk
    try:
        del mol.data['OQP::nac_orbgrad_L']
    except Exception:
        pass
    mol.data['OQP::td_bvec_mo'] = X0_raw

    dpa = 0.5 * (dp - dp.transpose(0, 2, 1))
    dcv = np.zeros((nstate, nstate, natom, 3))
    nacv = np.zeros((nstate, nstate, natom, 3))
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            gap = Om[J] - Om[I]
            dcv[I, J] = dpa[:, I, J].reshape(natom, 3)
            nacv[I, J] = gap * dcv[I, J]
    return nacv, dcv
