"""Production analytic MRSF-TDDFT NAC -- the nac-lagrangian assembly (v2).

Per pair (I < J), all terms from certified components (see
tools/nac_lagrangian/MRSF_NAC_DERIVATION.md, Secs. 4 and 7.24-7.41):

  d_IJ = antisym[ T1 - seam(X) + X:V + gamma:Sk ]
    T1    = [amp2e + esum](ytil_IJ, X_J)      slot-injected engines
    X     = MT_frozen + MT_response + gamma
    seam  = direct orbital-gradient injection into the Z-vector
    X:V   = explicit symmetric-U / overlap-derivative contraction
    gamma = closed-form gamma^formula;  Sk from NAC_DUMP_DS
    ytil  = (om_J - A)^{-1}|_perp G_met        (MINRES on the matvec)

MT_frozen currently comes from the resident Fortran ``mrsf_nac_wpair``
reference harvest.  Its internal orbital-generator loop is numerically
certified but will be replaced by the closed-form bilinear adjoint without
changing this Python orchestration layer.
"""
import os
import numpy as np


def analytic_nac(mol):
    import oqp
    from scipy.sparse.linalg import minres, LinearOperator
    from oqp.library import nac_kernel as FK

    os.environ.setdefault('NAC_DUMP_DS', '1')
    os.environ.setdefault('NAC_DUMP_RHS', '1')
    os.environ.setdefault('NAC_DUMP_PIJ', '1')

    if mol.config['tdhf']['multiplicity'] != 1:
        raise NotImplementedError(
            'analytic MRSF NAC v2 currently implements the singlet fold only'
        )

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

    def pack_sym(matrix):
        packed = np.zeros(nbf * (nbf + 1) // 2)
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                packed[idx] = matrix[p, q]
                idx += 1
        return packed

    def unpack_sym(packed):
        matrix = np.zeros((nbf, nbf))
        idx = 0
        for q in range(nbf):
            for p in range(q + 1):
                matrix[p, q] = packed[idx]
                matrix[q, p] = packed[idx]
                idx += 1
        return matrix

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
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
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
    dsf_raw = np.array(mol.data['OQP::dbg_dsfull'], copy=True).reshape(-1)
    Sk_an = np.zeros((ncoord, nbf, nbf))
    Sx_MO = np.zeros((ncoord, nbf, nbf))
    for c in range(ncoord):
        Sk_an[c] = C.T @ dsk_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T @ C
        dsf = dsf_raw[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T
        Sx_MO[c] = C.T @ dsf @ C

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

    def grad_now():
        oqp.tdhf_mrsf_gradient(mol)
        return np.array(mol.get_grad(), copy=True).reshape(-1)

    def frozen_orbital_gradient(I, J, y):
        """Call the resident Fortran reference harvest for one ordered pair."""
        mol.data['OQP::nac_ytil'] = np.array(y, copy=True)
        mol.data['OQP::nac_xstate'] = np.array(Xf[:, J], copy=True)
        oqp.mrsf_nac_wpair(mol, I + 1, J + 1)
        return np.array(
            mol.data['OQP::nac_mt_frozen'], copy=True
        ).reshape(nbf, nbf).T

    occ_a = np.zeros(nbf)
    occ_a[:noca] = 1.0
    occ_b = np.zeros(nbf)
    occ_b[:nocb] = 1.0

    def space_of(iorb):
        if iorb < nocb:
            return 0
        if iorb < noca:
            return 1
        return 2

    def symmetric_u_contraction(matrix):
        """Contract a response matrix with U_sym=-Sx/2 plus OV elimination."""
        out = np.zeros(ncoord)
        for c in range(ncoord):
            value = 0.0
            for p in range(nbf):
                for q in range(p):
                    if space_of(p) == space_of(q):
                        value -= 0.5 * (
                            matrix[p, q] + matrix[q, p]
                        ) * Sx_MO[c, p, q]
                    else:
                        hi, lo = (p, q) if space_of(p) > space_of(q) else (q, p)
                        value -= matrix[lo, hi] * Sx_MO[c, hi, lo]
            out[c] = value
        return out

    ytil = {
        (I, J): ytil_of(I, J)
        for I in range(nstate)
        for J in range(nstate)
        if I != J
    }
    mt_frozen = {
        (I, J): frozen_orbital_gradient(I, J, ytil[I, J])
        for I in range(nstate)
        for J in range(nstate)
        if I != J
    }

    dp = np.zeros((ncoord, nstate, nstate))
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            y = ytil[I, J]
            # T1 skeleton (slot injection)
            inject(I, y)
            oqp.mrsf_nac_amp(mol)
            a = np.array(mol.data['OQP::nac_amp'], copy=True
                         ).reshape((nstate, nstate, natom, 3))
            oqp.mrsf_nac_esum(mol, I + 1, J + 1)
            es = np.array(mol.data['OQP::nac_esum'], copy=True).reshape(-1)

            # Capture the interstate density while the ordered pair is active.
            pij_a = np.array(
                mol.data['OQP::dbg_pij_a'], copy=True
            ).reshape(nbf, nbf).T
            pij_b = np.array(
                mol.data['OQP::dbg_pij_b'], copy=True
            ).reshape(nbf, nbf).T
            mol.data['OQP::td_bvec_mo'] = X0_raw
            T1 = a[J, I].reshape(-1) + es
            mt = mt_frozen[I, J]

            # Full ground-state Fock response of the interstate density.
            mol.data['OQP::nac_dm1_a'] = pack_sym(pij_a)
            mol.data['OQP::nac_dm1_b'] = pack_sym(pij_b)
            oqp.mrsf_nac_response(mol)
            ga = unpack_sym(np.array(
                mol.data['OQP::nac_v1_a'], copy=True
            ).reshape(-1))
            gb = unpack_sym(np.array(
                mol.data['OQP::nac_v1_b'], copy=True
            ).reshape(-1))
            gma = C.T @ ga @ C
            gmb = C.T @ gb @ C
            mtg = 2.0 * (
                gma * occ_a[None, :] + gmb * occ_b[None, :]
            )

            xmat = mt + mtg + gam[I, J]

            # Direct-injection seam: seam(X) = -X:U_cross.
            Lmat = xmat
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
            seam = gZ - gS
            gsk = np.array([float(np.sum(gam[I, J] * Sk_an[c]))
                            for c in range(ncoord)])
            vmask = symmetric_u_contraction(xmat)
            dp[:, I, J] = T1 - seam + vmask + gsk
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
