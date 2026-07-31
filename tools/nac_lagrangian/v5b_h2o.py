"""v5b: STATE-POLARIZATION OF THE FULL PRODUCTION GRADIENT CHAIN.

The gradient Lagrangian g(X) (z-vector + relaxed densities + W + all
contractions) is a QUADRATIC form in the target-slot amplitudes, so

  g_pol^c[I,J] = 1/2 [ g(ytil+X_J) - g(ytil) - g(X_J) ]^c   (target J)

is the complete interstate ytil^T (dA)_fully-relaxed X_J including every
orthonormality/W/G[dD] bookkeeping the production gradient already has.
Constant (SCF/nuclear) and linear parts cancel identically in the
polarization. The full coupling is then

  d^c[I,J] = g_pol^c + gamma:Sk_an^c + z_gamma.B^x^c
  (z_gamma: seam solve with RHS = gamma_a via the orbgrad hook)

Run: python v5b_h2o.py <input.inp> <dnum.npz>
"""
import os
import sys
import numpy as np

EPSA = 1e-5


def main():
    import oqp
    from oqp.pyoqp import Runner
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp, dnum_npz = sys.argv[1], sys.argv[2]
    dn_f = np.load(dnum_npz)
    dcv_n = dn_f['dcv' if 'dcv' in dn_f.files else dn_f.files[0]]

    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v5b.log'))
    r.run()
    mol = r.mol
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
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    C = W0.T
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

    gam = FK.gamma_closed(ctx)
    flatF = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            flatF[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = \
                gam[I, J].T.reshape(-1)
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

    print(f'building full A ({nij}x{nij})...', flush=True)
    A = np.zeros((nij, nij))
    for k in range(nij):
        e = np.zeros(nij)
        e[k] = 1.0
        A[:, k] = matvec(e)
    mol.data['OQP::td_bvec_mo'] = X0_raw
    phys = [k for k in range(nij) if k != ijlr2 - 1]
    Ap = 0.5 * (A[np.ix_(phys, phys)] + A[np.ix_(phys, phys)].T)
    w_full, Vp = np.linalg.eigh(Ap)
    V = np.zeros((nij, nij - 1))
    V[phys, :] = Vp
    ks = [int(np.argmax(np.abs(V.T @ Xf[:, s]))) for s in range(nstate)]

    sij0 = FK.s_ij_of(ctx, np.eye(nbf))
    sab0 = FK.s_ab_of(ctx, np.eye(nbf))
    sia0 = FK.s_ia_of(ctx, np.eye(nbf))

    def ampdir_unf(J, dXt):
        Xp = [x.copy() for x in Xt0]
        Xm = [x.copy() for x in Xt0]
        Xp[J] = Xt0[J] + EPSA * dXt
        Xm[J] = Xt0[J] - EPSA * dXt
        Sp = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xp)
        Sm = FK.contraction(ctx, sij0, sab0, sia0, Xt0, Xm)
        return (Sp[:, J] - Sm[:, J]) / (2 * EPSA)

    print('G_met sweep...', flush=True)
    G_met = np.zeros((nstate, nstate, nij))
    for J in range(nstate):
        for k in range(nij):
            if k == ijlr2 - 1:
                continue
            e = np.zeros(nij)
            e[k] = 1.0
            g = ampdir_unf(J, unfold_vec(e))
            for I in range(nstate):
                if I != J:
                    G_met[I, J, k] = g[I]
    ytil = np.zeros((nstate, nstate, nij))
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            coef = V.T @ G_met[I, J]
            acc = np.zeros(nij)
            for k in range(nij - 1):
                if k == ks[J]:
                    continue
                den = Om[J] - w_full[k]
                if abs(den) < 1e-9:
                    continue
                acc += V[:, k] * (coef[k] / den)
            ytil[I, J] = acc

    # ---- full production gradient as a function of the target slot -------
    def grad_full(J, vec):
        rr = X0_raw.copy().reshape(-1)
        rr[J * nij:(J + 1) * nij] = vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        mol.data.set_tdhf_target(J + 1)
        oqp.tdhf_mrsf_z_vector(mol)
        conv = bool(mol.mol_energy.Z_Vector_converged)
        oqp.tdhf_mrsf_gradient(mol)
        g = np.array(mol.get_grad(), copy=True).reshape(-1)
        mol.data['OQP::td_bvec_mo'] = X0_raw
        return g, conv

    print('state-polarized gradient runs...', flush=True)
    g0, _ = grad_full(0, np.zeros(nij))     # constant (SCF/nuclear) part
    g_pol = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            gm, c1 = grad_full(J, ytil[I, J] + Xf[:, J])
            gy, c2 = grad_full(J, ytil[I, J])
            gx, c3 = grad_full(J, Xf[:, J])
            if not (c1 and c2 and c3):
                print(f'  ({I+1},{J+1}): z unconverged {c1} {c2} {c3}',
                      flush=True)
            g_pol[(I, J)] = 0.5 * (gm - gy - gx + g0)
            print(f'  ({I+1},{J+1}): |g_pol|={np.linalg.norm(g_pol[(I,J)]):.6f}',
                  flush=True)

    # ---- z_gamma seam (gamma_a alone through the hook) -------------------
    def unpack_dummy():
        pass

    def grad_now():
        oqp.tdhf_mrsf_gradient(mol)
        return np.array(mol.get_grad(), copy=True).reshape(-1)

    zg = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            mol.data['OQP::nac_orbgrad_L'] = gam[I, J].ravel(order='F').copy()
            mol.data['OQP::td_bvec_mo'] = X0_raw
            mol.data.set_tdhf_target(J + 1)
            oqp.set_mrsf_nac_cphf(mol, I + 1, J + 1)
            oqp.tdhf_mrsf_z_vector(mol)
            if not bool(mol.mol_energy.Z_Vector_converged):
                oqp.set_mrsf_nac_cphf(mol, 0, 0)
                zg[(I, J)] = None
                continue
            gZ = grad_now()
            mol.data['OQP::td_p'] = np.zeros_like(
                np.array(mol.data['OQP::td_p'], copy=True))
            mol.data['OQP::WAO'] = np.zeros_like(
                np.array(mol.data['OQP::WAO'], copy=True))
            gS = grad_now()
            oqp.set_mrsf_nac_cphf(mol, 0, 0)
            zg[(I, J)] = gZ - gS
    try:
        del mol.data['OQP::nac_orbgrad_L']
    except Exception:
        pass
    mol.data['OQP::td_bvec_mo'] = X0_raw
    print('z_gamma seams done.', flush=True)

    print('\n===== v5b POLARIZED-GRADIENT ASSEMBLY vs d_num =====')
    dp = np.zeros((ncoord, nstate, nstate))
    for I in range(nstate):
        for J in range(nstate):
            if I == J or zg[(I, J)] is None:
                continue
            gsk = np.array([float(np.sum(gam[I, J] * Sk_an[c]))
                            for c in range(ncoord)])
            dp[:, I, J] = g_pol[(I, J)] + zg[(I, J)] + gsk
    dpa = 0.5 * (dp - dp.transpose(0, 2, 1))
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            dn = dcv_n[I, J].reshape(-1)
            for tag, v in (('direct ', dp[:, I, J]), ('antisym', dpa[:, I, J])):
                cc = float(np.dot(dn, v)) / (np.linalg.norm(dn)
                                             * np.linalg.norm(v) + 1e-300)
                md = min(np.abs(v - dn).max(), np.abs(v + dn).max())
                print(f'pair ({I+1},{J+1}) [{tag}]: |d_num|={np.linalg.norm(dn):.6f} '
                      f'|pred|={np.linalg.norm(v):.6f} cos={cc:+.8f} '
                      f'sign-resolved maxdiff={md:.3e}')
    np.savez(inp.replace('.inp', '_v5b.npz'), dp=dp, dpa=dpa)


if __name__ == '__main__':
    main()
