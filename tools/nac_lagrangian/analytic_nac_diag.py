"""DIAGNOSTIC for the production-candidate assembly: dump every term
separately (single upper orientation, no JI mixing), check the symmetry
contracts, and least-squares-fit term coefficients against d_num to
localize wrong factors.  Ideal coefficients = [1,1,1,1,1] for
[amp2e/gap, esum/gap, zB(L)/gap, zB(gap*gam)/gap, ov].
"""
import os
import sys
import numpy as np


def main():
    import oqp
    from oqp.pyoqp import Runner
    from oqp.library.single_point import NAC
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_dg.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ncoord = 3 * natom
    E = list(mol.energies)
    Om = [E[k + 1] - E[0] for k in range(nstate)]

    ctx = FK.build_context(mol)
    noca, nocb, nbf, nij = ctx['noca'], ctx['nocb'], ctx['nbf'], ctx['nij']
    X0_raw = ctx['X0_raw']

    nac = NAC(mol)
    nacv_n, dcv_n, flags = nac.numerical_nac()
    print('numerical flags:', set(flags), flush=True)
    mol.data['OQP::td_bvec_mo'] = X0_raw

    oqp.mrsf_nac_amp(mol)
    amp2e = np.array(mol.data['OQP::nac_amp'], copy=True
                     ).reshape((nstate, nstate, natom, 3))
    esum_v = {}
    for i in range(1, nstate + 1):
        for j in range(1, nstate + 1):
            if i == j:
                continue
            oqp.mrsf_nac_esum(mol, i, j)
            esum_v[(i, j)] = np.array(mol.data['OQP::nac_esum'],
                                      copy=True).reshape(-1)

    gam = FK.gamma_closed(ctx)
    # symmetry contracts of gamma across pair orientation
    for I in range(nstate):
        for J in range(I + 1, nstate):
            a = gam[I, J]; b = gam[J, I]
            print(f'gam contract ({I+1},{J+1}): |g_IJ+g_JI|={np.abs(a+b).max():.2e} '
                  f'|g_IJ-g_JI.T|={np.abs(a-b.T).max():.2e} '
                  f'|g_IJ+g_IJ.T|={np.abs(a+a.T).max():.2e}', flush=True)

    flatF = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            flatF[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = \
                gam[I, J].T.reshape(-1)
    try:
        del mol.data['OQP::nac_gamma_tlf']
    except Exception:
        pass
    mol.data['OQP::nac_gamma_tlf'] = flatF.reshape((nbf * nbf, nstate, nstate))
    oqp.mrsf_nac_overlap(mol)
    ovF = np.array(mol.data['OQP::nac_overlap'], copy=True
                   ).reshape(-1).reshape((nstate, nstate, ncoord))
    ov = np.transpose(ovF, (1, 0, 2))

    os.environ['NAC_DUMP_RHS'] = '1'

    def rhs_of(target, amp_vec):
        rr = X0_raw.copy().reshape(-1)
        rr[(target - 1) * nij:target * nij] = amp_vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        mol.data.set_tdhf_target(target)
        oqp.tdhf_mrsf_z_vector(mol)
        out = np.array(mol.data['OQP::nac_zvec_rhs'], copy=True).reshape(-1)
        mol.data['OQP::td_bvec_mo'] = X0_raw
        return out

    def unpack_rot_to_mat(v):
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

    def zB_of(Lmat, target):
        mol.data['OQP::nac_orbgrad_L'] = Lmat.ravel(order='F').copy()
        mol.data['OQP::td_bvec_mo'] = X0_raw
        mol.data.set_tdhf_target(target)
        oqp.set_mrsf_nac_cphf(mol, target, target % nstate + 1)
        oqp.tdhf_mrsf_z_vector(mol)
        ok = bool(mol.mol_energy.Z_Vector_converged)
        if not ok:
            oqp.set_mrsf_nac_cphf(mol, 0, 0)
            return None
        gZ = grad_now()
        mol.data['OQP::td_p'] = np.zeros_like(
            np.array(mol.data['OQP::td_p'], copy=True))
        mol.data['OQP::WAO'] = np.zeros_like(
            np.array(mol.data['OQP::WAO'], copy=True))
        gS = grad_now()
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        return gZ - gS

    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    rhs_single = {s: rhs_of(s + 1, Xf[:, s]) for s in range(nstate)}

    print('\n===== TERM DIAGNOSTIC (upper orientation only) =====', flush=True)
    for I in range(nstate):
        for J in range(I + 1, nstate):
            gap = Om[J] - Om[I]
            dn = dcv_n[I, J].reshape(-1)
            Lpol = 0.5 * (rhs_of(I + 1, Xf[:, I] + Xf[:, J])
                          - rhs_single[I] - rhs_of(I + 1, Xf[:, J]))
            zL = zB_of(unpack_rot_to_mat(Lpol), I + 1)
            zG = zB_of(gap * gam[I, J], I + 1)
            if zL is None or zG is None:
                print(f'pair ({I+1},{J+1}): z NOT converged', flush=True)
                continue
            t_amp = amp2e[J, I].reshape(-1) / gap
            t_ampT = amp2e[I, J].reshape(-1) / gap
            t_es = esum_v[(I + 1, J + 1)] / gap
            t_esT = esum_v[(J + 1, I + 1)] / gap
            t_zL = zL / gap
            t_zG = zG / gap
            t_ov = ov[I, J]
            t_ovT = ov[J, I]
            print(f'\npair ({I+1},{J+1}) gap={gap:.6f} |d_num|={np.linalg.norm(dn):.6f}')
            print(f'  sym checks: |amp-ampT|={np.abs(t_amp-t_ampT).max():.2e} '
                  f'|es-esT|={np.abs(t_es-t_esT).max():.2e} '
                  f'|ov+ovT|={np.abs(t_ov+t_ovT).max():.2e} '
                  f'|ov-ovT|={np.abs(t_ov-t_ovT).max():.2e}')
            terms = {'amp': t_amp, 'esum': t_es, 'zL': t_zL,
                     'zG': t_zG, 'ov': t_ov}
            for k, v in terms.items():
                c = float(np.dot(dn, v)) / (np.linalg.norm(dn)
                                            * np.linalg.norm(v) + 1e-300)
                print(f'  {k:5s}: |v|={np.linalg.norm(v):.6f} cos(d_num)={c:+.4f}')
            A = np.stack(list(terms.values()), axis=1)
            coef, res, *_ = np.linalg.lstsq(A, dn, rcond=None)
            fit = A @ coef
            cc = float(np.dot(dn, fit)) / (np.linalg.norm(dn)
                                           * np.linalg.norm(fit) + 1e-300)
            print(f'  LSQ coef [amp,esum,zL,zG,ov] = '
                  + np.array2string(coef, precision=4, suppress_small=True)
                  + f'  fit-cos={cc:+.6f} resid={np.linalg.norm(dn-fit):.3e}')
            s1 = t_amp + t_es + t_zL + t_zG + t_ov
            print(f'  sum[1,1,1,1,1]: |v|={np.linalg.norm(s1):.6f} '
                  f'cos={float(np.dot(dn,s1))/(np.linalg.norm(dn)*np.linalg.norm(s1)+1e-300):+.6f} '
                  f'resid={np.linalg.norm(dn-s1):.3e}', flush=True)


if __name__ == '__main__':
    main()
