"""v7i: (1) DIRECT-INJECTION seams with X = Mt + Mt_G + gamma from the
v7h npz (production T2 form); (2) two-step displaced sweep (h, h/2) with
Richardson extrapolation to discriminate FD-referee precision from a
missing Mt channel; (3) the production-form assembly

  d = T1 - seam(X) + t_elim(X) + t_ss_sym(Mt) + gamma:Sk

Run: python v7i_h2o.py <input.inp> <v7h.npz> <dnum.npz>
"""
import os
import sys
import numpy as np


def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp, h_npz, dnum_npz = sys.argv[1], sys.argv[2], sys.argv[3]
    H = np.load(h_npz)
    dn = np.load(dnum_npz)
    dcv = dn['dcv']
    gam = H['gam']
    Sk_an, Sx = H['Sk_an'], H['Sx_MO']
    noca, nocb = int(H['noca']), int(H['nocb'])

    r = Runner(input_file=inp, log=inp.replace('.inp', '_v7i.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ncoord = 3 * natom

    ctx = FK.build_context(mol)
    nbf, nij = ctx['nbf'], ctx['nij']
    X0_raw = ctx['X0_raw']
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a_r = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b_r = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    def space_of(i):
        return 0 if i < nocb else (1 if i < noca else 2)

    # ---- direct-injection seams ------------------------------------------
    def grad_now():
        oqp.tdhf_mrsf_gradient(mol)
        return np.array(mol.get_grad(), copy=True).reshape(-1)

    def seam_of(Lmat, I, J):
        mol.data['OQP::nac_orbgrad_L'] = Lmat.ravel(order='F').copy()
        mol.data['OQP::td_bvec_mo'] = X0_raw
        mol.data.set_tdhf_target(J + 1)
        oqp.set_mrsf_nac_cphf(mol, I + 1, J + 1)
        oqp.tdhf_mrsf_z_vector(mol)
        conv = bool(mol.mol_energy.Z_Vector_converged)
        gZ = grad_now()
        mol.data['OQP::td_p'] = np.zeros_like(
            np.array(mol.data['OQP::td_p'], copy=True))
        mol.data['OQP::WAO'] = np.zeros_like(
            np.array(mol.data['OQP::WAO'], copy=True))
        gS = grad_now()
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        return gZ - gS, conv

    print('direct-injection seams...', flush=True)
    T2d = {}
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            X = H[f'MT_{I}{J}'] + H[f'MTG_{I}{J}'] + gam[I, J]
            s, conv = seam_of(X, I, J)
            T2d[(I, J)] = s
            print(f'  ({I+1},{J+1}) conv={conv} |seam|={np.linalg.norm(s):.6f}',
                  flush=True)
    try:
        del mol.data['OQP::nac_orbgrad_L']
    except Exception:
        pass

    # ---- two-step displaced sweep ----------------------------------------
    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False

    def unfold_vec(v):
        RS = ctx['RS']
        nvirb = ctx['nvirb']
        ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
        ijlr2 = (noca - nocb - 1) * noca + noca
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

    def refold(x):
        RS = ctx['RS']
        nvirb = ctx['nvirb']
        ijlr1 = (noca - 1 - nocb - 1) * noca + noca - 1
        ijlr2 = (noca - nocb - 1) * noca + noca
        v = np.zeros(nij)
        for i in range(1, noca + 1):
            for jj in range(nocb + 1, nbf + 1):
                ij = (jj - nocb - 1) * noca + i
                if ij == ijlr1:
                    v[ijlr1 - 1] = x[i - 1, jj - nocb - 1] / RS
                elif ij == ijlr2:
                    pass
                else:
                    v[ij - 1] = x[i - 1, jj - nocb - 1]
        return v

    def sg_apply(v, sg):
        x = unfold_vec(v)
        x = sg[:noca, None] * x * sg[None, nocb:]
        return refold(x)

    def matvec(v):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = v
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        oqp.mrsf_matvec_apply(mol)
        return np.array(mol.data['OQP::nac_mvax'],
                        copy=True).ravel()[:nij].copy()

    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A',
              'OQP::E_MO_B']}

    def displaced_all(coord):
        mol.update_system(coord)
        oqp.library.ints_1e(mol)
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        mol.data['OQP::xyz_old'] = xyz0.reshape((3, -1))
        mol.data['OQP::VEC_MO_A_old'] = W0
        mol.data['OQP::VEC_MO_B_old'] = Wb0
        mol.data['OQP::E_MO_A_old'] = e0a_r
        mol.data['OQP::E_MO_B_old'] = e0b_r
        oqp.get_structures_ao_overlap(mol)
        M_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        M_f = M_np.reshape(-1).reshape((nbf, nbf)).T
        sg = np.sign(np.diag(M_f))
        sg[sg == 0] = 1.0
        M_f = M_f * sg[None, :]
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        Ax_ref = np.zeros((nstate, nij))
        for s in range(nstate):
            Ax_ref[s] = sg_apply(matvec(sg_apply(Xf[:, s], sg)), sg)
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        oqp.get_structures_ao_overlap(mol)
        Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
        return M_f, Skf, Ax_ref

    res = {}
    for hh in (1.0e-3, 5.0e-4):
        print(f'displaced sweep h={hh}...', flush=True)
        Ux = np.zeros((ncoord, nbf, nbf))
        w_ref = np.zeros((ncoord, nstate, nij))
        for c in range(ncoord):
            cp = xyz0.copy(); cp[c] += hh
            Mp, Skp, Axp = displaced_all(cp)
            cm = xyz0.copy(); cm[c] -= hh
            Mm, Skm, Axm = displaced_all(cm)
            Ux[c] = (Mp - Mm) / (2 * hh) - (Skp - Skm) / (2 * hh)
            w_ref[c] = (Axp - Axm) / (2 * hh)
        res[hh] = (Ux, w_ref)
        mol.update_system(xyz0)
        oqp.library.ints_1e(mol)
        for k, v in SAVE0.items():
            mol.data[k] = v.copy()
        mol.data['OQP::td_bvec_mo'] = X0_raw
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass

    Ux1, w1 = res[1.0e-3]
    Ux2, w2 = res[5.0e-4]
    UxR = (4.0 * Ux2 - Ux1) / 3.0
    wR = (4.0 * w2 - w1) / 3.0

    # ---- verdicts ---------------------------------------------------------
    print('\n===== v7i VERDICTS =====')
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            X = H[f'MT_{I}{J}'] + H[f'MTG_{I}{J}'] + gam[I, J]
            y = H[f'ytil_{I}{J}']
            T1 = H[f'T1_{I}{J}']
            # (i) -seam vs t_pack (both h and Richardson)
            for tag, U in (('h1 ', Ux1), ('Rich', UxR)):
                t_pack = np.zeros(ncoord)
                for c in range(ncoord):
                    a = 0.0
                    for p in range(nbf):
                        for q in range(p):
                            sp, sq = space_of(p), space_of(q)
                            if sp != sq:
                                hi, lo = (p, q) if sp > sq else (q, p)
                                a += (X[hi, lo] - X[lo, hi]) * U[c][hi, lo]
                    t_pack[c] = a
                md = np.abs(-T2d[(I, J)] - t_pack).max()
                if tag == 'Rich':
                    print(f'({I+1},{J+1}) -seam vs pack.U[{tag}]: '
                          f'|seam|={np.linalg.norm(T2d[(I,J)]):.5f} '
                          f'|pack|={np.linalg.norm(t_pack):.5f} maxdiff={md:.3e}')
            # (ii) Mt-completeness with the Richardson referee
            uch = np.array([float(np.sum(X * UxR[c])) for c in range(ncoord)])
            rhs = np.array([float(np.dot(y, wR[c, J]))
                            for c in range(ncoord)]) - T1
            print(f'({I+1},{J+1}) Mt:U[Rich] vs exact[Rich]: '
                  f'{np.linalg.norm(uch):.5f} / {np.linalg.norm(rhs):.5f} '
                  f'maxdiff={np.abs(uch-rhs).max():.3e}')

    print('\n===== PRODUCTION-FORM ASSEMBLY vs d_num =====')
    dp = np.zeros((ncoord, nstate, nstate))
    for I in range(nstate):
        for J in range(nstate):
            if I == J:
                continue
            Mt = H[f'MT_{I}{J}'] + H[f'MTG_{I}{J}']
            X = Mt + gam[I, J]
            T1 = H[f'T1_{I}{J}']
            t_elim = np.zeros(ncoord)
            t_ss = np.zeros(ncoord)
            for c in range(ncoord):
                e = s1 = 0.0
                for p in range(nbf):
                    for q in range(p):
                        sp, sq = space_of(p), space_of(q)
                        if sp != sq:
                            hi, lo = (p, q) if sp > sq else (q, p)
                            e += -X[lo, hi] * Sx[c][hi, lo]
                        else:
                            s1 += (Mt[p, q] + Mt[q, p]) * (-0.5) * Sx[c][p, q]
                t_elim[c] = e
                t_ss[c] = s1
            gsk = np.array([float(np.sum(gam[I, J] * Sk_an[c]))
                            for c in range(ncoord)])
            dp[:, I, J] = T1 - T2d[(I, J)] + t_elim + t_ss + gsk
    dpa = 0.5 * (dp - dp.transpose(0, 2, 1))
    for I in range(nstate):
        for J in range(nstate):
            if I >= J:
                continue
            v = dpa[:, I, J]
            ref = dcv[I, J].reshape(-1)
            cc = float(np.dot(ref, v)) / (np.linalg.norm(ref)
                                          * np.linalg.norm(v) + 1e-300)
            md = min(np.abs(v - ref).max(), np.abs(v + ref).max())
            print(f'pair ({I+1},{J+1}): |d_num|={np.linalg.norm(ref):.6f} '
                  f'|pred|={np.linalg.norm(v):.6f} cos={cc:+.8f} '
                  f'sign-res maxdiff={md:.3e}')
    np.savez(inp.replace('.inp', '_v7i.npz'), Ux1=Ux1, Ux2=Ux2, w1=w1, w2=w2,
             **{f'T2d_{I}{J}': T2d[(I, J)] for I in range(nstate)
                for J in range(nstate) if I != J})


if __name__ == '__main__':
    main()
