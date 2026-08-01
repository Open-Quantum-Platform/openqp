"""Fetch phase-anchored context for the ETH v7o npz: gamma, Sk_an, Sx,
Xf + a one-coordinate w-probe for phase alignment with the saved npz.
Run: python eth_ctx_fetch.py <inp>
"""
import os, sys
import numpy as np

def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK
    inp = sys.argv[1]
    os.environ['NAC_DUMP_DS'] = '1'
    r = Runner(input_file=inp, log=inp.replace('.inp', '_ctx.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    natom = mol.data['natom']
    ncoord = 3 * natom
    ctx = FK.build_context(mol)
    nbf, nij = ctx['nbf'], ctx['nij']
    noca, nocb = ctx['noca'], ctx['nocb']
    X0_raw = ctx['X0_raw']
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    W0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
    Wb0 = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
    e0a = np.array(mol.data['OQP::E_MO_A'], copy=True)
    e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
    xyz0 = np.array(mol.get_system(), copy=True)
    C = W0.T
    gam = FK.gamma_closed(ctx)
    flatF = np.zeros(nbf * nbf * nstate * nstate)
    for I in range(nstate):
        for J in range(nstate):
            flatF[(I + J * nstate) * nbf * nbf:(I + J * nstate + 1) * nbf * nbf] = \
                gam[I, J].T.reshape(-1)
    mol.data['OQP::nac_gamma_tlf'] = flatF.reshape((nbf * nbf, nstate, nstate))
    oqp.mrsf_nac_overlap(mol)
    nbf2 = nbf * nbf
    dsk = np.array(mol.data['OQP::dbg_dsket'], copy=True).reshape(-1)
    dsf = np.array(mol.data['OQP::dbg_dsfull'], copy=True).reshape(-1)
    Sk_an = np.zeros((ncoord, nbf, nbf))
    Sx = np.zeros((ncoord, nbf, nbf))
    for c in range(ncoord):
        Sk_an[c] = C.T @ dsk[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T @ C
        Sx[c] = C.T @ dsf[c * nbf2:(c + 1) * nbf2].reshape(nbf, nbf).T @ C
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    # one-coordinate displaced w probe (for phase anchoring)
    RS = ctx['RS']; nvirb = ctx['nvirb']
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

    def refold(x):
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
        return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()[:nij].copy()

    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False
    HD = 1.0e-3

    def disp(coord):
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
        sg = np.sign(np.diag(M_f))
        sg[sg == 0] = 1.0
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        Ax = np.zeros((nstate, nij))
        for s in range(nstate):
            Ax[s] = sg_apply(matvec(sg_apply(Xf[:, s], sg)), sg)
        return Ax

    c0 = 0
    cp = xyz0.copy(); cp[c0] += HD
    Axp = disp(cp)
    cm = xyz0.copy(); cm[c0] -= HD
    Axm = disp(cm)
    wprobe = (Axp - Axm) / (2 * HD)
    np.savez(inp.replace('.inp', '_ctx.npz'), gam=gam, Sk_an=Sk_an, Sx=Sx,
             Xf=Xf, noca=np.array(noca), nocb=np.array(nocb),
             wprobe=wprobe, c0=np.array(c0))
    print('ctx saved.', flush=True)

if __name__ == '__main__':
    main()
