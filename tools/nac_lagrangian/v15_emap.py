"""v15 (FINAL MEASUREMENT): entrywise h-scaling of the reconstruction
error e(h) = W_displaced_signfixed - (1 + h u^T) W0. Entries whose
ratio e(h/2)/e(h) is ~1 (not ~1/4) scale as h^1 -> they name the
missing direction. Also: the response check A(C0 + e-only).
Run: python v15_emap.py <inp> <coord>
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
    c0 = int(sys.argv[2])
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v15.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
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
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    def matvec(v):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = v
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        oqp.mrsf_matvec_apply(mol)
        return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()[:nij].copy()

    J = 2
    XJ = Xf[:, J]
    base = matvec(XJ)
    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A', 'OQP::E_MO_B']}
    mol.save_data()
    cfg = mol.config
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = mol.log.replace('.log', '.json')
    cfg['guess']['continue_geom'] = False

    def restore():
        mol.update_system(xyz0)
        oqp.library.ints_1e(mol)
        for k, v in SAVE0.items():
            mol.data[k] = v.copy()
        mol.data['OQP::td_bvec_mo'] = X0_raw
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass

    def get_e(h):
        cp = xyz0.copy(); cp[c0] += h
        mol.update_system(cp)
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
        M_f = M_f * sg[None, :]
        Wd = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
        Wd_sf = sg[:, None] * Wd
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        oqp.get_structures_ao_overlap(mol)
        Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
        u = (M_f - Skf) / h
        W_rec = (np.eye(nbf) + h * u.T) @ W0
        e = Wd_sf - W_rec
        restore()
        return e, u, Wd_sf

    e1, u1, Wd1 = get_e(1e-3)
    e2, u2, Wd2 = get_e(5e-4)
    print('===== v15 entrywise scaling map (W-rows = MOs) =====')
    print(f'|e(h)|max={np.abs(e1).max():.3e} |e(h/2)|max={np.abs(e2).max():.3e}')
    # find entries NOT halving-by-4 (h^1 candidates): ratio e2/e1 ~ 0.5
    mask = np.abs(e1) > 1e-7
    ratio = np.where(mask, np.abs(e2) / np.maximum(np.abs(e1), 1e-30), 0.0)
    lin = (ratio > 0.35) & (ratio < 0.75) & mask       # ~h^1
    quad = (ratio <= 0.35) & mask                       # ~h^2
    print(f'entries |e|>1e-7: total={int(mask.sum())}, ~h^1={int(lin.sum())}, '
          f'~h^2={int(quad.sum())}')
    idx = np.argsort(-np.abs(e1) * lin.ravel()[np.argsort(-np.abs(e1).ravel())].astype(float)[0:1].sum() if False else -np.abs(e1 * lin).ravel())[:12]
    print('largest ~h^1 entries (MO-row p, AO-col mu):')
    for f in idx:
        p, mu = divmod(int(f), nbf)
        if not lin[p, mu]:
            continue
        print(f'  W[{p},{mu}] e(h)={e1[p,mu]:+.3e} e(h/2)={e2[p,mu]:+.3e} '
              f'ratio={ratio[p,mu]:.3f}  (MO {p}: '
              f'{"doc" if p<nocb else ("socc" if p<noca else "virt")})')
    # response attribution: perturb by the h^1 part only
    e_lin = np.where(lin, e1, 0.0)
    mol.data['OQP::VEC_MO_A'] = W0 + e_lin
    mol.data['OQP::VEC_MO_B'] = (W0 + e_lin).copy()
    a = matvec(XJ)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    resp = (a - base) / 1e-3
    print(f'\nresponse to h^1-part alone: |dA/h| = {np.linalg.norm(resp):.5f} '
          f'(target ~0.0196)')
    e_quad = np.where(quad, e1, 0.0)
    mol.data['OQP::VEC_MO_A'] = W0 + e_quad
    mol.data['OQP::VEC_MO_B'] = (W0 + e_quad).copy()
    a2 = matvec(XJ)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    print(f'response to h^2-part alone: |dA/h| = {np.linalg.norm((a2-base)/1e-3):.5f}')

if __name__ == '__main__':
    main()
