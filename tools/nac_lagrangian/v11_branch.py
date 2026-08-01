"""v11: branch/one-sided forensics at one coordinate.
 (B1) u+ vs -u- blockwise (branch kink check)
 (B2) one-sided FD [A(C+)-A(C0)]/h  vs  staged along u+   (linearity on
      the actual branch; also vs staged along Ux for comparison)
Run: python v11_branch.py <inp> <coord>
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
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v11.log'))
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
    HD = 1e-3

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
        M_f = M_f * sg[None, :]
        try:
            mol.data._data.control.int2e_cutoff = 1e-20
        except Exception:
            pass
        Wd = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
        Wd = sg[:, None] * Wd
        mol.data['OQP::VEC_MO_A'] = Wd
        mol.data['OQP::VEC_MO_B'] = Wd.copy()
        Ax = matvec(XJ)
        # frozen-C at the same geometry
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        Ax0 = matvec(XJ)
        oqp.get_structures_ao_overlap(mol)
        Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
        Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
        return M_f, Skf, Ax, Ax0

    cp = xyz0.copy(); cp[c0] += HD
    Mp, Skp, Axp, Ax0p = disp(cp)
    cm = xyz0.copy(); cm[c0] -= HD
    Mm, Skm, Axm, Ax0m = disp(cm)
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    up = (Mp - Skp * 0 - np.eye(nbf)) / HD - (Skp - np.eye(nbf)) / HD
    um = (Mm - np.eye(nbf)) / (-HD) - (Skm - np.eye(nbf)) / (-HD)
    Ux = 0.5 * (up + um)

    def space_of(i):
        return 0 if i < nocb else (1 if i < noca else 2)

    print('===== B1: u+ vs u- blockwise =====')
    dU = np.abs(up - um)
    blocks = {}
    for p in range(nbf):
        for q in range(nbf):
            key = f'{space_of(p)}{space_of(q)}'
            blocks.setdefault(key, 0.0)
            blocks[key] = max(blocks[key], dU[p, q])
    print('  max|u+ - u-| per block (00=dd,11=ss,22=vv):', 
          {k: round(v, 4) for k, v in sorted(blocks.items())})
    print(f'  |u+|max={np.abs(up).max():.4f} |Ux|max={np.abs(Ux).max():.4f}')
    idx = np.argsort(-dU.ravel())[:6]
    for f in idx:
        p, q = divmod(int(f), nbf)
        print(f'  worst ({p},{q}) [{space_of(p)}{space_of(q)}]: '
              f'u+={up[p,q]:+.4f} u-={um[p,q]:+.4f}')

    # B2: one-sided response
    lhs_p = (Axp - Ax0p) / HD          # one-sided C-channel on the +h branch
    t = 1e-5
    Wp1 = (np.eye(nbf) + t * up.T) @ W0
    mol.data['OQP::VEC_MO_A'] = Wp1
    mol.data['OQP::VEC_MO_B'] = Wp1
    a_p = matvec(XJ)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    staged_up = (a_p - base) / t
    WpX = (np.eye(nbf) + t * Ux.T) @ W0
    mol.data['OQP::VEC_MO_A'] = WpX
    mol.data['OQP::VEC_MO_B'] = WpX
    a_x = matvec(XJ)
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    staged_ux = (a_x - base) / t
    print('\\n===== B2: one-sided FD vs staged =====')
    print(f'  |lhs+|={np.linalg.norm(lhs_p):.5f} |staged(u+)|={np.linalg.norm(staged_up):.5f} '
          f'|staged(Ux)|={np.linalg.norm(staged_ux):.5f}')
    print(f'  |lhs+ - staged(u+)| = {np.linalg.norm(lhs_p - staged_up):.5f} '
          f'maxdiff={np.abs(lhs_p - staged_up).max():.4e}')
    print(f'  |lhs+ - staged(Ux)| = {np.linalg.norm(lhs_p - staged_ux):.5f}')

if __name__ == '__main__':
    main()
