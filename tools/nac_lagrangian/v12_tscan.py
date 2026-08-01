"""v12: t-scan of the staged response along u+ (one coordinate).
staged(t) = [A(C0(1+t u+)) - A(C0)] X / t   for t in {1e-5..1e-3}
vs the one-sided FD lhs+ = [A(C+) - A(C0, INTS+)]... NOTE: lhs+ mixes
displaced integrals; the pure-C comparison at REFERENCE integrals is
staged(t=h) vs staged(t->0) -- the t-scan itself is the verdict:
strong nonlinearity between 1e-5 and 1e-3 = the hidden premise.
Run: python v12_tscan.py <inp> <coord>
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
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v12.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']
    ctx = FK.build_context(mol)
    nbf, nij = ctx['nbf'], ctx['nij']
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

    # u+ from one displaced run
    mol.save_data()
    cfg = mol.config
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = mol.log.replace('.log', '.json')
    cfg['guess']['continue_geom'] = False
    HD = 1e-3
    SAVE0 = {k: np.array(mol.data[k], copy=True) for k in
             ['OQP::DM_A', 'OQP::DM_B', 'OQP::FOCK_A', 'OQP::FOCK_B',
              'OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::E_MO_A', 'OQP::E_MO_B']}
    cp = xyz0.copy(); cp[c0] += HD
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
    mol.data['OQP::VEC_MO_A'] = W0
    mol.data['OQP::VEC_MO_B'] = Wb0
    oqp.get_structures_ao_overlap(mol)
    Sk_np = np.array(mol.data['OQP::overlap_mo_non_orthogonal'], copy=True)
    Skf = Sk_np.reshape(-1).reshape((nbf, nbf)).T
    up = (M_f - Skf) / HD
    mol.update_system(xyz0)
    oqp.library.ints_1e(mol)
    for k, v in SAVE0.items():
        mol.data[k] = v.copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    print('===== t-scan of staged response along u+ =====')
    ref = None
    for t in (1e-5, 1e-4, 5e-4, 1e-3):
        Wp = (np.eye(nbf) + t * up.T) @ W0
        mol.data['OQP::VEC_MO_A'] = Wp
        mol.data['OQP::VEC_MO_B'] = Wp
        a = matvec(XJ)
        mol.data['OQP::VEC_MO_A'] = W0
        mol.data['OQP::VEC_MO_B'] = Wb0
        st = (a - base) / t
        if ref is None:
            ref = st
            print(f'  t={t:.0e}: |staged|={np.linalg.norm(st):.5f} (reference)')
        else:
            print(f'  t={t:.0e}: |staged|={np.linalg.norm(st):.5f} '
                  f'|staged - staged(1e-5)|={np.linalg.norm(st - ref):.5f} '
                  f'maxdiff={np.abs(st - ref).max():.4e}')

if __name__ == '__main__':
    main()
