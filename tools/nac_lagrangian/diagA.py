"""Diagnose the matvec-PT failure: eigen-identity checks for mrsf_matvec_apply."""
import sys
import numpy as np


def main():
    import oqp
    import oqp.library
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint

    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_dA.log'))
    r.run()
    mol = r.mol

    nstate = mol.config['tdhf']['nstate']
    noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
    nocb = noca - 2
    nbf = np.array(mol.data['OQP::VEC_MO_A']).shape[0]
    nvirb = nbf - nocb
    nij = noca * nvirb

    X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
    Xshape = X0_raw.shape
    X0 = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    xyz0 = np.array(mol.get_system(), copy=True)
    E = list(mol.energies)
    Om = np.array([E[k + 1] - E[0] for k in range(nstate)])

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

    def apply_A(cols):
        out = np.zeros((nij, cols.shape[1]))
        for s in range(cols.shape[1]):
            set_trial(cols[:, s])
            oqp.mrsf_matvec_apply(mol)
            out[:, s] = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
        return out

    np.set_printoptions(precision=8, suppress=True)
    print('Omega =', Om)
    print('X0 norms:', [f'{np.linalg.norm(X0[:, s]):.6f}' for s in range(nstate)])

    AX0 = apply_A(X0)
    G = X0.T @ AX0
    print('\nX0^T A(x0) X0 (should be diag(Omega)):')
    print(G)
    print('asymmetry |G - G^T|max =', np.abs(G - G.T).max())

    # displaced check
    mol.data['OQP::td_bvec_mo'] = X0_raw
    mol.save_data()
    cfg = mol.config
    json0 = mol.log.replace('.log', '.json')
    cfg['guess']['type'] = 'json'
    cfg['guess']['file'] = json0
    cfg['guess']['continue_geom'] = False
    h = 1e-3
    coord = xyz0.copy()
    coord[3] += h
    mol.update_system(coord)
    oqp.library.ints_1e(mol)
    oqp.library.guess(mol)
    SinglePoint(mol).energy()
    Ed = list(mol.energies)
    Omd = np.array([Ed[k + 1] - Ed[0] for k in range(nstate)])
    Xd_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
    Xd = Xd_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    AXd = apply_A(Xd)
    Gd = Xd.T @ AXd
    print('\nOmega(+h) =', Omd)
    print('Xd^T A(x+h) Xd (should be diag(Omega(+h))):')
    print(Gd)
    resid = AXd - Xd * Omd[None, :]
    print('eigen-residual |A Xd - Om Xd| per state:',
          [f'{np.linalg.norm(resid[:, s]):.2e}' for s in range(nstate)])

    # dA numerator symmetry
    AX0d = apply_A(X0)          # A(x+h) on reference vectors
    N = X0.T @ AX0d
    print('\nX0^T A(x+h) X0:')
    print(N)
    print('asymmetry |N - N^T|max =', np.abs(N - N.T).max())

    # WHO IS LYING: basis proximity and operator linearity
    print('\nX0^T Xd (state-overlap of folded vectors):')
    print(X0.T @ Xd)
    v = X0[:, 0]
    A1 = apply_A(v.reshape(-1, 1))[:, 0]
    A2 = apply_A((2.0 * v).reshape(-1, 1))[:, 0]
    print('linearity |A(2v) - 2A(v)| / |A(v)| =',
          np.linalg.norm(A2 - 2 * A1) / (np.linalg.norm(A1) + 1e-300))
    vv = 0.5 * (X0[:, 0] + X0[:, 1]) * np.sqrt(2)
    A3 = apply_A(vv.reshape(-1, 1))[:, 0]
    A12 = apply_A(X0[:, :2])
    lin2 = np.linalg.norm(A3 - np.sqrt(2) * 0.5 * (A12[:, 0] + A12[:, 1]))
    print('additivity |A(a+b) - A(a) - A(b)| (scaled) =', lin2)


if __name__ == '__main__':
    main()
