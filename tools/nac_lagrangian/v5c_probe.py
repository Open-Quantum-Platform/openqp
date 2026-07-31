"""v5c PROBE: locate the non-bilinear slot dependence in the production
z-vector + gradient chain.

  Q-test: purely-quadratic(+const) chain must satisfy
          g(2X) - 4 g(X) + 3 g(0) == 0
  Record test: after the z-solve, every density record must scale 4x
          from slot X -> slot 2X:  td_p, WAO, td_abxc, td_mrsf_density.

Run: python v5c_probe.py <input.inp>
"""
import os
import sys
import numpy as np


def main():
    import oqp
    from oqp.pyoqp import Runner
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK

    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v5c.log'))
    r.run()
    mol = r.mol
    nstate = mol.config['tdhf']['nstate']

    ctx = FK.build_context(mol)
    nij = ctx['nij']
    X0_raw = ctx['X0_raw']
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass

    RECS = ['OQP::td_p', 'OQP::WAO', 'OQP::td_abxc', 'OQP::td_mrsf_density']

    def chain(J, vec):
        rr = X0_raw.copy().reshape(-1)
        rr[J * nij:(J + 1) * nij] = vec
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        mol.data.set_tdhf_target(J + 1)
        oqp.tdhf_mrsf_z_vector(mol)
        conv = bool(mol.mol_energy.Z_Vector_converged)
        recs = {}
        for k in RECS:
            try:
                recs[k] = np.array(mol.data[k], copy=True).ravel().copy()
            except Exception:
                recs[k] = None
        oqp.tdhf_mrsf_gradient(mol)
        g = np.array(mol.get_grad(), copy=True).reshape(-1)
        mol.data['OQP::td_bvec_mo'] = X0_raw
        return g, recs, conv

    J = 1                       # state 2 as target
    X = Xf[:, J]
    g0, r0, c0 = chain(J, np.zeros(nij))
    g1, r1, c1 = chain(J, X)
    g2, r2, c2 = chain(J, 2.0 * X)
    print(f'z conv: {c0} {c1} {c2}')
    q = g2 - 4.0 * g1 + 3.0 * g0
    print(f'Q-TEST gradient: |g(2X)-4g(X)+3g(0)| max = {np.abs(q).max():.3e} '
          f'(|g(X)-g(0)| max = {np.abs(g1-g0).max():.3e})')
    for k in RECS:
        if r1[k] is None:
            print(f'  {k}: MISSING')
            continue
        a0 = r0[k]
        d1 = r1[k] - a0
        d2 = r2[k] - a0
        dev = np.abs(d2 - 4.0 * d1).max()
        print(f'  {k}: max|D(2X)-4D(X)| = {dev:.3e}  (|D(X)|max={np.abs(d1).max():.3e})')

    # cross-slot bilinearity with a fixed second vector: y = X of state 3
    Y = Xf[:, 2] if nstate > 2 else np.roll(X, 1)
    gy, ry, cy = chain(J, Y)
    gxy, rxy, cxy = chain(J, X + Y)
    print(f'z conv (y, x+y): {cy} {cxy}')
    cross = gxy - g1 - gy + g0
    # compare the same cross term from the scaled route:
    gxy2, _, _ = chain(J, X + 2.0 * Y)
    cross2 = 0.5 * (gxy2 - g1 - (chain(J, 2.0 * Y)[0]) + g0)
    print(f'CROSS consistency: |cross(X,Y) - cross2(X,Y)| max = '
          f'{np.abs(cross - cross2).max():.3e} (|cross|max={np.abs(cross).max():.3e})')


if __name__ == '__main__':
    main()
