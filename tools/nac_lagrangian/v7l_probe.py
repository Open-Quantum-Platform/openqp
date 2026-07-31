"""v7l: does the MRSF sigma depend on the DM records DIRECTLY (Fock
frozen)? Perturb DM_A/DM_B, keep FOCK_A/B, re-apply the matvec.
Run: python v7l_probe.py <input.inp>
"""
import os, sys
import numpy as np

def main():
    import oqp
    from oqp.pyoqp import Runner
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nac_formula_kernel as FK
    inp = sys.argv[1]
    r = Runner(input_file=inp, log=inp.replace('.inp', '_v7l.log'))
    r.run()
    mol = r.mol
    ctx = FK.build_context(mol)
    nij = ctx['nij']
    X0_raw = ctx['X0_raw']
    nstate = mol.config['tdhf']['nstate']
    try:
        mol.data._data.control.int2e_cutoff = 1e-20
    except Exception:
        pass
    Xf = X0_raw.reshape(-1).reshape((nstate, nij)).T.copy()

    def matvec(v):
        rr = X0_raw.copy().reshape(-1)
        rr[0:nij] = v
        mol.data['OQP::td_bvec_mo'] = rr.reshape(X0_raw.shape)
        oqp.mrsf_matvec_apply(mol)
        return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()[:nij].copy()

    d0a = np.array(mol.data['OQP::DM_A'], copy=True)
    d0b = np.array(mol.data['OQP::DM_B'], copy=True)
    base = matvec(Xf[:, 1])
    eps = 1e-4
    rng = np.sin(np.arange(d0a.size) * 1.7) * np.abs(d0a).max()
    mol.data['OQP::DM_A'] = (d0a.ravel() + eps * rng).reshape(d0a.shape)
    mol.data['OQP::DM_B'] = (d0b.ravel() + eps * rng[:d0b.size]).reshape(d0b.shape)
    pert = matvec(Xf[:, 1])
    mol.data['OQP::DM_A'] = d0a
    mol.data['OQP::DM_B'] = d0b
    print(f'|dAx|/eps = {np.abs(pert - base).max() / eps:.6e}  '
          f'(|Ax|max={np.abs(base).max():.4f})')
    # also SM probe
    s0 = np.array(mol.data['OQP::SM'], copy=True)
    mol.data['OQP::SM'] = (s0.ravel() + eps * rng[:s0.size] * 0.1).reshape(s0.shape)
    pert2 = matvec(Xf[:, 1])
    mol.data['OQP::SM'] = s0
    print(f'SM probe: |dAx|/eps = {np.abs(pert2 - base).max() / eps:.6e}')

if __name__ == '__main__':
    main()
