"""Channel-resolved amplitude term: reweight G-slot cross-bilinears by -2."""
import math
import oqp
from oqp.pyoqp import Runner
import numpy as np

inp = '/tmp/nactest/H2O_tight_dx0.001.inp'
r = Runner(input_file=inp, log='/tmp/nactest/chfix.log')
r.run()
mol = r.mol
nstate = 3
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca - 2
nbf = np.array(mol.data['OQP::VEC_MO_A']).shape[0]
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
nij = X0_raw.size // nstate
X = X0_raw.reshape(-1).reshape((nstate, nij))

# G-channel slots (1-based): ijg = (O2->O1), ijd = (O1->O2)
ijg = (noca-1-nocb-1)*noca + noca        # = 6 for noca=6
ijd = (noca-nocb-1)*noca + noca - 1      # = 11
PG = np.zeros(nij); PG[ijg-1] = 1.0; PG[ijd-1] = 1.0

def Gfull(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(X0_raw.shape)
    mol.data.set_tdhf_target(1)
    oqp.tdhf_mrsf_z_vector(mol)
    oqp.tdhf_mrsf_gradient(mol)
    g = mol.get_grad().reshape(-1).copy()
    mol.data['OQP::td_bvec_mo'] = X0_raw
    return g

def B(u, v):                              # bilinear via polarization (exact)
    return 0.25*(Gfull(u+v) - Gfull(u-v))

d = np.load('/tmp/nactest/gt2.npz')
damp_PT = d['damp_PT']
print('\n===== channel-resolved amplitude term vs damp_PT =====')
for (i, j) in ((1, 2), (1, 3), (2, 3)):
    XI, XJ = X[i-1], X[j-1]
    gI, bI = XI*PG, XI*(1-PG)
    gJ, bJ = XJ*PG, XJ*(1-PG)
    Bgg = B(gI, gJ); Bgb = B(gI, bJ); Bbg = B(bI, gJ); Bbb = B(bI, bJ)
    gap = mol.energies[j] - mol.energies[i]
    for (wgg, wx, name) in ((1.0, -2.0, 'w_GG=+1 w_Gx=-2'),
                            (-2.0, -2.0, 'w_GG=-2 w_Gx=-2'),
                            (1.0, 1.0,  'no fix          ')):
        amp = (wgg*Bgg + wx*(Bgb + Bbg) + Bbb) / gap
        t = damp_PT[:, i-1, j-1]
        c = np.dot(amp, t)/(np.linalg.norm(amp)*np.linalg.norm(t)+1e-30)
        print(f'({i},{j}) {name}: cos={c:+.6f} ratio={np.linalg.norm(amp)/np.linalg.norm(t):.4f} '
              f'resid={np.linalg.norm(amp-t):.5f}  (|t|={np.linalg.norm(t):.5f})')
    print()
