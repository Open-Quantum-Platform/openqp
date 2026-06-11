"""Fit block-resolved gamma weights:  dn - d_amp = sum_b c_b * (G_b . theta).

Blocks of gamma (doc d, socc s, virt v): dd, ds, sd, ss, sv, vs, vv
(dv/vd are zero at TDA). Fit jointly over HF and BHHLYP datasets.
"""
import numpy as np

BLOCKS = ['dd', 'ds', 'sd', 'ss', 'sv', 'vs', 'vv']


def masks(noca, nocb, nbf):
    sp = np.zeros(nbf, dtype=int)
    sp[nocb:noca] = 1
    sp[noca:] = 2
    out = {}
    for bi, b0 in enumerate('dsv'):
        for bj, b1 in enumerate('dsv'):
            out[b0 + b1] = (sp[:, None] == bi) & (sp[None, :] == bj)
    return out


rows, targets, labels = [], [], []
per_set = {}
for tag in ('hf', 'bhh'):
    d = np.load(f'/tmp/nactest/theta_data_{tag}.npz')
    noca, nocb, nbf = int(d['noca']), int(d['nocb']), int(d['nbf'])
    gam = d['tdF'][2, 1].reshape((nbf, nbf)).T      # gamma^{23}
    theta = d['theta']
    t = d['dn_al'] - d['damp']                      # true orbital target
    mk = masks(noca, nocb, nbf)
    cols = []
    for b in BLOCKS:
        gb = np.where(mk[b], gam, 0.0)
        cols.append(np.einsum('pq,kpq->k', gb, theta))
    A = np.stack(cols, axis=1)                      # (9, nblocks)
    rows.append(A)
    targets.append(t)
    per_set[tag] = (A, t)
    print(f'{tag}: block contraction norms:',
          {b: round(np.linalg.norm(c), 4) for b, c in zip(BLOCKS, cols)})

A = np.vstack(rows)
t = np.concatenate(targets)
c, res, rank, sv = np.linalg.lstsq(A, t, rcond=None)
fit = A @ c
print('\njoint fit over both functionals:')
for b, cb in zip(BLOCKS, c):
    print(f'  c_{b} = {cb:+.5f}')
print('fit resid =', np.linalg.norm(fit - t), ' |t| =', np.linalg.norm(t))
for tag in ('hf', 'bhh'):
    Ai, ti = per_set[tag]
    fi = Ai @ c
    cos = np.dot(fi, ti) / (np.linalg.norm(fi) * np.linalg.norm(ti) + 1e-30)
    print(f'{tag}: cos={cos:+.6f} |fit|={np.linalg.norm(fi):.5f} '
          f'|t|={np.linalg.norm(ti):.5f} resid={np.linalg.norm(fi-ti):.5f}')

# also per-dataset independent fits, to see if c_b are transferable
print('\nper-dataset fits:')
for tag in ('hf', 'bhh'):
    Ai, ti = per_set[tag]
    ci, *_ = np.linalg.lstsq(Ai, ti, rcond=None)
    print(f'{tag}: ', {b: round(x, 4) for b, x in zip(BLOCKS, ci)},
          ' resid=%.5f' % np.linalg.norm(Ai @ ci - ti))
