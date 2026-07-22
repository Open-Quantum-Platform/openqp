"""
The one test (from the derivation): closed-form d_amp assembly with a CLEAN
explicit_esum extracted by epsilon-linearisation.

explicit_esum = Tr(P^IJ dF_AO/dR) is LINEAR in P^IJ. The p20 injection used full
P (eps=1), so the gradient also captured the P*P quadratic term (~10x too big).
Here P is scaled by a small eps and the linear part isolated:
    explicit_esum = [ grad(td_p = eps*P) - grad(td_p = 0) ] / eps ,  eps -> 0.
The P*P term is O(eps^2) and vanishes; the P*rho(ref) cross + XC-kernel cross
(both = pieces of dF/dR) survive.

Full assembly test:  oracle*gap  ?=  explicit_2e + explicit_esum + s*d_L
with a SINGLE sign s (= the z-vector orbital response). If the residual is small
with one constant s on all pairs, the closed form is assembled.
"""
import math
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
P19 = '/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p19_zvecL.npz'
SQ = 1.0 / math.sqrt(2.0)
EPS = 1.0e-3          # linearisation strength for P

r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p21.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
Craw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = Craw.shape[0]
mo_a = Craw.T.copy()
Crb = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
mo_b = Crb.T.copy()
nvirb = nbf - nocb
nij = noca * nvirb
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
ijlr1 = (noca - 1 - nocb - 1) * noca + (noca - 1) - 1
ijlr2 = (noca - nocb - 1) * noca + (noca) - 1
pairs = [(i, j) for i in range(1, nstate + 1) for j in range(1, nstate + 1) if i < j]
d = np.load(P19)
Om = d['Om']


def unfold(col):
    x = np.zeros((noca, nvirb))
    for i in range(noca):
        for a in range(nvirb):
            ij = a * noca + i
            if ij == ijlr1:
                x[i, a] = col[ijlr1] * SQ
            elif ij == ijlr2:
                x[i, a] = -col[ijlr1] * SQ
            else:
                x[i, a] = col[ij]
    return x


Xt = [unfold(X0[k]) for k in range(nstate)]


def pack(M):
    out = np.zeros(nbf * (nbf + 1) // 2)
    idx = 0
    for i in range(nbf):
        for j in range(i + 1):
            out[idx] = M[i, j]
            idx += 1
    return out


def Pij(I, J):
    GA = np.zeros((nbf, nbf))
    GB = np.zeros((nbf, nbf))
    GA[0:noca, 0:noca] = -0.5 * (Xt[I] @ Xt[J].T + Xt[J] @ Xt[I].T)
    GB[nocb:nbf, nocb:nbf] = 0.5 * (Xt[I].T @ Xt[J] + Xt[J].T @ Xt[I])
    return mo_a @ GA @ mo_a.T, mo_b @ GB @ mo_b.T


tp = np.array(mol.data['OQP::td_p'], copy=True)
tp_shape = tp.shape
wao = np.array(mol.data['OQP::WAO'], copy=True)


def grad_tdp(pa, pb):
    newtp = np.zeros(tp_shape)
    if tp_shape[0] == 2:
        newtp[0, :] = pa
        newtp[1, :] = pb
    else:
        newtp[:, 0] = pa
        newtp[:, 1] = pb
    mol.data['OQP::td_p'] = newtp
    mol.data['OQP::WAO'] = np.zeros_like(wao)
    try:
        mol.data['OQP::td_abxc'] = np.zeros_like(np.array(mol.data['OQP::td_abxc'], copy=True))
    except Exception:
        pass
    oqp.tdhf_mrsf_gradient(mol)
    return mol.get_grad().reshape((natom, 3)).copy()


z = np.zeros(nbf * (nbf + 1) // 2)
esum = {}
for (i, j) in pairs:
    mol.data.set_tdhf_target(i)
    Pa, Pb = Pij(i - 1, j - 1)
    gP = grad_tdp(EPS * pack(Pa), EPS * pack(Pb))
    g0 = grad_tdp(z, z)
    esum[(i, j)] = (gP - g0) / EPS      # linear part = Tr(P dF/dR)
mol.data['OQP::td_p'] = tp
mol.data['OQP::WAO'] = wao

print(f'# eps-linearised explicit_esum; full assembly '
      f'oracle*gap ?= e2e + esum + s*d_L')
print(f'# {"pair":>7} {"|orc*g|":>9} {"|e2e|":>9} {"|esum|":>9} {"|d_L|":>9}'
      f' {"s(fit)":>8} {"|resid|":>9} {"resid/o":>8} {"cos(res,dL)":>12}')
for (i, j) in pairs:
    gap = Om[j - 1] - Om[i - 1]
    orc = d[f'oracle_{i}{j}'].reshape(-1) * gap
    e2 = d[f'ana2e_{i}{j}'].reshape(-1)
    dl = d[f'dL_{i}{j}'].reshape(-1)
    es = esum[(i, j)].reshape(-1)
    tgt = orc - e2 - es                 # what d_L must supply
    s = (tgt @ dl) / (dl @ dl + 1e-300)
    resid = tgt - s * dl
    c = tgt @ dl / (np.linalg.norm(tgt) * np.linalg.norm(dl) + 1e-300)
    print(f'  {str((i,j)):>7} {np.linalg.norm(orc):9.4f} {np.linalg.norm(e2):9.4f} '
          f'{np.linalg.norm(es):9.4f} {np.linalg.norm(dl):9.4f} {s:8.3f} '
          f'{np.linalg.norm(resid):9.5f} {np.linalg.norm(resid)/(np.linalg.norm(orc)+1e-30):8.4f} {c:+12.4f}')
