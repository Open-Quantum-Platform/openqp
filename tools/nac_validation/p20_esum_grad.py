"""
Step 3, final term: explicit_esum = Tr(P^IJ dF_AO/dR), the explicit Fock-derivative
gradient contracted with the interstate density P^IJ = C Gam C^T. Computed by
injecting P^IJ as OQP::td_p and running the existing gradient seam (gradient with
P^IJ minus SCF baseline), reusing the whole 1e/2e gradient engine.

Gam (esum trace identity, verified):
  Gam_A = -sym(Xt_I Xt_J^T) on alpha occ-occ,  Gam_B = +sym(Xt_I^T Xt_J) on beta
  virt-virt, Xt = SOMO-unfolded amplitude. P^IJ_a = mo_a Gam_A mo_a^T (AO).

Full assembly under test:
  oracle*gap  ?=  explicit_2e  +  f*explicit_esum  +  s*d_L
with explicit_2e = mrsf_nac_amp, d_L from p19 (z-vector, bare-L RHS). Reads the
saved p19 data for oracle / ana2e / d_L. If a constant (f, s) closes all pairs,
the closed-form d_amp is assembled.
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

r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p20.log')
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
    """lower-triangular pack, same convention as OQP::FOCK_A/td_p."""
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
    Pa = mo_a @ GA @ mo_a.T
    Pb = mo_b @ GB @ mo_b.T
    return Pa, Pb


# gradient seam: inject P^IJ as td_p, others zeroed, grad - baseline
tp = np.array(mol.data['OQP::td_p'], copy=True)
tp_shape = tp.shape
wao = np.array(mol.data['OQP::WAO'], copy=True)


def grad_with_tdp(pa_pack, pb_pack):
    newtp = np.zeros(tp_shape)
    # td_p is (npack, 2): column 0 alpha, column 1 beta
    if tp_shape[0] == 2:            # (2, npack)
        newtp[0, :] = pa_pack
        newtp[1, :] = pb_pack
    else:                            # (npack, 2)
        newtp[:, 0] = pa_pack
        newtp[:, 1] = pb_pack
    mol.data['OQP::td_p'] = newtp
    mol.data['OQP::WAO'] = np.zeros_like(wao)
    try:
        ab = np.array(mol.data['OQP::td_abxc'], copy=True)
        mol.data['OQP::td_abxc'] = np.zeros_like(ab)
    except Exception:
        pass
    oqp.tdhf_mrsf_gradient(mol)
    return mol.get_grad().reshape((natom, 3)).copy()


print(f'# td_p shape = {tp_shape}')
esum = {}
for (i, j) in pairs:
    mol.data.set_tdhf_target(i)
    Pa, Pb = Pij(i - 1, j - 1)
    gE = grad_with_tdp(pack(Pa), pack(Pb))
    gS = grad_with_tdp(np.zeros(nbf * (nbf + 1) // 2), np.zeros(nbf * (nbf + 1) // 2))
    esum[(i, j)] = gE - gS
mol.data['OQP::td_p'] = tp
mol.data['OQP::WAO'] = wao

print('# full assembly:  oracle*gap  ?=  explicit_2e + f*esum + s*d_L')
print(f'# {"pair":>7} {"|orc*g|":>9} {"|e2e|":>9} {"|esum|":>9} {"|d_L|":>9}'
      f' {"best f":>7} {"best s":>7} {"|resid|":>9} {"resid/o":>8}')
for (i, j) in pairs:
    gap = Om[j - 1] - Om[i - 1]
    orc = d[f'oracle_{i}{j}'].reshape(-1) * gap
    e2 = d[f'ana2e_{i}{j}'].reshape(-1)
    dl = d[f'dL_{i}{j}'].reshape(-1)
    es = esum[(i, j)].reshape(-1)
    # least-squares over (f, s): orc - e2 = f*es + s*dl
    tgt = orc - e2
    A = np.vstack([es, dl]).T
    coef, *_ = np.linalg.lstsq(A, tgt, rcond=None)
    f, s = coef
    resid = tgt - A @ coef
    print(f'  {str((i,j)):>7} {np.linalg.norm(orc):9.4f} {np.linalg.norm(e2):9.4f} '
          f'{np.linalg.norm(es):9.4f} {np.linalg.norm(dl):9.4f} {f:7.3f} {s:7.3f} '
          f'{np.linalg.norm(resid):9.5f} {np.linalg.norm(resid)/(np.linalg.norm(orc)+1e-30):8.4f}')
np.savez('/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p20_esum.npz',
         **{f'esum_{i}{j}': esum[(i, j)] for (i, j) in pairs})
print('\nsaved -> data_snapshots/p20_esum.npz')
