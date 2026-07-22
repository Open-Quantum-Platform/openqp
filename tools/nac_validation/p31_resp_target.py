"""
Clean `resp` target, now that esum is FD-validated.

  X_I^T (dA/dR) X_J  =  ana2e  +  esum  +  resp
  =>  resp_target = oracle*gap - ana2e - esum

Every earlier resp diagnostic compared L:U^x against
  missing = oracle*gap - ana2e  =  esum + resp,
which is a NEARLY-CANCELLED residue (esum ~ -resp), so it was the wrong target.
This builds the right one and characterises it:
  - |resp_target| (expect ~|esum|, i.e. ~2-3, opposing esum)
  - cos(resp_target, esum)  (expect ~ -1 if they nearly cancel)
  - cos(resp_target, oracle*gap)
MRSF = ROHF + functional: run on the production BHHLYP input only.
ORDERING: ana2e/esum MUST be computed BEFORE the oracle (it re-SCFs at displaced
geometries and does not restore MOs/geometry).
"""
import sys
import numpy as np
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import NAC

INP = sys.argv[1] if len(sys.argv) > 1 else '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p31.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]
pairs = [(i, j) for i in range(1, nstate + 1) for j in range(1, nstate + 1) if i < j]
bkey = 'OQP::td_bvec_mo'
raw0 = np.array(mol.data[bkey], copy=True)

# ---- clean-orbital quantities FIRST ----
oqp.mrsf_nac_amp(mol)
ana2e = np.array(mol.data['OQP::nac_amp'], copy=True).reshape(-1).reshape(
    (nstate, nstate, natom, 3))
esum = {}
wsx = {}
for (i, j) in pairs:
    oqp.mrsf_nac_esum(mol, i, j)
    esum[(i, j)] = np.array(mol.data['OQP::nac_esum'], copy=True).reshape(-1)
    wsx[(i, j)] = np.array(mol.data['OQP::nac_wsx'], copy=True).reshape(-1)

# ---- oracle LAST (perturbs geometry/MOs) ----
mol.data[bkey] = raw0
oracle = NAC(mol)._compute_amp_damp(dx=1.0e-3)
mol.data[bkey] = raw0


def cos(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-300))


print('# CLEAN resp target = oracle*gap - ana2e - esum   (BHHLYP/ROHF)')
print(f'{"pair":>6} {"|orc*g|":>9} {"|ana2e|":>9} {"|esum|":>9} {"|resp_t|":>9} '
      f'{"cos(rt,esum)":>13} {"cos(rt,orc)":>12} {"|esum+rt|":>10}')
out = {}
for (i, j) in pairs:
    gap = Om[j - 1] - Om[i - 1]
    orc = oracle[(i, j)].reshape(-1) * gap
    a2 = ana2e[i - 1, j - 1].reshape(-1)
    es = esum[(i, j)]
    rt = orc - a2 - es
    out[(i, j)] = rt
    print(f'{str((i,j)):>6} {np.linalg.norm(orc):9.5f} {np.linalg.norm(a2):9.5f} '
          f'{np.linalg.norm(es):9.5f} {np.linalg.norm(rt):9.5f} '
          f'{cos(rt, es):13.5f} {cos(rt, orc):12.5f} {np.linalg.norm(es+rt):10.5f}')
print()
print('# Does the Fock-weighted overlap term -Tr[W.S^x] account for resp?')
print(f'{"pair":>6} {"|resp_t|":>9} {"|wsx|":>9} {"cos(wsx,rt)":>12} {"ratio":>8} {"|rt-wsx|":>9}')
for (i, j) in pairs:
    rt = out[(i, j)]; w = wsx[(i, j)]
    print(f'{str((i,j)):>6} {np.linalg.norm(rt):9.5f} {np.linalg.norm(w):9.5f} '
          f'{cos(w, rt):12.5f} {np.linalg.norm(w)/(np.linalg.norm(rt)+1e-30):8.4f} '
          f'{np.linalg.norm(rt-w):9.5f}')
print('# cos~+1 and ratio~1 => -Tr[W.S^x] IS the dominant resp piece.')

np.savez('/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p31_resp_target.npz',
         **{f'resp_{i}{j}': out[(i, j)] for (i, j) in pairs},
         **{f'esum_{i}{j}': esum[(i, j)] for (i, j) in pairs},
         **{f'ana2e_{i}{j}': ana2e[i - 1, j - 1] for (i, j) in pairs},
         **{f'oracle_{i}{j}': oracle[(i, j)] for (i, j) in pairs},
         Om=np.array(Om))
print('saved -> data_snapshots/p31_resp_target.npz')
