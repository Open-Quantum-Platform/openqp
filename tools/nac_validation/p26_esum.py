"""
Validate the CORRECTED closed-form interstate esum (mrsf_nac_esum).

esum = Tr(P^IJ_a dFOCK_A^skel/dR) + Tr(P^IJ_b dFOCK_B^skel/dR), built from the
PURE skeleton Fock derivative at frozen reference density (1e kinetic+nuclear via
grd1, 2e via fock_deriv_contract_os), deliberately EXCLUDING the W.dS^x overlap
term and P^IJ's own 2e response -- the two contaminations that made the p20 esum
~20x too large.

Checks, per pair:
  1. magnitude sanity vs the old p20 esum (expect MUCH smaller, ~3.5x for (1,3));
  2. the residual that the orbital response must still supply:
       missing = oracle*gap - ana2e         (before esum)
       resid   = missing - esum             (after esum)
     If esum is right, |resid| should drop well below |missing|.
Run on both the BHHLYP input (comparable to the p20 snapshot; note esum omits the
XC-potential derivative there) and the HF input (no XC -> esum is COMPLETE).
"""
import os
import sys
import numpy as np
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import NAC

INP = sys.argv[1] if len(sys.argv) > 1 else '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
TAG = os.path.basename(INP).replace('.inp', '')
r = Runner(input_file=INP, log=f'/bighome/alireza/openqp-nac/tools/nac_validation/p26_{TAG}.log')
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

# IMPORTANT ORDERING: _compute_amp_damp does 6N re-SCF at DISPLACED geometries and
# does NOT restore the MOs/geometry (p18 has to call set_mo(mo0,mo0b) after it).
# So every quantity that depends on the unperturbed orbitals MUST be computed
# BEFORE the oracle, otherwise it is silently evaluated on perturbed orbitals
# (this produced ~30% run-to-run drift in esum before the reorder).

# ---- explicit 2e (clean orbitals) ----
oqp.mrsf_nac_amp(mol)
ana2e = np.array(mol.data['OQP::nac_amp'], copy=True).reshape(-1).reshape(
    (nstate, nstate, natom, 3))

# ---- NEW closed-form esum (clean orbitals) ----
esum = {}
esum1e = {}
for (i, j) in pairs:
    oqp.mrsf_nac_esum(mol, i, j)
    raw = np.array(mol.data['OQP::nac_esum'], copy=True)
    esum[(i, j)] = raw.reshape(-1).reshape((natom, 3))
    esum1e[(i, j)] = np.array(mol.data['OQP::nac_esum_1e'], copy=True).reshape(-1)

# ---- oracle LAST (it perturbs geometry/MOs) ----
mol.data[bkey] = raw0
oracle = NAC(mol)._compute_amp_damp(dx=1.0e-3)
mol.data[bkey] = raw0

# ---- old p20 esum for comparison (BHHLYP snapshot only) ----
p20 = None
try:
    p20 = np.load('/bighome/alireza/openqp-nac/tools/nac_validation/'
                  'data_snapshots/p20_esum.npz')
except Exception:
    pass

print(f'\n### {TAG} ###')
print(f'{"pair":>6} {"|orc*gap|":>10} {"|ana2e|":>9} {"|esum_new|":>10} '
      f'{"|esum_p20|":>10} {"|missing|":>9} {"|resid|":>9} {"resid/miss":>10}')
for (i, j) in pairs:
    gap = Om[j - 1] - Om[i - 1]
    orc = oracle[(i, j)].reshape(-1) * gap
    a2 = ana2e[i - 1, j - 1].reshape(-1)
    es = esum[(i, j)].reshape(-1)
    missing = orc - a2
    resid = missing - es
    e20 = (np.linalg.norm(p20[f'esum_{i}{j}'])
           if (p20 is not None and f'esum_{i}{j}' in p20.files) else float('nan'))
    print(f'{str((i,j)):>6} {np.linalg.norm(orc):10.5f} {np.linalg.norm(a2):9.5f} '
          f'{np.linalg.norm(es):10.5f} {e20:10.5f} {np.linalg.norm(missing):9.5f} '
          f'{np.linalg.norm(resid):9.5f} {np.linalg.norm(resid)/(np.linalg.norm(missing)+1e-30):10.4f}')
print('# resid/miss < 1 => esum removes real physics from the residual.')
print('# esum_new should be FAR below esum_p20 (p20 was ~20x too big).')
np.savez(f'/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p26_esum_{TAG}.npz',
         **{f'esum_{i}{j}': esum[(i, j)] for (i, j) in pairs},
         **{f'oracle_{i}{j}': oracle[(i, j)] for (i, j) in pairs},
         **{f'ana2e_{i}{j}': ana2e[i - 1, j - 1] for (i, j) in pairs},
         Om=np.array(Om))
print(f'saved -> data_snapshots/p26_esum_{TAG}.npz')

# ---- 1e/(2e+xc) split diagnostic (from the CLEAN pre-oracle pass) ----
print('\n# 1e/(2e+xc) split of esum')
print(f'{"pair":>6} {"|esum|":>10} {"|1e|":>10} {"|2e+xc|":>10}')
for (i, j) in pairs:
    tot = esum[(i, j)].reshape(-1)
    e1 = esum1e[(i, j)]
    print(f'{str((i,j)):>6} {np.linalg.norm(tot):10.5f} {np.linalg.norm(e1):10.5f} '
          f'{np.linalg.norm(tot-e1):10.5f}')
