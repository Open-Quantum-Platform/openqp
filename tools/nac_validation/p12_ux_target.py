"""
Phase 12 step 0: regenerate the U^x TARGET for the closed-form amplitude term.

The amplitude term of the analytical MRSF NAC is
    d_amp(I,J) = X_I^T (dA/dR) X_J / (Om_J - Om_I)
and A depends on R both explicitly (through the AO integrals) and implicitly
(through the MO coefficients).  The established decomposition is

    oracle  =  ana2e_explicit  +  2e_orbital_response  +  esum_full

`ana2e_explicit` is already closed-form in Fortran (mrsf_nac_amp, tag
OQP::nac_amp).  Everything else is dominated by the INTERSTATE ORBITAL RESPONSE
U^x, which is the one genuinely missing piece.  Its target is therefore

    missing(I,J) = oracle(I,J) - ana2e(I,J)/gap

where `oracle` is the semi-numerical transported-matvec damp that the production
analytical NAC currently uses (and which is validated against both the numerical
NAC and a GAMESS-tracked reference).

This regenerates that target (the previous copy lived in /tmp and was swept by
the system tmp cleaner) and saves it next to the harnesses, NOT in /tmp.

Run under the OpenQP NAC fork env; see tools/nac_validation/README or
~/gamess/h2nac/run_oqp.sh.
"""
import os
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import NAC

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
OUTNPZ = sys.argv[2] if len(sys.argv) > 2 else \
    '/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p12_ux_target.npz'

r = Runner(input_file=INP, log=os.path.splitext(OUTNPZ)[0] + '.log')
r.run()
mol = r.mol

nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]

nac = NAC(mol)

# ---------------------------------------------------------------- oracle
# semi-numerical, but it IS the validated production amplitude term
damp = nac._compute_amp_damp(dx=1.0e-3)
if damp is None:
    raise SystemExit('FAIL: _compute_amp_damp returned None')

# ------------------------------------------------- closed-form 2e explicit
# mrsf_nac_amp fills OQP::nac_amp with the RAW G_IJ = X_I^T (d_x A_2e) X_J
# for every pair, shape (3, natom, nstate, nstate) in Fortran order.
oqp.mrsf_nac_amp(mol)
raw = np.array(mol.data['OQP::nac_amp'], copy=True)
ana2e = raw.reshape(-1).reshape((nstate, nstate, natom, 3))

out = {'Om': np.array(Om), 'natom': natom, 'nstate': nstate}
print(f'# H2O-like target, nstate={nstate}, natom={natom}')
print(f'# {"pair":>7} {"|oracle|":>11} {"|ana2e/gap|":>12} {"|missing|":>11} '
      f'{"cos(ana2e,oracle)":>18} {"|miss|/|oracle|":>15}')
for I in range(nstate):
    for J in range(nstate):
        if I == J:
            continue
        key = (I + 1, J + 1)
        if damp.get(key) is None:
            continue
        gap = Om[J] - Om[I]
        orc = damp[key].reshape(-1)
        a2e = (ana2e[I, J] / gap).reshape(-1)
        miss = orc - a2e
        out[f'oracle_{I+1}{J+1}'] = orc
        out[f'ana2e_{I+1}{J+1}'] = a2e
        out[f'missing_{I+1}{J+1}'] = miss
        c = (a2e @ orc) / (np.linalg.norm(a2e) * np.linalg.norm(orc) + 1e-300)
        print(f'  {str(key):>7} {np.linalg.norm(orc):11.6f} {np.linalg.norm(a2e):12.6f} '
              f'{np.linalg.norm(miss):11.6f} {c:+18.6f} '
              f'{np.linalg.norm(miss)/(np.linalg.norm(orc)+1e-300):15.3f}')

os.makedirs(os.path.dirname(OUTNPZ), exist_ok=True)
np.savez(OUTNPZ, **out)
print(f'\nsaved -> {OUTNPZ}')
