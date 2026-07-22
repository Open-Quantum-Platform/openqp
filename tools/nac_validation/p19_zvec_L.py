"""
Step 3 decisive test: drive the z-vector with the BARE orbital gradient L as RHS
(via the new OQP::nac_orbgrad_L injection) and test whether the resulting
orbital-response contribution closes d_amp against the oracle.

Assembly under test:   oracle*gap  ?=  explicit_2e  +  a * d_L
  explicit_2e = mrsf_nac_amp (frozen-C explicit 2e).
  d_L = gZ - gS from the CPHF seam with RHS = L (bare orbital gradient, correct
        layout mo_a = VEC_MO_A.T).
If a single sign a = +-1 makes the residual small on all pairs, the bare-L RHS is
validated and the closed form is one explicit-Fock term away. We report a via a
per-pair least-squares fit and the residual after removing a*d_L.
"""
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint, NAC

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
EPS = 1.0e-4
bkey = 'OQP::td_bvec_mo'

r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p19.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
Craw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = Craw.shape[0]
mo0 = Craw.T.copy()
Crb = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
mo0b = Crb.T.copy()
nvirb = nbf - nocb
nij = noca * nvirb
X0_raw = np.array(mol.data[bkey], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]
pairs = [(i, j) for i in range(1, nstate + 1) for j in range(1, nstate + 1) if i < j]


def set_mo(moa, mob):
    mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(moa.T)
    mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(mob.T)


def set_bvec(col):
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = col
    mol.data[bkey] = rr.reshape(Xshape)


def Eij(I, J):
    set_bvec(X0[J])
    oqp.mrsf_matvec_apply(mol)
    return float(X0[I] @ np.array(mol.data['OQP::nac_mvax'], copy=True).ravel())


# ---- oracle + explicit 2e ----
oracle = NAC(mol)._compute_amp_damp(dx=1.0e-3)
set_mo(mo0, mo0b)
mol.data[bkey] = X0_raw
oqp.mrsf_nac_amp(mol)
ana2e = np.array(mol.data['OQP::nac_amp'], copy=True).reshape(-1).reshape((nstate, nstate, natom, 3))

# ---- bare orbital gradient L for each pair (0-based storage) ----
Lp = {(i, j): np.zeros((nbf, nbf)) for (i, j) in pairs}
for p in range(nbf):
    for q in range(nbf):
        if p == q:
            continue
        for sgn in (+1, -1):
            moa, mob = mo0.copy(), mo0b.copy()
            moa[:, q] += sgn * EPS * mo0[:, p]
            mob[:, q] += sgn * EPS * mo0b[:, p]
            set_mo(moa, mob)
            for (i, j) in pairs:
                Lp[(i, j)][p, q] += sgn * Eij(i - 1, j - 1) / (2.0 * EPS)
set_mo(mo0, mo0b)

# ---- CPHF seam driven by L (via OQP::nac_orbgrad_L) ----
def push_L(L):
    mol.data.remove_records(['OQP::nac_orbgrad_L']) if hasattr(mol.data, 'remove_records') else None
    mol.data['OQP::nac_orbgrad_L'] = np.ascontiguousarray(L.T).reshape(-1)


def clear_L():
    try:
        mol.data.remove_records(['OQP::nac_orbgrad_L'])
    except Exception:
        pass


d_L = {}
mol.data[bkey] = X0_raw
for (i, j) in pairs:
    push_L(Lp[(i, j)])
    mol.data.set_tdhf_target(i)
    oqp.set_mrsf_nac_cphf(mol, i, j)
    oqp.tdhf_mrsf_z_vector(mol)
    if not mol.mol_energy.Z_Vector_converged:
        oqp.set_mrsf_nac_cphf(mol, 0, 0)
        clear_L()
        d_L[(i, j)] = None
        continue
    oqp.tdhf_mrsf_gradient(mol)
    gZ = mol.get_grad().reshape((natom, 3)).copy()
    mol.data['OQP::td_p'] = np.zeros_like(np.array(mol.data['OQP::td_p'], copy=True))
    mol.data['OQP::WAO'] = np.zeros_like(np.array(mol.data['OQP::WAO'], copy=True))
    oqp.tdhf_mrsf_gradient(mol)
    gS = mol.get_grad().reshape((natom, 3)).copy()
    oqp.set_mrsf_nac_cphf(mol, 0, 0)
    clear_L()
    d_L[(i, j)] = gZ - gS

# ---- assembly test ----
print('# z-vector with bare-L RHS: oracle*gap ?= explicit_2e + a*d_L')
print(f'# {"pair":>7} {"|oracle*g|":>11} {"|exp2e|":>10} {"|d_L|":>10} '
      f'{"a(fit)":>8} {"|resid|":>10} {"resid/oracle":>12} {"cos(res)":>9}')
for (i, j) in pairs:
    gap = Om[j - 1] - Om[i - 1]
    orc = (oracle[(i, j)].reshape(-1)) * gap
    e2 = ana2e[i - 1, j - 1].reshape(-1)
    if d_L[(i, j)] is None:
        print(f'  {str((i,j)):>7}  z-vector not converged')
        continue
    dl = d_L[(i, j)].reshape(-1)
    tgt = orc - e2                     # what d_L must supply
    a = (tgt @ dl) / (dl @ dl + 1e-300)
    resid = tgt - a * dl
    c = tgt @ dl / (np.linalg.norm(tgt) * np.linalg.norm(dl) + 1e-300)
    print(f'  {str((i,j)):>7} {np.linalg.norm(orc):11.5f} {np.linalg.norm(e2):10.5f} '
          f'{np.linalg.norm(dl):10.5f} {a:8.3f} {np.linalg.norm(resid):10.5f} '
          f'{np.linalg.norm(resid)/(np.linalg.norm(orc)+1e-300):12.4f} {c:+9.4f}')

np.savez('/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p19_zvecL.npz',
         **{f'oracle_{i}{j}': oracle[(i, j)] for (i, j) in pairs},
         **{f'ana2e_{i}{j}': ana2e[i - 1, j - 1] for (i, j) in pairs},
         **{f'dL_{i}{j}': (d_L[(i, j)] if d_L[(i, j)] is not None else np.zeros((natom, 3))) for (i, j) in pairs},
         Om=np.array(Om))
print('\nsaved -> data_snapshots/p19_zvecL.npz')
