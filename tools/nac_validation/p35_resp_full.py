"""
resp, tested against the VALIDATED esum -- everything in ONE run.

  X_I^T (dA/dR) X_J = ana2e + esum + resp
  resp_target = oracle*gap - ana2e - esum        (esum now FD-validated: 1e/2e/XC)

Structure under test:  resp = L : U^x
  L_pq  = d(X_I^T A X_J)/d theta_pq   -- bare orbital gradient of the MATVEC
          coupling (orbital-rotation FD of oqp.mrsf_matvec_apply). This is the
          CORRECT operator (see the homogeneity proof: the gradient chain's
          operator differs off-diagonal, the matvec's does not).
  U^x   = dM/dR, M = C_ref^T S C_disp -- orbital response (nuclear FD, re-SCF).
          Split into antisymmetric (true rotation) and symmetric (orthonormality,
          = -1/2 S^x_MO) parts; BOTH are tested, since the symmetric part carries
          the -Tr[W S^x]-like piece.

ORDERING (critical): ana2e/esum/L need UNPERTURBED orbitals and must be computed
BEFORE anything that re-SCFs; U^x and the oracle perturb geometry/MOs, so they go
last. MRSF = ROHF + functional -> production BHHLYP input only.
"""
import sys
import numpy as np
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint, NAC

INP = sys.argv[1] if len(sys.argv) > 1 else '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
EPS = 1.0e-4
DX = 1.0e-3
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p35.log')
r.run()
mol = r.mol
nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
Craw = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
nbf = Craw.shape[0]
mo0 = Craw.T.copy()
Crb = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
mo0b = Crb.T.copy()
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
nvirb = nbf - nocb
nij = noca * nvirb
bkey = 'OQP::td_bvec_mo'
X0_raw = np.array(mol.data[bkey], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
xyz0 = np.array(mol.get_system(), copy=True)
E = list(mol.energies)
Om = [E[k + 1] - E[0] for k in range(nstate)]
nc = 3 * natom
pairs = [(i, j) for i in range(1, nstate + 1) for j in range(1, nstate + 1) if i < j]


def set_mo(moa, mob):
    mol.data['OQP::VEC_MO_A'] = np.ascontiguousarray(moa.T)
    mol.data['OQP::VEC_MO_B'] = np.ascontiguousarray(mob.T)


def set_bvec(col):
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = col
    mol.data[bkey] = rr.reshape(Xshape)


def Eij(I, J):
    """X_I^T A X_J via the MATVEC (0-based state indices)."""
    set_bvec(X0[J])
    oqp.mrsf_matvec_apply(mol)
    return float(X0[I] @ np.array(mol.data['OQP::nac_mvax'], copy=True).ravel())


def cos(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-300))


# ---------- (1) clean-orbital quantities ----------
mol.data[bkey] = X0_raw
oqp.mrsf_nac_amp(mol)
ana2e = np.array(mol.data['OQP::nac_amp'], copy=True).reshape(-1).reshape(
    (nstate, nstate, natom, 3))
esum, wsx = {}, {}
for (i, j) in pairs:
    oqp.mrsf_nac_esum(mol, i, j)
    esum[(i, j)] = np.array(mol.data['OQP::nac_esum'], copy=True).reshape(-1)
    wsx[(i, j)] = np.array(mol.data['OQP::nac_wsx'], copy=True).reshape(-1)

# ---------- (2) L = bare orbital gradient of the matvec coupling ----------
L = {p: np.zeros((nbf, nbf)) for p in pairs}
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
                L[(i, j)][p, q] += sgn * Eij(i - 1, j - 1) / (2.0 * EPS)
set_mo(mo0, mo0b)
mol.data[bkey] = X0_raw

# ---------- (3) U^x = dM/dR (perturbs: re-SCF) ----------
cfg = mol.config
json0 = mol.log.replace('.log', '.json')
mol.save_data()
gb = dict(cfg['guess'])
cfg['guess']['type'] = 'json'
cfg['guess']['file'] = json0
cfg['guess']['continue_geom'] = False


def mo_overlap(coord):
    mol.update_system(coord)
    oqp.library.ints_1e(mol)
    oqp.library.guess(mol)
    SinglePoint(mol).energy()
    mol.data['OQP::xyz_old'] = xyz0.reshape((3, natom))
    mol.data['OQP::VEC_MO_A_old'] = Craw
    mol.data['OQP::E_MO_A_old'] = e0
    mol.data['OQP::VEC_MO_B_old'] = Crb
    mol.data['OQP::E_MO_B_old'] = e0b
    oqp.get_structures_ao_overlap(mol)
    return np.array(mol.data[OVTAG], copy=True).reshape(-1).reshape((nbf, nbf)).T


Ux = np.zeros((nc, nbf, nbf))
for k in range(nc):
    Ux[k] = (mo_overlap(xyz0 + DX * np.eye(nc)[k]) -
             mo_overlap(xyz0 - DX * np.eye(nc)[k])) / (2.0 * DX)
cfg['guess'].update(gb)
mol.update_system(xyz0)
oqp.library.ints_1e(mol)
oqp.library.guess(mol)
SinglePoint(mol).energy()
mol.data[bkey] = X0_raw

# ---------- (4) oracle (last) ----------
oracle = NAC(mol)._compute_amp_damp(dx=1.0e-3)
mol.data[bkey] = X0_raw

# ---------- (5) compare ----------
Ua = 0.5 * (Ux - np.transpose(Ux, (0, 2, 1)))
Us = 0.5 * (Ux + np.transpose(Ux, (0, 2, 1)))
print('# resp vs L:U^x, against the VALIDATED-esum target (one run, BHHLYP/ROHF)')
print(f'{"pair":>6} {"|resp_t|":>9} {"proj":>9} {"cos":>9} {"ratio":>8}')
res = {}
for (i, j) in pairs:
    gap = Om[j - 1] - Om[i - 1]
    orc = oracle[(i, j)].reshape(-1) * gap
    a2 = ana2e[i - 1, j - 1].reshape(-1)
    rt = orc - a2 - esum[(i, j)]
    res[(i, j)] = rt
    for tag, U in [('Ux_full', Ux), ('Ua_anti', Ua), ('Us_sym', Us)]:
        v = np.array([np.sum(L[(i, j)] * U[k]) for k in range(nc)])
        print(f'{str((i,j)):>6} {np.linalg.norm(rt):9.5f} {tag:>9} '
              f'{cos(v, rt):+9.4f} {np.linalg.norm(v)/(np.linalg.norm(rt)+1e-30):8.4f}')
    w = wsx[(i, j)]
    print(f'{"":>6} {"":>9} {"wsx":>9} {cos(w, rt):+9.4f} '
          f'{np.linalg.norm(w)/(np.linalg.norm(rt)+1e-30):8.4f}')
    print()

np.savez('/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p35_resp.npz',
         Ux=Ux, **{f'L_{i}{j}': L[(i, j)] for (i, j) in pairs},
         **{f'resp_{i}{j}': res[(i, j)] for (i, j) in pairs},
         **{f'esum_{i}{j}': esum[(i, j)] for (i, j) in pairs},
         **{f'wsx_{i}{j}': wsx[(i, j)] for (i, j) in pairs},
         **{f'ana2e_{i}{j}': ana2e[i - 1, j - 1] for (i, j) in pairs},
         **{f'oracle_{i}{j}': oracle[(i, j)] for (i, j) in pairs},
         Om=np.array(Om), noca=noca, nocb=nocb, nbf=nbf, natom=natom)
print('saved -> data_snapshots/p35_resp.npz')
