"""
Phase 12 step 3: isolate ONE defect at a time by testing the MOVING-FRAME,
FROZEN-FOCK chain-rule identity

    dE/dR  =  [dE/dR]_explicit  +  sum_pq R_pq U^x_pq                    (*)

with  E(C,R) = X_I^T A(C,R) X_J  evaluated at FIXED amplitude COMPONENTS and
with the AO Fock held frozen at its reference value in every evaluation.

WHY THIS AND NOT THE ORACLE.  The production oracle differentiates the
TRANSPORTED operator A_ref = T^T A T (T = amplitude-space transport built from
the MO overlap), which keeps the PHYSICAL state fixed.  Using A X_J = Om_J X_J
and antisymmetry of tau = dT/dR,

    X_I^T (dA_ref/dR) X_J = X_I^T (dA/dR) X_J + (Om_I - Om_J) X_I.tau.X_J

so the transported and moving-frame derivatives differ by a term that is LARGE
(it is exactly the amplitude gauge term).  Comparing sum R.U^x directly against
`missing` therefore mixes three separate defects at once - which is what the
first attempt (p12_orbgrad.py) did, giving cos -0.78/+0.06/+0.28.

Here every quantity is measured in ONE consistent world (moving frame, frozen
AO Fock), so (*) must hold to FD accuracy if and only if my R and U^x are
correctly defined.  This tests ONLY the R/U^x bookkeeping, with the transport
term and the density-response channel deliberately excluded from both sides.

PASS  -> R and U^x are right; remaining work is (a) the density-response channel
         in R, (b) the explicit esum part, (c) the transport term - each of which
         can then be added and tested separately.
FAIL  -> R or U^x is mis-defined; fix that before anything else. Write no Fortran.
"""
import sys

import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint

import numpy as np

INP = sys.argv[1] if len(sys.argv) > 1 else \
    '/bighome/alireza/gamess/h2nac/ana/h2o_ana.inp'
OUT = '/bighome/alireza/openqp-nac/tools/nac_validation/data_snapshots/p12_identity3.npz'
EPS = 1.0e-4     # MO-mixing step for R
DX = 1.0e-3      # nuclear step
OVTAG = 'OQP::overlap_mo_non_orthogonal'

r = Runner(input_file=INP,
           log='/bighome/alireza/openqp-nac/tools/nac_validation/p12_identity3.log')
r.run()
mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20

nstate = mol.config['tdhf']['nstate']
natom = mol.data['natom']
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0])
nocb = noca - 2
C0 = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
F0a = np.array(mol.data['OQP::FOCK_A'], copy=True)
F0b = np.array(mol.data['OQP::FOCK_B'], copy=True)
nbf = C0.shape[0]
nvirb = nbf - nocb
nij = noca * nvirb
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
xyz0 = np.array(mol.get_system(), copy=True)
ncoord = 3 * natom
pairs = [(I, J) for I in range(nstate) for J in range(nstate) if I < J]


def freeze_fock():
    mol.data['OQP::FOCK_A'] = F0a
    mol.data['OQP::FOCK_B'] = F0b


def set_bvec(col):
    rr = X0_raw.copy().reshape(-1)
    rr[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = rr.reshape(Xshape)


def E_all():
    """{(I,J): X_I^T A X_J} with the CURRENT MOs/geometry and the frozen AO Fock."""
    out = {}
    cache = {}
    for (I, J) in pairs:
        if J not in cache:
            set_bvec(X0[J])
            oqp.mrsf_matvec_apply(mol)
            cache[J] = np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()
        out[(I, J)] = float(X0[I] @ cache[J])
    return out


def put_C(Ca, Cb):
    mol.data['OQP::VEC_MO_A'] = Ca
    mol.data['OQP::VEC_MO_B'] = Cb


# ============================================================ 1. R_pq (frozen Fock)
freeze_fock()
R = {p: np.zeros((nbf, nbf)) for p in pairs}
for p in range(nbf):
    for q in range(nbf):
        # NB: the diagonal p==q is INCLUDED -- V_pp != 1 in general (column
        # rescaling) and dE/dV_pp != 0.  Skipping it was the second bug.
        for sgn in (+1, -1):
            Ca, Cb = C0.copy(), C0b.copy()
            Ca[:, q] += sgn * EPS * C0[:, p]
            Cb[:, q] += sgn * EPS * C0b[:, p]
            put_C(Ca, Cb)
            freeze_fock()
            e = E_all()
            for pr in pairs:
                R[pr][p, q] += sgn * e[pr] / (2.0 * EPS)
put_C(C0, C0b)
freeze_fock()

# ================================================== 2. displaced quantities
cfg = mol.config
json0 = mol.log.replace('.log', '.json')
mol.save_data()
guess_bak = dict(cfg['guess'])
cfg['guess']['type'] = 'json'
cfg['guess']['file'] = json0
cfg['guess']['continue_geom'] = False


def at_geom(coord, relax):
    """Move nuclei. relax=True -> re-SCF (C relaxes); False -> keep reference C.
    In BOTH cases the AO Fock is forced back to its reference value, so the whole
    identity is tested in one consistent frozen-Fock world."""
    mol.update_system(coord)
    oqp.library.ints_1e(mol)
    if relax:
        oqp.library.guess(mol)
        SinglePoint(mol).energy()
        Cd = np.array(mol.data['OQP::VEC_MO_A'], copy=True)
        Cdb = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
        ed = np.array(mol.data['OQP::E_MO_A'], copy=True)
        edb = np.array(mol.data['OQP::E_MO_B'], copy=True)
        # CONJUGATE VARIABLE.  R was measured by perturbing the RAW coefficient
        # matrix in the direction C0[:,p]:  dC[:,q] = eps*C0[:,p].  The variable
        # conjugate to that is therefore V = C0^{-1} C(R) -- NOT the MO overlap.
        # dM/dR = U^x + C^T S^(ket) C carries an extra one-sided overlap term, so
        # contracting R with dM/dR is not a chain rule at all (that was the bug).
        V = np.linalg.solve(C0, Cd)
        # ---- ORBITAL PHASE/ORDER FIX (this was bug #3, and it invalidated every
        # earlier test).  An independent SCF at the displaced geometry returns
        # orbitals with an ARBITRARY per-orbital sign, and near-degenerate ones can
        # come back swapped.  Raw V = C0^-1 C(R) then has O(1) off-diagonal or
        # negative-diagonal entries, so dV/dR ~ O(1/dx) ~ 1e3 (measured: 3.3e3 on
        # every coordinate) and the finite difference is meaningless.
        # Fix: greedy max-|overlap| assignment of displaced -> reference orbitals,
        # then force a positive diagonal.  This is the same continuity problem the
        # production _compute_amp_damp solves with a per-block Loewdin transport.
        perm = [-1] * nbf
        W = np.abs(V).copy()
        for _ in range(nbf):
            i, j = np.unravel_index(W.argmax(), W.shape)   # i=ref, j=displaced
            perm[j] = i
            W[i, :] = -1.0
            W[:, j] = -1.0
        order = np.argsort(perm)              # displaced index -> slot of ref orbital
        Cd, Cdb = Cd[:, order], Cdb[:, order]
        V = np.linalg.solve(C0, Cd)
        sgn = np.sign(np.diag(V))
        sgn[sgn == 0] = 1.0
        Cd, Cdb = Cd * sgn, Cdb * sgn
        M = np.linalg.solve(C0, Cd)
        dev = np.abs(M - np.eye(nbf)).max()
        if dev > 0.2:
            print(f'  WARN: |V-I|={dev:.3e} after phase/order fix (still discontinuous)')
        put_C(Cd, Cdb)
    else:
        put_C(C0, C0b)
        M = None
    freeze_fock()
    return E_all(), M


dE = {p: np.zeros(ncoord) for p in pairs}          # total, moving frame
dExp = {p: np.zeros(ncoord) for p in pairs}        # explicit (C fixed)
Ux = np.zeros((ncoord, nbf, nbf))
for k in range(ncoord):
    d = DX * np.eye(ncoord)[k]
    ep, Mp = at_geom(xyz0 + d, True)
    em, Mm = at_geom(xyz0 - d, True)
    xp, _ = at_geom(xyz0 + d, False)
    xm, _ = at_geom(xyz0 - d, False)
    Ux[k] = (Mp - Mm) / (2.0 * DX)   # = dV/dR, conjugate to R
    for pr in pairs:
        dE[pr][k] = (ep[pr] - em[pr]) / (2.0 * DX)
        dExp[pr][k] = (xp[pr] - xm[pr]) / (2.0 * DX)

cfg['guess'].update(guess_bak)
mol.update_system(xyz0)
oqp.library.ints_1e(mol)
oqp.library.guess(mol)
SinglePoint(mol).energy()
mol.data['OQP::td_bvec_mo'] = X0_raw

# ============================================================ 3. the identity
print('# MOVING-FRAME FROZEN-FOCK identity:  dE/dR  ==  explicit + sum R.U^x')
print(f'# {"pair":>7} {"|dE|":>11} {"|explicit|":>11} {"|R.Ux|":>11} '
      f'{"|resid|":>11} {"resid/|dE|":>11} {"cos(lhs,rhs)":>13}')
out = {}
for pr in pairs:
    rux = np.array([np.sum(R[pr] * Ux[k]) for k in range(ncoord)])
    rhs = dExp[pr] + rux
    res = dE[pr] - rhs
    c = dE[pr] @ rhs / (np.linalg.norm(dE[pr]) * np.linalg.norm(rhs) + 1e-300)
    print(f'  {str((pr[0]+1, pr[1]+1)):>7} {np.linalg.norm(dE[pr]):11.6f} '
          f'{np.linalg.norm(dExp[pr]):11.6f} {np.linalg.norm(rux):11.6f} '
          f'{np.linalg.norm(res):11.6f} '
          f'{np.linalg.norm(res)/(np.linalg.norm(dE[pr])+1e-300):11.4f} {c:+13.6f}')
    out[f'dE_{pr[0]+1}{pr[1]+1}'] = dE[pr]
    out[f'dExp_{pr[0]+1}{pr[1]+1}'] = dExp[pr]
    out[f'RUx_{pr[0]+1}{pr[1]+1}'] = rux
    out[f'R_{pr[0]+1}{pr[1]+1}'] = R[pr]
out['Ux'] = Ux
np.savez(OUT, **out)
print(f'\nsaved -> {OUT}')
