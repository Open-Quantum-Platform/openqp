"""BP#1 split: separate the missing damp term into
  (1eW) explicit nuclear derivative of the matvec at FROZEN orbitals C0
        (only AO Fock + 1e/2e ints move; the orbital coeffs are NOT updated,
         no transport) -> the Hellmann-Feynman-like 1e/W + ERI explicit deriv.
  (Ux)  the remainder = oracle - (frozen explicit FD damp).

We FD the full matvec with orbitals pinned to C0 (we rebuild the AO Fock at each
displaced geometry from the FROZEN density D0 so the matvec sees the moved 1e/2e
operators but unchanged orbitals).  Comparing:
  full(oracle)            = transported full-matvec FD  (= eigenvector deriv)
  frozenC explicit FD     = dA/dR at fixed C0 (no orbital response)
  ana2e (Fortran)         = explicit ERI 2e part only
tells us how much of `miss` is the explicit mrsfesum(Fock)-derivative (1e/W) vs
the orbital-response Ux.
"""
import os
import math
import oqp
import oqp.library
from oqp.pyoqp import Runner
from oqp.library.single_point import SinglePoint
import numpy as np

INP = '/tmp/nactest/H2O_tight_dx0.001.inp'
ORACLE = '/tmp/nactest/p11_damp_oracle.npz'
DELTA = 1e-3

r = Runner(input_file=INP, log='/tmp/nactest/p11_split.log'); r.run(); mol = r.mol
mol.data._data.control.int2e_cutoff = 1e-20
nstate = 3
noca = int(np.asarray(mol.data['nelec_A']).ravel()[0]); nocb = noca-2
C0raw = np.array(mol.data['OQP::VEC_MO_A'], copy=True); nbf = C0raw.shape[0]
nvirb = nbf-nocb; nij = noca*nvirb
C0b = np.array(mol.data['OQP::VEC_MO_B'], copy=True)
e0 = np.array(mol.data['OQP::E_MO_A'], copy=True)
e0b = np.array(mol.data['OQP::E_MO_B'], copy=True)
DMA0 = np.array(mol.data['OQP::DM_A'], copy=True)
DMB0 = np.array(mol.data['OQP::DM_B'], copy=True)
FA0 = np.array(mol.data['OQP::FOCK_A'], copy=True)
FB0 = np.array(mol.data['OQP::FOCK_B'], copy=True)
natom = mol.data['natom']; ncoord = 3*natom
xyz0 = np.array(mol.get_system(), copy=True)
X0_raw = np.array(mol.data['OQP::td_bvec_mo'], copy=True); Xshape = X0_raw.shape
X0 = X0_raw.reshape(-1).reshape((nstate, nij))
E = list(mol.energies); Om = [E[k+1]-E[0] for k in range(nstate)]

def set_bvec(col):
    raw = X0_raw.copy().reshape(-1); raw[0:nij] = col
    mol.data['OQP::td_bvec_mo'] = raw.reshape(Xshape)

def Acol(c):
    set_bvec(c); oqp.mrsf_matvec_apply(mol)
    return np.array(mol.data['OQP::nac_mvax'], copy=True).ravel()

def Amat():
    A = np.zeros((nij, nij)); e = np.zeros(nij)
    for j in range(nij):
        e[:] = 0.0; e[j] = 1.0; A[:, j] = Acol(e)
    return A

mol.save_data()

def frozenC_A(coord):
    """matvec at displaced geometry but FROZEN orbitals C0 and FROZEN density D0:
    rebuild the AO Fock from D0 at the new geometry (Hartree+XC ints move), keep
    the MO coefficients = C0. No SCF, no transport."""
    mol.update_system(coord)
    oqp.library.ints_1e(mol)
    # pin the orbitals + density to the reference
    mol.data['OQP::VEC_MO_A'] = C0raw; mol.data['OQP::VEC_MO_B'] = C0b
    mol.data['OQP::E_MO_A'] = e0;      mol.data['OQP::E_MO_B'] = e0b
    mol.data['OQP::DM_A'] = DMA0;      mol.data['OQP::DM_B'] = DMB0
    # rebuild the AO Fock at the new geometry from the FROZEN density
    oqp.library.hf_energy(mol)   # builds FOCK_A/B from current DM at new ints
    return Amat()

print("building frozen-C explicit dA...")
dA_fz = np.zeros((ncoord, nij, nij))
for k in range(ncoord):
    Ap = frozenC_A(xyz0+DELTA*np.eye(ncoord)[k])
    Am = frozenC_A(xyz0-DELTA*np.eye(ncoord)[k])
    dA_fz[k] = (Ap-Am)/(2*DELTA)
mol.update_system(xyz0); oqp.library.ints_1e(mol)
mol.data['OQP::VEC_MO_A'] = C0raw; mol.data['OQP::DM_A'] = DMA0
mol.data['OQP::FOCK_A'] = FA0; mol.data['OQP::FOCK_B'] = FB0

oracle = np.load(ORACLE)
def cs(a, b): return float(np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-30))
print("\n=== frozen-C explicit damp vs oracle (Ux = oracle - frozenC) ===")
for (I, J) in [(0, 1), (0, 2), (1, 2)]:
    gap = Om[J]-Om[I]
    fz = np.array([X0[I] @ dA_fz[k] @ X0[J] for k in range(ncoord)])/gap
    o = oracle[f'{I+1}{J+1}_damp']
    ux = o - fz
    print(f"\n({I+1},{J+1}): |oracle|={np.linalg.norm(o):.5f}  |frozenC|={np.linalg.norm(fz):.5f}  cos(frozenC,oracle)={cs(fz,o):+.5f}  ratio={np.linalg.norm(fz)/(np.linalg.norm(o)+1e-30):.3f}")
    print(f"  |Ux=o-frozenC|={np.linalg.norm(ux):.5f}  cos(Ux,oracle)={cs(ux,o):+.5f}")
    print(f"  oracle ={np.array2string(o, precision=5, suppress_small=True)}")
    print(f"  frozenC={np.array2string(fz, precision=5, suppress_small=True)}")
mol.update_system(xyz0)
