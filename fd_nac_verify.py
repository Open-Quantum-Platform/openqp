#!/usr/bin/env python
"""Verify the factor-2/sign in OpenQP's FD NAC path on a single LiH coordinate.

Approach: run `runtype=nacme` twice (S+ and S-), each computing the MRSF state
overlap td_states_overlap(i,j) = <Psi_i(origin) | Psi_j(displaced)> (orientation per
get_states_overlap.F90:915, get_dcv comment). Then compute, per the user's spec:

  (1) hand    d_ij     = (S+(i,j) - S-(i,j)) / (2 dR)
  (2) current dcv_ij   = [ (S+ - S+^T) - (S- - S-^T) ]_ij / (2 dR)   (== numerical_nac net)
  (3) corrected_ij     = same / (4 dR)                              (HST 1/2 restored)

and report ratios + sign/symmetry diagnostics.
"""
import os, sys, numpy as np

np.set_printoptions(precision=6, suppress=True, linewidth=140)

from oqp.pyoqp import Runner

WORK = "/tmp/fd_nac_lih"
os.makedirs(WORK, exist_ok=True)

ANG2BOHR = 1.8897259886
DR_ANG = 0.01                 # one-coordinate displacement of H along z (Angstrom)
DR_BOHR = DR_ANG * ANG2BOHR   # NAC is per Bohr
LI_Z, H_Z = 0.0, 1.595        # LiH near-equilibrium (Angstrom), far from avoided crossing
NSTATE = 4

INP = """[input]
system=
 3   0.000000000   0.000000000   {liz:.9f}
 1   0.000000000   0.000000000   {cur_hz:.9f}
system2=
 3   0.000000000   0.000000000   {liz:.9f}
 1   0.000000000   0.000000000   {prev_hz:.9f}
charge=0
method=tdhf
basis=6-31g*
runtype={runtype}
functional=bhhlyp
d4=False

[guess]
type=huckel
save_mol=False

[scf]
type=rohf
maxit=100
maxdiis=5
multiplicity=3
conv=1.0e-9
save_molden=False

[dftgrid]
rad_npts=96
ang_npts=302

[properties]
grad=1

[nac]
align=reorder

[tdhf]
type=mrsf
maxit=100
multiplicity=1
conv=1.0e-9
nstate={nstate}
zvconv=1.0e-9
"""

def run(tag, runtype, cur_hz, prev_hz):
    inp = os.path.join(WORK, f"lih_{tag}.inp")
    log = os.path.join(WORK, f"lih_{tag}.log")
    with open(inp, "w") as f:
        f.write(INP.format(liz=LI_Z, cur_hz=cur_hz, prev_hz=prev_hz,
                           runtype=runtype, nstate=NSTATE))
    r = Runner(project=f"lih_{tag}", input_file=inp, log=log, silent=1)
    r.run()
    return r.mol

# --- origin energies (clean reference gap) ---
mol0 = run("origin", "energy", H_Z, H_Z)
E = np.array(mol0.energies)            # [E_scf, E_1, E_2, ...]
print("# origin absolute energies (Ha):", E)

# --- S+ : current = origin + dR ; previous = origin ---
molp = run("plus", "nacme", H_Z + DR_ANG, H_Z)
Sp = np.array(molp.data["OQP::td_states_overlap"])

# --- S- : current = origin - dR ; previous = origin ---
molm = run("minus", "nacme", H_Z - DR_ANG, H_Z)
Sm = np.array(molm.data["OQP::td_states_overlap"])

n = Sp.shape[0]
print(f"\n# nstate matrix dim = {n};  dR = {DR_ANG} Ang = {DR_BOHR:.6f} Bohr")
print("\n# S+  = <Psi_i(origin)|Psi_j(origin+dR)>\n", Sp)
print("\n# S-  = <Psi_i(origin)|Psi_j(origin-dR)>\n", Sm)
print("\n# diag(S+) (should be ~1, no reorder):", np.diag(Sp))
print("# diag(S-) (should be ~1, no reorder):", np.diag(Sm))

# --- the three quantities ---
hand      = (Sp - Sm) / (2.0 * DR_BOHR)                       # (1)
antisym_p = Sp - Sp.T
antisym_m = Sm - Sm.T
current   = (antisym_p - antisym_m) / (2.0 * DR_BOHR)          # (2) numerical_nac net
corrected = (antisym_p - antisym_m) / (4.0 * DR_BOHR)          # (3) HST 1/2 restored

def sym_diag(name, M):
    asym = np.linalg.norm(M + M.T) / (np.linalg.norm(M) + 1e-30)
    sym  = np.linalg.norm(M - M.T) / (np.linalg.norm(M) + 1e-30)
    print(f"  {name:10s}  ||M+M^T||/||M||={asym:.2e} (antisym test)   "
          f"||M-M^T||/||M||={sym:.2e} (sym test)")

print("\n# symmetry diagnostics:")
sym_diag("hand",      hand)
sym_diag("current",   current)
sym_diag("corrected", corrected)

print("\n# per-pair table (i<j), values are the coupling scalar for this 1 coordinate:")
print(f"{'pair':>6} {'hand d_ij':>14} {'current dcv':>14} {'corrected':>14} "
      f"{'cur/hand':>10} {'corr/hand':>10} {'dE=Ej-Ei':>12}")
for i in range(n):
    for j in range(i + 1, n):
        h_ij = hand[i, j]
        c_ij = current[i, j]
        k_ij = corrected[i, j]
        dE = E[j + 1] - E[i + 1]      # E indexed with SCF ref at [0]
        rc = c_ij / h_ij if abs(h_ij) > 1e-12 else float('nan')
        rk = k_ij / h_ij if abs(h_ij) > 1e-12 else float('nan')
        print(f"{i+1}->{j+1:<3} {h_ij:14.6f} {c_ij:14.6f} {k_ij:14.6f} "
              f"{rc:10.4f} {rk:10.4f} {dE:12.6f}")

# explicit raw off-diagonal to expose orientation/sign
print("\n# raw off-diagonals for pair (1,2):")
print(f"  S+[0,1]=<1o|2+>={Sp[0,1]: .6f}   S+[1,0]=<2o|1+>={Sp[1,0]: .6f}")
print(f"  S-[0,1]=<1o|2->={Sm[0,1]: .6f}   S-[1,0]=<2o|1->={Sm[1,0]: .6f}")
print(f"  hand d_12 = (S+[0,1]-S-[0,1])/(2dR) = {hand[0,1]: .6f}")
print(f"  hand d_21 = (S+[1,0]-S-[1,0])/(2dR) = {hand[1,0]: .6f}  (antisym check vs -d_12={-hand[0,1]: .6f})")
print("\nDONE")
