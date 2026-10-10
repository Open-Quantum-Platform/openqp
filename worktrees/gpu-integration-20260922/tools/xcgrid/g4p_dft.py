#!/usr/bin/env python3
"""gpu4pyscf DF-RKS timing bar (spherical, same molecule/func/aux as ours)."""
import sys, time, itertools, numpy as np, cupy as cp
from pyscf import gto
from gpu4pyscf import dft as gdft
def water_cluster(n):
    r=0.9572; th=np.deg2rad(104.52); a=2.8; at=[]; c=0
    for i,j,k in itertools.product(range(4),repeat=3):
        if c>=n: break
        o=(i*a,j*a,k*a)
        at+=[("O",o),("H",(o[0]+r,o[1],o[2])),("H",(o[0]+r*np.cos(th),o[1]+r*np.sin(th),o[2]))]; c+=1
    return at
nw=int(sys.argv[1]); FUNC=sys.argv[2] if len(sys.argv)>2 else "blyp"; NREP=int(sys.argv[3]) if len(sys.argv)>3 else 4
mol=gto.M(atom=water_cluster(nw),basis="cc-pvdz",unit="Angstrom",cart=False,verbose=0)
def sync(): cp.cuda.runtime.deviceSynchronize()
_a=cp.random.rand(4096,4096)
for _ in range(40): _c=_a@_a
sync()
loops=[]; E=None; nit=None
for rep in range(NREP+1):
    mf=gdft.RKS(mol).density_fit(auxbasis="def2-universal-jkfit"); mf.xc=FUNC
    mf.conv_tol=1e-8; mf.max_cycle=100; mf.verbose=0
    sync(); t0=time.time(); E=mf.kernel(); sync(); t1=time.time()
    nit=mf.cycles if hasattr(mf,"cycles") else None
    if rep>0: loops.append((t1-t0)*1e3)
loops.sort()
print(f"[g4p DFT w{nw} {FUNC}] E={E:.8f} iters={nit} loop_med={loops[len(loops)//2]:.1f}ms  reps={[round(x,1) for x in loops]}")
