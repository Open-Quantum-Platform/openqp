#!/usr/bin/env python3
"""Full DFT test inputs for openqp-gpu, with a spherical SCF + Cartesian XC bridge.

DF-SCF inputs (B, H/S, guess) are SPHERICAL (well-conditioned, the driver's basis);
the XC grid/basis are CARTESIAN (what the XC kernel needs); c2s.bin maps between
them (AO_sph_i = sum_a c2s[a,i] AO_cart_a). Reference EREF = pyscf spherical DF-RKS.
Writes to <outdir>: B_w{n}.bin H_w{n}.npy S_w{n}.npy guess_w{n}.bin (spherical);
grid.bin basis.bin cart.bin (cartesian); c2s.bin; meta.env.
Usage: prep_dft.py <nwat> <outdir> [func=blyp]
"""
import numpy as np, struct, os, sys, itertools
from pyscf import gto, dft, df, scf
from pyscf.gto import NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF

n=int(sys.argv[1]); OUT=sys.argv[2]; FUNC=sys.argv[3] if len(sys.argv)>3 else "blyp"
os.makedirs(OUT, exist_ok=True)
def water_cluster(n):
    r=0.9572; th=np.deg2rad(104.52); a=2.8; at=[]; c=0
    for i,j,k in itertools.product(range(4),repeat=3):
        if c>=n: break
        o=(i*a,j*a,k*a)
        at+=[("O",o),("H",(o[0]+r,o[1],o[2])),("H",(o[0]+r*np.cos(th),o[1]+r*np.sin(th),o[2]))]; c+=1
    return at
ATOM=water_cluster(n)
mol =gto.M(atom=ATOM,basis="cc-pvdz",unit="Angstrom",cart=False,verbose=0)  # SCF (spherical)
molc=gto.M(atom=ATOM,basis="cc-pvdz",unit="Angstrom",cart=True, verbose=0)  # XC (cartesian)
nsph=mol.nao_nr(); ncart=molc.nao_nr(); nocc=mol.nelectron//2; enuc=mol.energy_nuc()
ntri=nsph*(nsph+1)//2; hfscale=dft.libxc.hybrid_coeff(FUNC)
def wi(f,*v): f.write(struct.pack("<%dq"%len(v),*v))
def wd(f,a): f.write(np.ascontiguousarray(a,dtype="<f8").tobytes())
def pack(M):
    nb=M.shape[0]; o=np.zeros(nb*(nb+1)//2)
    for i in range(nb): b=i*(i+1)//2; o[b:b+i+1]=M[i,:i+1]
    return o
CART=lambda l:[(i,j,l-i-j) for i in range(l,-1,-1) for j in range(l-i,-1,-1)]

# ---- spherical 1e + DF tensor + guess (the SCF side) ----
S=mol.intor("int1e_ovlp"); H=mol.intor("int1e_kin")+mol.intor("int1e_nuc")
np.save(f"{OUT}/H_w{n}.npy",H); np.save(f"{OUT}/S_w{n}.npy",S)
auxmol=df.addons.make_auxmol(mol,"def2-universal-jkfit")
T3=df.incore.aux_e2(mol,auxmol); Vm=auxmol.intor("int2c2e")
w,U=np.linalg.eigh(Vm); keep=w>1e-13*w.max(); Wm=U[:,keep]/np.sqrt(w[keep])
B=np.ascontiguousarray(np.moveaxis(np.tensordot(T3,Wm,axes=(2,0)),2,0)); naux=B.shape[0]
with open(f"{OUT}/B_w{n}.bin","wb") as f: f.write(struct.pack("ii",naux,nsph)); B.astype("<f8").tofile(f)
del T3,B
pack(scf.RHF(mol).get_init_guess(key="minao")).astype("<f8").tofile(open(f"{OUT}/guess_w{n}.bin","wb"))

# ---- c2s: AO_sph_i = sum_a c2s[a,i] AO_cart_a  ({ncart, nsph}) ----
c2s=mol.cart2sph_coeff(normalized="sp")   # (ncart, nsph)
with open(f"{OUT}/c2s.bin","wb") as f: wi(f,ncart,nsph); wd(f,c2s.reshape(-1))

# ---- Cartesian XC grid + basis (mol._env coeffs; general contractions segmented) ----
grids=dft.gen_grid.Grids(molc); grids.level=3; grids.build()
coords,weights=grids.coords,grids.weights; npts=len(weights)
with open(f"{OUT}/grid.bin","wb") as f: wi(f,npts); wd(f,coords.reshape(-1)); wd(f,weights)
am=[];g0=[];nc=[];ao=[];cx=[];cy=[];cz=[];ex=[];cc=[]; aoidx=0
for ish in range(molc.nbas):
    l=molc.bas_angular(ish); npr=molc._bas[ish,NPRIM_OF]; nct=molc._bas[ish,NCTR_OF]
    pe=molc._bas[ish,PTR_EXP]; pc=molc._bas[ish,PTR_COEFF]
    e=molc._env[pe:pe+npr]; coef=molc._env[pc:pc+npr*nct].reshape(nct,npr).T; ctr=molc.bas_coord(ish)
    for jc in range(nct):
        am.append(l); g0.append(len(ex)+1); nc.append(npr); ao.append(aoidx+1)
        cx.append(ctr[0]); cy.append(ctr[1]); cz.append(ctr[2])
        for p in range(npr): ex.append(float(e[p])); cc.append(float(coef[p,jc]))
        aoidx+=(l+1)*(l+2)//2
NSH=len(am); exA=np.array(ex); ccA=np.array(cc)
raw=np.zeros((npts,ncart))
for ish in range(NSH):
    l=am[ish]; p0=g0[ish]-1; a0=ao[ish]-1; c=np.array([cx[ish],cy[ish],cz[ish]])
    d=coords-c; rr=np.einsum("pi,pi->p",d,d); rad=np.zeros(npts)
    for p in range(nc[ish]): rad+=ccA[p0+p]*np.exp(-exA[p0+p]*rr)
    for k,(ix,iy,iz) in enumerate(CART(l)): raw[:,a0+k]=rad*(d[:,0]**ix)*(d[:,1]**iy)*(d[:,2]**iz)
phi_py=molc.eval_gto("GTOval_cart",coords)
bfnrm=np.einsum("p,pi,pi->i",weights,phi_py,raw)/(np.einsum("p,pi,pi->i",weights,raw,raw)+1e-300)
resid=np.sqrt(np.einsum("p,pi->i",weights,(bfnrm*raw-phi_py)**2))/(np.sqrt(np.einsum("p,pi->i",weights,phi_py**2))+1e-300)
with open(f"{OUT}/basis.bin","wb") as f:
    wi(f,NSH,ncart,len(ex))
    for s in range(NSH): wi(f,am[s],0,g0[s],nc[s],ao[s],0); wd(f,[cx[s],cy[s],cz[s]]); wd(f,[1e30])
    wd(f,ex); wd(f,cc); wd(f,[1e30]*len(ex)); wd(f,bfnrm)
maxang=max(am); maxcart=(maxang+1)*(maxang+2)//2
tx=np.zeros(maxcart*(maxang+1),dtype=np.int64); ty=tx.copy(); tz=tx.copy()
for l in range(maxang+1):
    for k,(ix,iy,iz) in enumerate(CART(l)): tx[k+maxcart*l]=ix; ty[k+maxcart*l]=iy; tz[k+maxcart*l]=iz
with open(f"{OUT}/cart.bin","wb") as f: wi(f,maxang,maxcart); f.write(tx.tobytes()); f.write(ty.tobytes()); f.write(tz.tobytes())

# ---- reference: pyscf spherical DF-RKS ----
mf=dft.RKS(mol).density_fit(auxbasis="def2-universal-jkfit"); mf.xc=FUNC; mf.grids=grids; mf.conv_tol=1e-9
E=mf.kernel()
open(f"{OUT}/meta.env","w").write(f"NAO={nsph}\nNOCC={nocc}\nENUC={enuc:.12f}\nNAUX={naux}\n"
                                  f"HFSCALE={hfscale:.6f}\nEREF={E:.10f}\nFUNC={FUNC.upper()}\nNCART={ncart}\n")
print(f"[dft w{n} {FUNC}] nsph={nsph} ncart={ncart} naux={naux} npts={npts} hfscale={hfscale} "
      f"bfnrm-resid={resid.max():.1e} pyscf_sph_DFRKS_E={E:.8f} -> {OUT}")
