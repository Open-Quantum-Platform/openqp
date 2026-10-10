#!/usr/bin/env python3
"""Generate a self-contained, frame-matched XC test case for routec_vxc.

Everything from ONE pyscf Cartesian molecule so the grid, orbitals, density and
reference XC are guaranteed consistent:
  grid.bin  npts:int64, xyz[3*npts]:f64, w[npts]:f64
  basis.bin nshell,nbf,nprim:int64; per shell {am,_,g0,nc,ao,_, cx,cy,cz, md2};
            ex[nprim],cc[nprim],pmd2[nprim],bfnrm[nbf]  (all f64)
  cart.bin  maxang,maxcart:int64; tx,ty,tz[maxcart*(maxang+1)]:int64
  dens.bin  ntri:int64, D_packed[ntri]:f64      (pyscf Cartesian density)
  ref.bin   ntri:int64, Vxc_packed[ntri]:f64, Exc:f64, totele:f64  (pyscf numint)

bfnrm is calibrated NUMERICALLY: bfnrm_i = pyscf_AO_i(r0) / phi_raw_i(r0), so
that bfnrm * (raw collocation AO) == pyscf's AO exactly, whatever the convention.
Usage: make_xc_testcase.py <outdir> [func=svwn]
"""
import numpy as np, struct, os, sys, itertools
from pyscf import gto, dft
from pyscf.gto import NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF

OUT = sys.argv[1]; FUNC = sys.argv[2] if len(sys.argv) > 2 else "svwn"
os.makedirs(OUT, exist_ok=True)
# pyscf functional string per requested name (must match the XC kernel's functional)
PYXC = {"svwn":"svwn", "slater":"slater,", "blyp":"blyp",
        "b3lyp":"b3lypg", "bhhlyp":"bhandhlyp"}[FUNC]

mol = gto.M(atom="O 0 0 0; H 0.9572 0 0; H -0.239987656 0.926627480 0",
            basis="cc-pvdz", unit="Angstrom", cart=True, verbose=0)
nbf = mol.nao_nr(); ntri = nbf*(nbf+1)//2
def wi(f,*v): f.write(struct.pack("<%dq"%len(v),*v))
def wd(f,a): f.write(np.ascontiguousarray(a,dtype="<f8").tobytes())
def pack(M):
    o=np.zeros(ntri)
    for i in range(nbf):
        b=i*(i+1)//2; o[b:b+i+1]=M[i,:i+1]
    return o

# ---- grid ----
grids = dft.gen_grid.Grids(mol); grids.level = 3; grids.build()
coords, weights = grids.coords, grids.weights; npts = len(weights)
with open(f"{OUT}/grid.bin","wb") as f:
    wi(f, npts); wd(f, coords.reshape(-1)); wd(f, weights)

# ---- shells: ex, cc(with primitive norm), centers; bfnrm calibrated to pyscf ----
# bfnrm_i = pyscf_AO_i(probe)/raw_i(probe), probe placed NEAR EACH SHELL'S center
# so both are O(1) (calibrating far from an atom divides two ~0 numbers).
am=[]; g0=[]; nc=[]; ao=[]; cx=[]; cy=[]; cz=[]; md2=[]
ex=[]; cc=[]; pmd2=[]
bfnrm = np.zeros(nbf)
CART = lambda l: [(i,j,l-i-j) for i in range(l,-1,-1) for j in range(l-i,-1,-1)]
off = np.array([0.15, 0.20, 0.25])          # distinct comps so x^i y^j z^k != 0
aoidx = 0
for ish in range(mol.nbas):
    l = mol.bas_angular(ish)
    npr = mol._bas[ish, NPRIM_OF]; nct = mol._bas[ish, NCTR_OF]
    pe = mol._bas[ish, PTR_EXP]; pcoef = mol._bas[ish, PTR_COEFF]
    e = mol._env[pe:pe+npr]                                  # exponents GTOval uses
    coeff = mol._env[pcoef:pcoef+npr*nct].reshape(nct, npr).T  # (nprim,nctr) GTOval uses
    ctr = mol.bas_coord(ish)
    for jc in range(coeff.shape[1]):        # SEGMENT: one shell per contraction column
        c = coeff[:, jc]                    # already for normalized primitives
        am.append(l); g0.append(len(ex)+1); nc.append(len(e))
        ao.append(aoidx+1); cx.append(ctr[0]); cy.append(ctr[1]); cz.append(ctr[2]); md2.append(1e30)
        for p in range(len(e)):
            ex.append(float(e[p])); cc.append(float(c[p])); pmd2.append(1e30)
        aoidx += (l+1)*(l+2)//2
NSH = len(am)                               # total segmented shells
exA=np.array(ex); ccA=np.array(cc)

# bfnrm by least-squares over the WHOLE grid: bfnrm_i = <phi_py_i, raw_i>_w /
# <raw_i, raw_i>_w so bfnrm_i * raw_i best-matches pyscf's AO_i everywhere.
raw = np.zeros((npts, nbf))
for ish in range(NSH):
    l=am[ish]; p0=g0[ish]-1; a0=ao[ish]-1; c=np.array([cx[ish],cy[ish],cz[ish]])
    d=coords-c; rr=np.einsum("pi,pi->p",d,d)
    rad=np.zeros(npts)
    for p in range(nc[ish]): rad+=ccA[p0+p]*np.exp(-exA[p0+p]*rr)
    for k,(ix,iy,iz) in enumerate(CART(l)):
        raw[:,a0+k]=rad*(d[:,0]**ix)*(d[:,1]**iy)*(d[:,2]**iz)
phi_py = mol.eval_gto("GTOval_cart", coords)
num = np.einsum("p,pi,pi->i", weights, phi_py, raw)
den = np.einsum("p,pi,pi->i", weights, raw, raw) + 1e-300
bfnrm = num/den
resid = np.sqrt(np.einsum("p,pi->i", weights, (bfnrm*raw-phi_py)**2)) / \
        (np.sqrt(np.einsum("p,pi->i", weights, phi_py**2))+1e-300)
print(f"  bfnrm fit: max AO shape-residual = {resid.max():.2e} (AO {int(resid.argmax())})")
with open(f"{OUT}/basis.bin","wb") as f:
    wi(f, NSH, nbf, len(ex))
    for s in range(NSH):
        wi(f, am[s], 0, g0[s], nc[s], ao[s], 0); wd(f, [cx[s],cy[s],cz[s]]); wd(f,[md2[s]])
    wd(f, ex); wd(f, cc); wd(f, pmd2); wd(f, bfnrm)

# ---- cart.bin: monomial powers per (component, l) ----
maxang = max(am); maxcart = (maxang+1)*(maxang+2)//2
tx=np.zeros(maxcart*(maxang+1),dtype=np.int64); ty=tx.copy(); tz=tx.copy()
for l in range(maxang+1):
    for k,(ix,iy,iz) in enumerate(CART(l)):
        tx[k+maxcart*l]=ix; ty[k+maxcart*l]=iy; tz[k+maxcart*l]=iz
with open(f"{OUT}/cart.bin","wb") as f:
    wi(f, maxang, maxcart); f.write(tx.tobytes()); f.write(ty.tobytes()); f.write(tz.tobytes())

# ---- density + reference Vxc/Exc/totele (pyscf numint, SAME grid) ----
mf = dft.RKS(mol); mf.xc = PYXC; mf.grids = grids; E = mf.kernel()
D = mf.make_rdm1()
ni = mf._numint
totele, exc_e, vxc = ni.nr_rks(mol, grids, mf.xc, D)
with open(f"{OUT}/dens.bin","wb") as f: wi(f, ntri); wd(f, pack(D))
with open(f"{OUT}/ref.bin","wb") as f:
    wi(f, ntri); wd(f, pack(vxc)); wd(f,[exc_e]); wd(f,[totele])

# self-check: replicate the kernel's density-on-grid from OUR basis.bin arrays
# (phi_i = bfnrm_i * raw_i) and integrate -> must give nelec BEFORE the GPU run.
ex=np.array(ex); cc=np.array(cc); bf=bfnrm
phi=np.zeros((npts,nbf))
for ish in range(NSH):
    l=am[ish]; p0=g0[ish]-1; a0=ao[ish]-1; ctr=np.array([cx[ish],cy[ish],cz[ish]])
    d=coords-ctr; rr=np.einsum("pi,pi->p",d,d)
    rad=np.zeros(npts)
    for p in range(nc[ish]): rad+=cc[p0+p]*np.exp(-ex[p0+p]*rr)
    for k,(ix,iy,iz) in enumerate(CART(l)):
        phi[:,a0+k]=bf[a0+k]*rad*(d[:,0]**ix)*(d[:,1]**iy)*(d[:,2]**iz)
rho=np.einsum("pi,ij,pj->p",phi,D,phi)
selfte=float(np.sum(weights*rho))
# diagnostic: pyscf's own AOs on the grid + which of our AOs differ
phi_py = mol.eval_gto("GTOval_cart", coords)
te_py = float(np.sum(weights*np.einsum("pi,ij,pj->p",phi_py,D,phi_py)))
colnorm = np.sqrt(np.sum((phi-phi_py)**2,axis=0)) / (np.sqrt(np.sum(phi_py**2,axis=0))+1e-30)
bad = np.argsort(-colnorm)[:6]
print(f"[testcase {FUNC}] nbf(cart)={nbf} npts={npts} nelec={mol.nelectron} "
      f"RKS_E={E:.8f} Exc_ref={exc_e:.8f} totele_ref={totele:.6f}")
print(f"  totele: ours={selfte:.4f}  pyscf-phi={te_py:.4f}  target={mol.nelectron}")
print(f"  worst AO rel-diffs: " + ", ".join(f"ao{int(i)}(am{am[np.searchsorted([a-1 for a in ao],i,side='right')-1]}):{colnorm[i]:.2f}" for i in bad))
print(f"  ratio ours/pyscf per worst AO: " + ", ".join(f"{(phi[:,i]/(phi_py[:,i]+1e-30))[np.argmax(np.abs(phi_py[:,i]))]:.3f}" for i in bad))
