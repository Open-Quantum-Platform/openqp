#!/usr/bin/env python3
"""A2-FP32 prep: dump shells/aux for the C++ bench AND dump the frozen
device-side transform operands (aux-conv Cag, orbital perm+scale, Coulomb whiten
Linv, lower-tri indices) as raw binaries so the single-process C++ B-build can
load them and do assembly+whiten on device with NO Python on the critical path.

Writes into <workdir>:
  shells_c_w{n}.txt, aux_c_w{n}.txt     (3c kernel inputs)
  meta_w{n}.bin     : int32 nao,naux,npair
  cag_w{n}.bin      : float64 (naux,naux) aux-conv = blockdiag inv(M[l_aux])  [row-major]
  perm_w{n}.bin     : int32 (nao,)   orbital perm
  dscale_w{n}.bin   : float64 (nao,) orbital scale
  linv_w{n}.bin     : float64 (naux,naux) frozen Cholesky inverse (lower) row-major
  iu_w{n}.bin,ju_w{n}.bin : int32 (npair,) lower-tri row/col indices (i>=j)

Usage: python a2fp32_prep.py <nwat> <workdir>"""
import sys, os, itertools, numpy as np, scipy.linalg as sla
GEN=os.path.expanduser("~/routec_gen_20260610")
sys.path.insert(0, GEN)
from routec_pyscf_calibration import M as CAL
from pyscf import gto, df
from pyscf.gto.mole import gto_norm

def water_cluster(n_edge):
    r=0.9572; th=np.deg2rad(104.52); a=2.8; atoms=[]; cnt=0
    for (i,j,k) in itertools.product(range(4),repeat=3):
        if cnt>=n_edge: break
        ox,oy,oz=i*a,j*a,k*a
        atoms+=[('O',(ox,oy,oz)),('H',(ox+r,oy,oz)),('H',(ox+r*np.cos(th),oy+r*np.sin(th),oz))]
        cnt+=1
    return atoms

nwat=int(sys.argv[1]); wd=sys.argv[2]; os.makedirs(wd,exist_ok=True)
mol=gto.M(atom=water_cluster(nwat),basis='cc-pvdz',unit='Angstrom',cart=False,verbose=0)
auxmol=df.addons.make_auxmol(mol,'def2-universal-jkfit')
nao=mol.nao; naux=auxmol.nao; npair=nao*(nao+1)//2

def dump(m, path):
    lines, nsh = [], 0
    for ib in range(m.nbas):
        l=m.bas_angular(ib); exps=m.bas_exp(ib); cc=m.bas_ctr_coeff(ib)
        x,y,z=m.bas_coord(ib)
        for col in range(cc.shape[1]):
            sel=np.abs(cc[:,col])>0
            es,cs=exps[sel],cc[sel,col]*gto_norm(l,exps[sel])
            lines.append(f"{l} {x:.12f} {y:.12f} {z:.12f} {len(es)}")
            for e,c in zip(es,cs): lines.append(f"{e:.12g} {c:.12g}")
            nsh+=1
    open(path,"w").write(f"{nsh}\n"+"\n".join(lines)+"\n")
dump(mol,f"{wd}/shells_c_w{nwat}.txt"); dump(auxmol,f"{wd}/aux_c_w{nwat}.txt")

def Mi(l):
    n=2*l+1; Mm=np.array(CAL[l]) if l<=2 else np.sqrt(4*np.pi/(2*l+1))*np.eye(n); return np.linalg.inv(Mm)
def bd(m):
    n=m.nao; C=np.zeros((n,n)); o=0
    for ib in range(m.nbas):
        l=m.bas_angular(ib); k=2*l+1
        for _ in range(m.bas_nctr(ib)): C[o:o+k,o:o+k]=Mi(l); o+=k
    assert o==n; return C
# orbital perm+scale (inv(M[l]) is scalar*perm per row)
perm=np.arange(nao,dtype=np.int32); dscale=np.ones(nao); o=0
for ib in range(mol.nbas):
    l=mol.bas_angular(ib); k=2*l+1; Mi_=Mi(l)
    for _ in range(mol.bas_nctr(ib)):
        for a in range(k):
            row=Mi_[a]; j=int(np.argmax(np.abs(row)))
            perm[o+a]=o+j; dscale[o+a]=row[j]
        o+=k
assert o==nao
# aux perm+scale (inv(M[l]) is scalar*perm per row) -- folds the aux-conv into the
# pack kernel, eliminating the dense block-diagonal aux-conv GEMM (mostly zeros).
aperm=np.arange(naux,dtype=np.int32); ascale=np.ones(naux); o=0
for ib in range(auxmol.nbas):
    l=auxmol.bas_angular(ib); k=2*l+1; Mi_=Mi(l)
    for _ in range(auxmol.bas_nctr(ib)):
        for a in range(k):
            row=Mi_[a]; j=int(np.argmax(np.abs(row)))
            aperm[o+a]=o+j; ascale[o+a]=row[j]
        o+=k
assert o==naux
Cag=bd(auxmol)                                            # (naux,naux)
# verify Cag == perm+scale (so the GEMM-free aux-conv is exact)
Cchk=np.zeros((naux,naux)); Cchk[np.arange(naux),aperm]=ascale
aux_permscale_err=np.abs(Cchk-Cag).max()
V=auxmol.intor('int2c2e',hermi=1); L=sla.cholesky(V,lower=True)
Linv=sla.solve_triangular(L,np.eye(naux),lower=True)      # (naux,naux)
iu,ju=np.tril_indices(nao)                                # i>=j, len npair

np.array([nao,naux,npair],dtype=np.int32).tofile(f"{wd}/meta_w{nwat}.bin")
np.ascontiguousarray(Cag,dtype=np.float64).tofile(f"{wd}/cag_w{nwat}.bin")
perm.tofile(f"{wd}/perm_w{nwat}.bin")
np.ascontiguousarray(dscale,dtype=np.float64).tofile(f"{wd}/dscale_w{nwat}.bin")
aperm.tofile(f"{wd}/aperm_w{nwat}.bin")
np.ascontiguousarray(ascale,dtype=np.float64).tofile(f"{wd}/ascale_w{nwat}.bin")
np.ascontiguousarray(Linv,dtype=np.float64).tofile(f"{wd}/linv_w{nwat}.bin")
iu.astype(np.int32).tofile(f"{wd}/iu_w{nwat}.bin")
ju.astype(np.int32).tofile(f"{wd}/ju_w{nwat}.bin")
print(f"w{nwat}: nao={nao} naux={naux} npair={npair} aux_permscale_err={aux_permscale_err:.1e} prepared in {wd}",flush=True)
