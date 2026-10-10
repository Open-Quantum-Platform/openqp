#!/usr/bin/env python3
"""Drive the openqp-gpu RKS DFT SCF and compare/ time vs the pyscf reference.
Env: LIB (driver .so), ART (prep_dft.py output dir), NW (n water).
Reads meta.env for NOCC/ENUC/HFSCALE/EREF/FUNC, sets the XC grid + B, runs
routec_scf_solve with se=HFSCALE.  N warm reps, median loop wall."""
import os, ctypes, time, numpy as np
LIB=os.environ["LIB"]; ART=os.environ["ART"]; NW=os.environ["NW"]
NREP=int(os.environ.get("NREP","5"))
meta={}
for ln in open(f"{ART}/meta.env"):
    k,v=ln.strip().split("="); meta[k]=v
nocc=int(meta["NOCC"]); enuc=float(meta["ENUC"]); hfscale=float(meta["HFSCALE"])
eref=float(meta["EREF"]); func=meta["FUNC"]
H=np.load(f"{ART}/H_w{NW}.npy"); S=np.load(f"{ART}/S_w{NW}.npy"); nbf=H.shape[0]
ntri=nbf*(nbf+1)//2
def pk(M):
    o=np.zeros(ntri)
    for i in range(nbf): b=i*(i+1)//2; o[b:b+i+1]=M[i,:i+1]
    return np.ascontiguousarray(o)
hpk,spk=pk(H),pk(S)
os.environ["OQP_ROUTEC_B"]=f"{ART}/B_w{NW}.bin"
os.environ["OQP_OWNXC_DIR"]=ART; os.environ["OQP_OWNXC_FUNC"]=func
os.environ["OQP_OWNXC_C2S"]=f"{ART}/c2s.bin"   # spherical SCF <-> cartesian XC bridge
os.environ["OQP_SCF_GUESS_D"]=f"{ART}/guess_w{NW}.bin"
lib=ctypes.CDLL(LIB); fn=lib.routec_scf_solve; fn.restype=None
P=ctypes.POINTER(ctypes.c_double); Pi=ctypes.POINTER(ctypes.c_int)
fn.argtypes=[P]*2+[Pi]*2+[P]*5+[Pi]+[P]*4+[Pi]*2
dp=lambda x:x.ctypes.data_as(P)
def solve():
    d=np.zeros(ntri); c=np.zeros(nbf*nbf); eps=np.zeros(nbf)
    nb=ctypes.c_int(nbf); no=ctypes.c_int(nocc); en=ctypes.c_double(enuc)
    sc=ctypes.c_double(1.0); se=ctypes.c_double(hfscale)
    ce=ctypes.c_double(1e-8); cd=ctypes.c_double(1e-7); mx=ctypes.c_int(100)
    e=ctypes.c_double(0); nc=ctypes.c_int(0); info=ctypes.c_int(1)
    t=time.time()
    # signature order after maxit is e_out, d_out, c_out, eps_out (energy FIRST)
    fn(dp(hpk),dp(spk),ctypes.byref(nb),ctypes.byref(no),ctypes.byref(en),
       ctypes.byref(sc),ctypes.byref(se),ctypes.byref(ce),ctypes.byref(cd),ctypes.byref(mx),
       ctypes.byref(e),dp(d),dp(c),dp(eps),ctypes.byref(nc),ctypes.byref(info))
    return e.value, nc.value, info.value, time.time()-t
E,ncyc,info,_=solve()
walls=[solve()[3] for _ in range(NREP)]
walls.sort()
print(f"[ours DFT w{NW} {func}] E={E:.8f}  ref={eref:.8f}  dE={abs(E-eref):.2e}  "
      f"cyc={ncyc} info={info}  loop_med={walls[len(walls)//2]*1e3:.1f}ms")
