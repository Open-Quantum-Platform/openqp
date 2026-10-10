import sys, numpy as np
def load_v1(p):
    with open(p,'rb') as f:
        h=np.fromfile(f,np.int32,2); B=np.fromfile(f,np.float64)
    a,b=int(h[0]),int(h[1])
    if a>0 and b<0: naux,nao=a,-b
    elif a<0 and b>0: nao,naux=-a,b
    else: raise SystemExit(f"v1 hdr {a},{b}")
    npair=nao*(nao+1)//2; return naux,nao,B.reshape(naux,npair)
def load_v2(p):
    with open(p,'rb') as f:
        magic=int(np.fromfile(f,np.int32,1)[0]); assert magic==0x32464443, hex(magic)
        naux,nao,ncp=[int(x) for x in np.fromfile(f,np.int32,3)]
        keep=np.fromfile(f,np.int32,ncp); cB=np.fromfile(f,np.float64).reshape(naux,ncp)
    return naux,nao,keep,cB
naux,nao,B1=load_v1(sys.argv[1]); n2,nao2,keep,cB=load_v2(sys.argv[2])
npair=nao*(nao+1)//2; full=np.zeros((naux,npair)); full[:,keep]=cB
d=float(np.abs(full-B1).max())
dropped=np.ones(npair,bool); dropped[keep]=False
maxdrop=float(np.abs(B1[:,dropped]).max()) if dropped.any() else 0.0
print(f"v1(naux={naux},npair={npair})  v2 ncp={len(keep)}/{npair}={len(keep)/npair:.4f}  saved={100*(1-len(keep)/npair):.1f}%")
print(f"round-trip max|v2_recon - v1| = {d:.3e}  (must be 0)")
print(f"max|v1 dropped cols|          = {maxdrop:.3e}  (must be 0 = lossless)")
print("PASS" if d==0.0 and maxdrop==0.0 else "FAIL")
