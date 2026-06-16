import numpy as np
d=np.load('classfit.npz')
noca=int(d['noca']); nocb=int(d['nocb']); nbf=int(d['nbf']); nvirb=nbf-nocb
lr1=noca-2; lr2=noca-1; sq=1.0/np.sqrt(2.0)
def iatogen(xv):
    w=np.zeros((nbf,nbf)); w[0:noca,nocb:nbf]=xv.reshape(nvirb,noca).T; return w
def comp(w):
    p=np.zeros(noca*nvirb); ij=0
    for j in range(nocb,nbf):
        for i in range(0,noca):
            p[ij]=w[i,j]; ij+=1
    return p
def esum_full(xraw,fij,fab):
    wrk=iatogen(xraw); scr=wrk.copy(); scr[lr1,lr1]=0.0; scr[lr2,lr2]=0.0
    tmp1=np.zeros((nbf,nbf))
    tmp1[0:noca,nocb:nbf]=scr[0:noca,nocb:nbf]@fab[nocb:nbf,nocb:nbf]
    tmp1[0:noca,nocb:nbf]+=-fij[0:noca,0:noca]@scr[0:noca,nocb:nbf]
    xlr=wrk[lr1,lr1]; wrk1=np.zeros((nbf,nbf))
    for j in range(nocb,nbf):
        for i in range(0,noca):
            wrk1[i,j]+=tmp1[i,j]
            if i==lr1: wrk1[i,j]+=fab[j,lr1]*xlr*sq
            if i==lr2: wrk1[i,j]-=fab[j,lr2]*xlr*sq
            if j==lr1: wrk1[i,j]-=fij[i,lr1]*xlr*sq
            if j==lr2: wrk1[i,j]+=fij[i,lr2]*xlr*sq
    dumn=(-fij[lr1,0:noca]@scr[0:noca,lr1]+fij[lr2,0:noca]@scr[0:noca,lr2]
          +fab[lr1,nocb:nbf]@scr[lr1,nocb:nbf]-fab[lr2,nocb:nbf]@scr[lr2,nocb:nbf])
    wrk1[lr1,lr1]=dumn*sq+xlr*(fab[lr1,lr1]+fab[lr2,lr2]-fij[lr1,lr1]-fij[lr2,lr2])*0.5
    wrk1[lr2,lr2]=0.0
    return wrk1
def mrsfxvec(xv):
    ijlr1=(noca-1-nocb-1)*noca+(noca-1)-1; ijlr2=(noca-nocb-1)*noca+noca-1
    out=xv.copy(); val=xv[ijlr1]; out[ijlr1]=val*sq; out[ijlr2]=-val*sq; return out
def Q2bil(xi,xj,fa,fb):
    xfi=mrsfxvec(xi).reshape(nvirb,noca).T; xfj=mrsfxvec(xj).reshape(nvirb,noca).T
    Tij=-0.5*(xfi@xfj.T+xfj@xfi.T); Tab=0.5*(xfi.T@xfj+xfj.T@xfi)
    return np.trace(fa[0:noca,0:noca]@Tij)+np.trace(fb[nocb:nbf,nocb:nbf]@Tab)

for seed in [3,7,11,42]:
    rng=np.random.RandomState(seed)
    A=rng.randn(nbf,nbf); fa=0.5*(A+A.T); B=rng.randn(nbf,nbf); fb=0.5*(B+B.T)
    M=np.zeros((noca*nvirb,noca*nvirb))
    for k in range(noca*nvirb):
        ek=np.zeros(noca*nvirb); ek[k]=1.0; M[:,k]=comp(esum_full(ek,fa,fb))
    N=np.zeros((noca*nvirb,noca*nvirb))
    for a in range(noca*nvirb):
        ea=np.zeros(noca*nvirb); ea[a]=1.0
        for b in range(noca*nvirb):
            eb=np.zeros(noca*nvirb); eb[b]=1.0; N[a,b]=Q2bil(ea,eb,fa,fb)
    print(f'seed={seed}: |M-N|max={np.abs(M-N).max():.2e}  |M-M^T|max={np.abs(M-M.T).max():.2e}  |N-N^T|max={np.abs(N-N.T).max():.2e}')
