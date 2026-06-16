import numpy as np
d=np.load('classfit.npz')
X0=d['X0']; noca=int(d['noca']); nocb=int(d['nocb']); nbf=int(d['nbf']); nvirb=nbf-nocb
lr1=noca-2; lr2=noca-1; sq=1.0/np.sqrt(2.0)
rng=np.random.RandomState(7)
A=rng.randn(nbf,nbf); fa0=0.5*(A+A.T); B=rng.randn(nbf,nbf); fb0=0.5*(B+B.T)

def mrsfxvec_mat(col):
    f=col.copy(); ijlr1=(noca-1-nocb-1)*noca+(noca-1)-1; ijlr2=(noca-nocb-1)*noca+noca-1
    f[ijlr1]=col[ijlr1]*sq; f[ijlr2]=-col[ijlr1]*sq
    return f.reshape((nvirb,noca)).T  # (noca,nvirb)

def Trelax(I,J,fa,fb):
    xi=mrsfxvec_mat(X0[:,I]); xj=mrsfxvec_mat(X0[:,J])
    Tij=-0.5*(xi@xj.T+xj@xi.T); Tab=0.5*(xi.T@xj+xj.T@xi)
    return np.trace(fa[0:noca,0:noca]@Tij)+np.trace(fb[nocb:nbf,nocb:nbf]@Tab), Tij, Tab

# sfrorhs 2FT projection (real orbital indices)
def project(hxa,hxb):
    prs=([(i,j) for i in range(nocb,noca) for j in range(0,nocb)]
        +[(k,j) for k in range(noca,nbf) for j in range(0,nocb)]
        +[(k,i) for k in range(noca,nbf) for i in range(nocb,noca)])
    rhs=np.zeros(len(prs)); n=0
    for i in range(nocb,noca):
        for j in range(0,nocb):
            rhs[n]=hxa[i,j]-hxa[j,i]-hxb[j,i]; n+=1
    for k in range(noca,nbf):
        for j in range(0,nocb):
            rhs[n]=hxa[k,j]-hxb[j,k]; n+=1
    for k in range(noca,nbf):
        for i in range(nocb,noca):
            rhs[n]=hxa[k,i]+hxb[k,i]-hxb[i,k]; n+=1
    return -rhs
def L_2ft(I,J,fa,fb):
    _,Tij,Tab=Trelax(I,J,fa,fb)
    xhxa=2.0*fa[:,0:noca]@Tij
    wrk=np.zeros((nbf,nbf)); wrk[nocb:,nocb:]=Tab
    xhxb=2.0*fb@wrk
    return project(xhxa,xhxb)

# Now: is sfrorhs's L_2ft == the true theta-gradient of Trelax?
# theta_pq rotates the orbitals -> fa,fb transform. In FROZEN-FOCK (skeleton): F_AO fixed,
# only C differentiated -> fa = C^T F C, dfa/dtheta_pq = [fa, kappa]_pq style commutator.
# For a single rotation generator E_pq-E_qp: fa -> fa + t*(kappa^T fa + fa kappa) with kappa
# the antisymmetric generator. d Trelax/dtheta_pq = Tr( (dfa/dtheta) Tij ) + Tr((dfb/dtheta) Tab).
# Build this analytically and compare to sfrorhs L_2ft.
def theta_grad(I,J,fa,fb):
    _,Tij,Tab=Trelax(I,J,fa,fb)
    # extend Tij,Tab to full nbf x nbf in proper blocks
    Tij_f=np.zeros((nbf,nbf)); Tij_f[0:noca,0:noca]=Tij
    Tab_f=np.zeros((nbf,nbf)); Tab_f[nocb:,nocb:]=Tab
    prs=([(i,j) for i in range(nocb,noca) for j in range(0,nocb)]
        +[(k,j) for k in range(noca,nbf) for j in range(0,nocb)]
        +[(k,i) for k in range(noca,nbf) for i in range(nocb,noca)])
    g=np.zeros(len(prs))
    for n,(p,q) in enumerate(prs):
        # generator kappa = E_pq - E_qp ; dfa = kappa^T fa + fa kappa  (rotation of fa, alpha for both? )
        # alpha rotation affects fa; beta rotation affects fb. The rotation blocks:
        # doc-socc & doc-virt & soc-virt are ALPHA orbital rotations in this MRSF convention?
        # Test both: build dF from kappa on a generic F.
        K=np.zeros((nbf,nbf)); K[p,q]=1.0; K[q,p]=-1.0
        dfa=K.T@fa+fa@K; dfb=K.T@fb+fb@K
        g[n]=np.trace(dfa@Tij_f)+np.trace(dfb@Tab_f)
    return g

for (I,J) in [(0,1),(0,2),(1,2)]:
    L=L_2ft(I,J,fa0,fb0)
    G=theta_grad(I,J,fa0,fb0)
    c=np.dot(L,G)/(np.linalg.norm(L)*np.linalg.norm(G)+1e-30)
    print(f'pair({I+1},{J+1}): |L_2ft|={np.linalg.norm(L):.4f} |theta_grad|={np.linalg.norm(G):.4f} cos={c:+.6f} ratio={np.linalg.norm(L)/np.linalg.norm(G):.4f} maxdiff={np.abs(L-G).max():.3e}')
