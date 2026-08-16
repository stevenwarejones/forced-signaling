import numpy as np
from scipy.optimize import linprog
from itertools import product

w3 = np.exp(2j*np.pi/3)
def fourier_basis(alpha):
    """basis vectors f_k(alpha) = (1/sqrt3) sum_j w^(j(k+alpha)) |j>"""
    B=[]
    for k in range(3):
        v=np.array([np.exp(2j*np.pi*j*(k+alpha)/3) for j in range(3)])/np.sqrt(3)
        B.append(v)
    return B
comp_basis=[np.eye(3)[:,k].astype(complex) for k in range(3)]

def qutrit_cluster():
    """4-qutrit linear cluster: CZ3 chain on |+3>^4, order A,B,C,D. Return 3x3x3x3 tensor."""
    psi=np.ones((3,3,3,3),dtype=complex)/9.0
    for (i,j) in [(0,1),(1,2),(2,3)]:
        for idxs in product(range(3),repeat=4):
            psi[idxs]*=w3**(idxs[i]*idxs[j])
    return psi

def sigma_hic_qutrit(basesA,basesB,basesC,basesD,psi):
    # quantum marginals
    def probs(x,y,z,w):
        vA=basesA[x]; vB=basesB[y]; vC=basesC[z]; vD=basesD[w]
        P=np.zeros((3,3,3,3))
        for a,b,c,d in product(range(3),repeat=4):
            amp=np.einsum('jklm,j,k,l,m->',psi,vA[a].conj(),vB[b].conj(),vC[c].conj(),vD[d].conj())
            P[a,b,c,d]=abs(amp)**2
        return P
    PQ={(x,y,z,w):probs(x,y,z,w) for x,y,z,w in product(range(2),repeat=4)}
    # model variables
    idx={}; k=0
    for x,ww,a,d,fb,fg in product(range(2),range(2),range(3),range(3),range(9),range(9)):
        idx[(x,ww,a,d,fb,fg)]=k; k+=1
    N=k
    fB=lambda fn,y: fn//3 if y==0 else fn%3
    fC=lambda fn,z: fn//3 if z==0 else fn%3
    eq_rows=[]; eq_rhs=[]
    for x,ww in product(range(2),range(2)):
        r=np.zeros(N)
        for a,d,fb,fg in product(range(3),range(3),range(9),range(9)):
            r[idx[(x,ww,a,d,fb,fg)]]=1
        eq_rows.append(r); eq_rhs.append(1.0)
    for x,y,ww,a,b,d in product(range(2),range(2),range(2),range(3),range(3),range(3)):
        r=np.zeros(N)
        for fb in range(9):
            if fB(fb,y)!=b: continue
            for fg in range(9): r[idx[(x,ww,a,d,fb,fg)]]+=1
        eq_rows.append(r); eq_rhs.append(PQ[(x,y,0,ww)][a,b,:,d].sum())
    for x,z,ww,a,c,d in product(range(2),range(2),range(2),range(3),range(3),range(3)):
        r=np.zeros(N)
        for fg in range(9):
            if fC(fg,z)!=c: continue
            for fb in range(9): r[idx[(x,ww,a,d,fb,fg)]]+=1
        eq_rows.append(r); eq_rhs.append(PQ[(x,0,z,ww)][a,:,c,d].sum())
    def mcols(x,y,z,ww,a,b,c,d):
        return [idx[(x,ww,a,d,fb,fg)] for fb in range(9) if fB(fb,y)==b for fg in range(9) if fC(fg,z)==c]
    contexts=[('x',)+t for t in product(range(2),repeat=3)]+[('w',)+t for t in product(range(2),repeat=3)]
    total=N; ctx_slack={}
    for ci,_ in enumerate(contexts):
        ctx_slack[ci]=list(range(total,total+27)); total+=27
    delta=total; total+=1
    Aub=[]; bub=[]
    for ci,ctx in enumerate(contexts):
        for oi,out in enumerate(product(range(3),repeat=3)):
            uv=ctx_slack[ci][oi]; r1=np.zeros(total); r2=np.zeros(total)
            if ctx[0]=='x':
                _,ww,y,z=ctx; b,c,d=out
                for a in range(3):
                    for col in mcols(0,y,z,ww,a,b,c,d): r1[col]+=1; r2[col]-=1
                    for col in mcols(1,y,z,ww,a,b,c,d): r1[col]-=1; r2[col]+=1
            else:
                _,x,y,z=ctx; a,b,c=out
                for d in range(3):
                    for col in mcols(x,y,z,0,a,b,c,d): r1[col]+=1; r2[col]-=1
                    for col in mcols(x,y,z,1,a,b,c,d): r1[col]-=1; r2[col]+=1
            r1[uv]-=1; r2[uv]-=1; Aub+=[r1,r2]; bub+=[0,0]
        r=np.zeros(total)
        for uv in ctx_slack[ci]: r[uv]=1
        r[delta]=-2; Aub.append(r); bub.append(0)
    Aeq=np.zeros((len(eq_rows),total))
    for i,r in enumerate(eq_rows): Aeq[i,:N]=r
    cobj=np.zeros(total); cobj[delta]=1
    res=linprog(cobj,A_ub=np.array(Aub),b_ub=np.array(bub),A_eq=Aeq,b_eq=np.array(eq_rhs),
                bounds=[(0,None)]*total,method='highs')
    return res.fun if res.success else None
