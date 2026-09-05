import numpy as np
from adversary import solve_lp, LPFailure
from itertools import product
import sys

I2=np.eye(2); X=np.array([[0,1],[1,0]]); Z=np.array([[1,0],[0,-1]])
def kron(*ops):
    out=np.array([[1.0]])
    for o in ops: out=np.kron(out,o)
    return out
plus=np.array([1,1])/np.sqrt(2)
psi=plus
for _ in range(3): psi=np.kron(psi,plus)
def CZ(n,i,j):
    M=np.eye(2**n)
    for s in range(2**n):
        if (s>>(n-1-i))&1 and (s>>(n-1-j))&1: M[s,s]=-1
    return M
psi=CZ(4,0,1)@CZ(4,1,2)@CZ(4,2,3)@psi   # A,B,C,D

def projs(O):
    vals,vecs=np.linalg.eigh(O); P={}
    for o in (0,1):
        t=1 if o==0 else -1
        P[o]=sum(np.outer(vecs[:,k],vecs[:,k]) for k in range(2) if abs(vals[k]-t)<1e-9)
    return P

def sigma_hic(nA, nB, thetasA, phisB):
    """A: nA settings cos(th)X+sin(th)Z ; B: nB settings cos(ph)Z+sin(ph)X ; C:[Z,X]; D:[X,Z]"""
    PA=[projs(np.cos(t)*X+np.sin(t)*Z) for t in thetasA]
    PB=[projs(np.cos(p)*Z+np.sin(p)*X) for p in phisB]
    PC=[projs(Z),projs(X)]; PD=[projs(X),projs(Z)]
    PQ={}
    for x in range(nA):
        for y in range(nB):
            for z,w in product(range(2),range(2)):
                for a,b,c,d in product(range(2),repeat=4):
                    M=kron(PA[x][a],PB[y][b],PC[z][c],PD[w][d])
                    PQ[(x,y,z,w,a,b,c,d)]=np.real(psi@M@psi)
    nFB = 2**nB   # B response functions
    idx={}; k=0
    for x,w,a,d,fb,fg in product(range(nA),range(2),range(2),range(2),range(nFB),range(4)):
        idx[(x,w,a,d,fb,fg)]=k; k+=1
    N=k
    def fB(fn,y): return (fn>>y)&1
    def fC(fn,z): return (fn>>z)&1
    eq_rows=[]; eq_rhs=[]
    for x,w in product(range(nA),range(2)):
        r=np.zeros(N)
        for a,d,fb,fg in product(range(2),range(2),range(nFB),range(4)):
            r[idx[(x,w,a,d,fb,fg)]]=1
        eq_rows.append(r); eq_rhs.append(1.0)
    for x,y,w,a,b,d in product(range(nA),range(nB),range(2),range(2),range(2),range(2)):
        r=np.zeros(N)
        for fb in range(nFB):
            if fB(fb,y)!=b: continue
            for fg in range(4): r[idx[(x,w,a,d,fb,fg)]]+=1
        q=sum(PQ[(x,y,0,w,a,b,c,d)] for c in range(2))
        eq_rows.append(r); eq_rhs.append(q)
    for x,z,w,a,c,d in product(range(nA),range(2),range(2),range(2),range(2),range(2)):
        r=np.zeros(N)
        for fg in range(4):
            if fC(fg,z)!=c: continue
            for fb in range(nFB): r[idx[(x,w,a,d,fb,fg)]]+=1
        q=sum(PQ[(x,0,z,w,a,b,c,d)] for b in range(2))
        eq_rows.append(r); eq_rhs.append(q)
    # signaling: all pairs of x values, and w-flip; remote marginal = all other outputs
    contexts=[]
    for x1,x2 in product(range(nA),range(nA)):
        if x1<x2:
            for w,y,z in product(range(2),range(nB),range(2)):
                contexts.append(('x',x1,x2,w,y,z))
    for x,y,z in product(range(nA),range(nB),range(2)):
        contexts.append(('w',x,y,z))
    total=N; ctx_slack={}
    for ci,_ in enumerate(contexts):
        ctx_slack[ci]=list(range(total,total+8)); total+=8
    delta=total; total+=1
    def mcols(x,y,z,w,a,b,c,d):
        return [idx[(x,w,a,d,fb,fg)] for fb in range(nFB) if fB(fb,y)==b for fg in range(4) if fC(fg,z)==c]
    Aub=[]; bub=[]
    for ci,ctx in enumerate(contexts):
        for oi,out in enumerate(product(range(2),repeat=3)):
            uv=ctx_slack[ci][oi]; r1=np.zeros(total); r2=np.zeros(total)
            if ctx[0]=='x':
                _,x1,x2,w,y,z=ctx; b,c,d=out
                for a in range(2):
                    for col in mcols(x1,y,z,w,a,b,c,d): r1[col]+=1; r2[col]-=1
                    for col in mcols(x2,y,z,w,a,b,c,d): r1[col]-=1; r2[col]+=1
            else:
                _,x,y,z=ctx; a,b,c=out
                for d in range(2):
                    for col in mcols(x,y,z,0,a,b,c,d): r1[col]+=1; r2[col]-=1
                    for col in mcols(x,y,z,1,a,b,c,d): r1[col]-=1; r2[col]+=1
            r1[uv]-=1; r2[uv]-=1; Aub+=[r1,r2]; bub+=[0,0]
        r=np.zeros(total)
        for uv in ctx_slack[ci]: r[uv]=1
        r[delta]=-2
        Aub.append(r); bub.append(0)
    Aeq=np.zeros((len(eq_rows),total))
    for i,r in enumerate(eq_rows): Aeq[i,:N]=r
    beq=np.array(eq_rhs)
    cobj=np.zeros(total); cobj[delta]=1
    res=solve_lp(cobj,A_ub=np.array(Aub),b_ub=np.array(bub),A_eq=Aeq,b_eq=beq,
                 bounds=[(0,None)]*total,context=f"chained nA={nA},nB={nB}")
    return res.fun

if __name__=="__main__":
    ref=(np.sqrt(2)-1)/4
    # n=2 sanity: LC4 angles A:[0,pi/2] (X,Z), B:[pi/4,-pi/4]
    s2=sigma_hic(2,2,[0,np.pi/2],[np.pi/4,-np.pi/4])
    print(f"n=2 (LC4 angles): Sigma = {s2:.12f}   ref (sqrt2-1)/4 = {ref:.12f}")
    # n=3 chained angles: A_j = j*pi/3 ; B_k = (2k+1)*pi/6
    s3=sigma_hic(3,3,[0,np.pi/3,2*np.pi/3],[np.pi/6,np.pi/2,5*np.pi/6])
    print(f"n=3 (chained):    Sigma = {s3:.12f}   chained bound guesses: (3cos(pi/6)*2-4)/? ")
    # n=4 chained: A_j=j*pi/4, B_k=(2k+1)*pi/8
    s4=sigma_hic(4,4,[j*np.pi/4 for j in range(4)],[(2*k+1)*np.pi/8 for k in range(4)])
    print(f"n=4 (chained):    Sigma = {s4:.12f}")
    for n,s in [(2,s2),(3,s3),(4,s4)]:
        if s is not None:
            print(f"  n={n}: Sigma={s:.9f}, ratio to ref={s/ref:.6f}")
