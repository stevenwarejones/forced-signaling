#!/usr/bin/env python3
"""Reproduce the core floating-point numerical claims of the manuscript.
Modes: QUICK (default) validates a documented subset; FULL=1 runs everything.
Solver: scipy linprog (HiGHS), feasibility tolerance ~1e-9. Seeds fixed below.
Claims covered here: LC4 numeric, tilt identity, GHZ/W/random states, chained
n=3(,4 FULL), single-copy Tsirelson distance, parallel flatness (FULL),
completion spectrum (1/8 subset; FULL for complete), Monte-Carlo cover, budgets.
See reproduce_extra.py for: setting supersets, random states + random
measurements, LC5 (FULL), locked mixed-flavor (FULL).
NOT covered (open items, see MANIFEST): CGLMP ladder script, qutrit/Barnea evaluations."""
import os, numpy as np
from itertools import product
from scipy.optimize import linprog
from scipy import sparse
from adversary import SigmaLP, cluster4, PA, PB, PC, PD, projs, kron
FULL = os.environ.get("FULL") == "1"
ceil = (np.sqrt(2)-1)/4
rng = np.random.default_rng(7)
lp22 = SigmaLP(2,2,2,PA,PB,PC,PD)

print("== Table 1 rows (Sigma_HIC via LP) ==")
s_lc4=lp22.solve(cluster4())
print(f"LC4 cluster:            {s_lc4:.9f}  (claim: saturates {ceil:.9f})")
assert abs(s_lc4-ceil)<1e-8
for th in (0.85,0.90,0.95):
    s = lp22.solve(cluster4(np.pi*th))
    print(f"tilted theta={th}: Sigma={s:.9f}")
ghz=np.zeros(16,dtype=complex); ghz[0]=ghz[15]=1/np.sqrt(2)
w4=np.zeros(16,dtype=complex)
for i in (1,2,4,8): w4[i]=0.5
print(f"GHZ4: {lp22.solve(ghz):.2e}   W4: {lp22.solve(w4):.2e}  (claim: <1e-9)")
worst=0.0
for _ in range(10 if not FULL else 40):
    v=rng.normal(size=16)+1j*rng.normal(size=16); v/=np.linalg.norm(v)
    worst=max(worst, lp22.solve(v) or 0)
print(f"random states, worst: {worst:.2e}  (claim: <1e-9; seed 7)")

print("== chained values (imports chained.sigma_hic) ==")
from chained import sigma_hic
print(f"n=3: {sigma_hic(3,3,[0,np.pi/3,2*np.pi/3],[np.pi/6,np.pi/2,5*np.pi/6]):.9f}  (claim 0.0748)")
if FULL:
    print(f"n=4: {sigma_hic(4,4,[j*np.pi/4 for j in range(4)],[(2*k+1)*np.pi/8 for k in range(4)]):.9f}  (claim 0.0766)")

print("== parallel flatness / CGLMP (max-context TV distance to local) ==")
def local_distance(nset,nout,Q):
    stratB=list(product(range(nout),repeat=nset)); stratC=stratB
    nB=len(stratB); nV=nB*nB
    Bm=np.array(stratB)
    Ar=[];Ac=[];Av=[];bub=[];rc=0; nU=nset*nset*nout*nout; ncol=nV+nU+1
    for ci,(y,z) in enumerate(product(range(nset),range(nset))):
        bofv=np.repeat(Bm[:,y],nB); cofv=np.tile(Bm[:,z],nB)
        for b,c in product(range(nout),range(nout)):
            u=nV+ci*nout*nout+b*nout+c
            sel=np.where((bofv==b)&(cofv==c))[0]; q=Q[y,z,b,c]
            Ar+=[rc]*len(sel)+[rc];Ac+=list(sel)+[u];Av+=[1.]*len(sel)+[-1.];bub.append(q);rc+=1
            Ar+=[rc]*len(sel)+[rc];Ac+=list(sel)+[u];Av+=[-1.]*len(sel)+[-1.];bub.append(-q);rc+=1
        us=[nV+ci*nout*nout+k for k in range(nout*nout)]
        Ar+=[rc]*len(us)+[rc];Ac+=us+[ncol-1];Av+=[1.]*len(us)+[-2.];bub.append(0.);rc+=1
    A=sparse.csr_matrix((Av,(Ar,Ac)),shape=(rc,ncol))
    Ae=sparse.csr_matrix((np.ones(nV),(np.zeros(nV),np.arange(nV))),shape=(1,ncol))
    c=np.zeros(ncol); c[-1]=1
    r=linprog(c,A_ub=A,b_ub=np.array(bub),A_eq=Ae,b_eq=[1.],bounds=[(0,None)]*ncol,method='highs')
    return r.fun
r2v=1/np.sqrt(2); E={(0,0):r2v,(0,1):r2v,(1,0):r2v,(1,1):-r2v}
Q1=np.zeros((2,2,2,2))
for y,z,b,c in product(range(2),repeat=4): Q1[y,z,b,c]=(1+((-1)**(b^c))*E[(y,z)])/4
d1=local_distance(2,2,Q1); print(f"single Tsirelson: {d1:.12f}")
assert abs(d1-ceil)<1e-9
if FULL:
    Q2=np.zeros((4,4,4,4))
    for y,z,b,c in product(range(4),repeat=4):
        Q2[y,z,b,c]=Q1[y//2,z//2,b//2,c//2]*Q1[y%2,z%2,b%2,c%2]
    print(f"parallel two-copy: {local_distance(4,4,Q2):.12f}  (claim: equals single to 12 digits)")

print("== completion spectrum (Theorem 1 context) ==")
idx={}; n=0
for x,w,a,d,fb,fg in product(range(2),range(2),range(2),range(2),range(4),range(4)):
    idx[(x,w,a,d,fb,fg)]=n; n+=1
N=n; f=lambda fn,y:(fn>>y)&1
def cols(x,y,z,w,a,b,c,d):
    return [idx[(x,w,a,d,fb,fg)] for fb in range(4) if f(fb,y)==b for fg in range(4) if f(fg,z)==c]
Aeq=np.zeros((4,N))
for i,(x,w) in enumerate(product(range(2),range(2))):
    for a,d,fb,fg in product(range(2),range(2),range(4),range(4)): Aeq[i,idx[(x,w,a,d,fb,fg)]]=1
ctxs=[('x',)+u for u in product(range(2),repeat=3)]+[('w',)+u for u in product(range(2),repeat=3)]
total=N; sl={}
for ci,_ in enumerate(ctxs): sl[ci]=list(range(total,total+8)); total+=8
rows=[]
for ci,ctx in enumerate(ctxs):
    for oi,out in enumerate(product(range(2),repeat=3)):
        u=sl[ci][oi]; r1=np.zeros(total); r2=np.zeros(total)
        if ctx[0]=='x':
            _,w,y,z=ctx; b,c,d=out
            for a in range(2):
                for cc in cols(0,y,z,w,a,b,c,d): r1[cc]+=1; r2[cc]-=1
                for cc in cols(1,y,z,w,a,b,c,d): r1[cc]-=1; r2[cc]+=1
        else:
            _,x,y,z=ctx; a,b,c=out
            for d in range(2):
                for cc in cols(x,y,z,0,a,b,c,d): r1[cc]+=1; r2[cc]-=1
                for cc in cols(x,y,z,1,a,b,c,d): r1[cc]-=1; r2[cc]+=1
        r1[u]-=1; r2[u]-=1; rows+=[r1,r2]
tv=[]
for ci,_ in enumerate(ctxs):
    r=np.zeros(total)
    for u in sl[ci]: r[u]=1
    tv.append(r)
Aub=np.vstack([rows,tv]); AeqT=np.hstack([Aeq,np.zeros((4,total-N))]); eps=0.01
def K_of(comp):
    (z1,w1),(z2,w2),z3,z4,(x5,y5),y6=comp
    obj=np.zeros(total)
    def add(cf,x,y,z,w,sg):
        for a,b,c,d in product(range(2),repeat=4):
            for cc in cols(x,y,z,w,a,b,c,d): obj[cc]+=cf*sg(a,b,c,d)
    add(1,0,0,z1,w1,lambda a,b,c,d:(-1)**(a^b)); add(1,0,1,z2,w2,lambda a,b,c,d:(-1)**(a^b))
    add(1,1,0,z3,0,lambda a,b,c,d:(-1)**(a^b^d)); add(-1,1,1,z4,0,lambda a,b,c,d:(-1)**(a^b^d))
    add(2,x5,y5,0,0,lambda a,b,c,d:(-1)**(c^d)); add(2,0,y6,1,1,lambda a,b,c,d:(-1)**(a^c^d))
    bub=np.concatenate([np.zeros(len(rows)),2*eps*np.ones(len(tv))])
    r=linprog(-obj,A_ub=Aub,b_ub=bub,A_eq=AeqT,b_eq=np.ones(4),bounds=[(0,None)]*total,method='highs')
    return round((-r.fun-6)/eps,6)
comps=list(product(product(range(2),range(2)),product(range(2),range(2)),range(2),range(2),product(range(2),range(2)),range(2)))
subset = comps if FULL else comps[::8]  # documented QUICK subset: every 8th completion
from collections import Counter
cnt=Counter(K_of(c) for c in subset)
print(f"completions evaluated: {len(subset)} ({'FULL' if FULL else 'QUICK 1/8 subset'}); spectrum: {dict(sorted(cnt.items()))}")
print("(FULL claim: {8:64, 10:128, 12:128, 14:128, 16:64})")

print("== Monte-Carlo delay cover (seed 1) ==")
B,Kv=1.34e-3,1e4
gB=1/np.sqrt(1-B*B); a=1/(gB*np.sqrt(Kv*Kv-1)); phi=np.arctan(a); U=np.arcsin(B*np.cos(phi))
Nc=int(np.floor(U/phi))+1
th=np.array([-U+(j+.5)*2*U/Nc for j in range(Nc)]); rhos=np.sin(th)
def blind(rho,q,beta,kappa): return (kappa*kappa-1)/(1-beta*beta)*(rho-q)**2 < 1-rho*rho
r1=np.random.default_rng(1); fails=0
Ntrial=500000 if FULL else 50000
for _ in range(Ntrial):
    beta=B*r1.uniform()**(1/3); kap=1+(Kv-1)*r1.uniform()**0.25; q=beta*r1.uniform(-1,1)
    if not any(blind(rr,q,beta,kap) for rr in rhos): fails+=1
print(f"N={Nc}, trials={Ntrial}, coverage failures: {fails} (claim: 0)")
assert fails==0

print("== event budgets ==")
import math
for nu in (0.95,0.90):
    m=nu*(4+2*math.sqrt(2))-6
    print(f"nu={nu}: N >= {(8*5/m)**2:,.0f} (claims ~6,746 / ~75,490)")
