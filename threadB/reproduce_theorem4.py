#!/usr/bin/env python3
"""Drivers for Theorem 4's numerical support and the conjecture's perturbation sweep.
(1) Exhaustive 24-vertex check of S <= 2 + 4q on the NS polytope (a finite proof step).
(2) 400 random quantum behaviors: marginal-constrained local distance == (S-2)/8 (seed 3).
(3) Perturbation sweep near the cluster point: perturbed states and measurement
    configurations never exceed (sqrt2-1)/4 (seeds 13/17; the larger 300/150 sweep
    reported in the manuscript was an external audit — this reproduces it in kind).
Solver: HiGHS with primal/dual feasibility tolerances set explicitly to 1e-9;
every solve is checked and its residuals recomputed from the returned solution."""
import numpy as np
from itertools import product
from scipy import sparse
from adversary import (SigmaLP, cluster4, projs, PA, PB, PC, PD,
                       solve_lp, LPFailure, dependency_report, diagnostics_report)
ceil=(np.sqrt(2)-1)/4
print(f"environment: {dependency_report()}")
X=np.array([[0,1],[1,0]],dtype=complex); Z=np.diag([1,-1]).astype(complex)
Y=np.array([[0,-1j],[1j,0]])

print("== (1) 24-vertex inequality S <= 2+4q ==")
def lv(fa,fb):
    q=np.zeros((2,2,2,2))
    for x,y in product(range(2),range(2)): q[x,y,fa[x],fb[y]]=1.
    return q
def pr(mu,nu,sg):
    q=np.zeros((2,2,2,2))
    for x,y,a,b in product(range(2),repeat=4):
        if (a^b)==((x&y)^(mu&x)^(nu&y)^sg): q[x,y,a,b]=.5
    return q
verts=[lv(fa,fb) for fa in product(range(2),repeat=2) for fb in product(range(2),repeat=2)]
verts+=[pr(*t) for t in product(range(2),repeat=3)]
pats=[(1,1,1,-1),(1,1,-1,1),(1,-1,1,1),(-1,1,1,1),(-1,-1,-1,1),(-1,-1,1,-1),(-1,1,-1,-1),(1,-1,-1,-1)]
viol=0; checked=0; minmargin=None
for q in verts:
    E=np.zeros((2,2))
    for x,y,a,b in product(range(2),repeat=4): E[x,y]+=((-1)**(a^b))*q[x,y,a,b]
    for s in pats:
        S=sum(s[2*x+y]*E[x,y] for x,y in product(range(2),range(2)))
        for x,y,a,b in product(range(2),repeat=4):
            if ((-1)**(a^b))==s[2*x+y]:
                checked+=1; m=2+4*q[x,y,a,b]-S
                if m<-1e-12: viol+=1
                minmargin=m if minmargin is None else min(minmargin,m)
print(f"triples checked: {checked}, violations: {viol}, min margin: {minmargin:.3g}")
assert viol==0

print("== (2) constrained distance == (S-2)/8 on 400 random quantum boxes (seed 3) ==")
def constrained_distance(Q):
    strat=list(product(range(2),repeat=2))
    verts=[(fB,fC) for fB in strat for fC in strat]
    nV=16; ncol=nV+16+1
    Ar=[];Ac=[];Av=[];bub=[];rc=0
    for ci,(y,z) in enumerate(product(range(2),range(2))):
        for b,c in product(range(2),range(2)):
            u=nV+ci*4+b*2+c
            sel=[vi for vi,(fB,fC) in enumerate(verts) if fB[y]==b and fC[z]==c]
            q=Q[y,z,b,c]
            Ar+=[rc]*len(sel)+[rc];Ac+=sel+[u];Av+=[1.]*len(sel)+[-1.];bub.append(q);rc+=1
            Ar+=[rc]*len(sel)+[rc];Ac+=sel+[u];Av+=[-1.]*len(sel)+[-1.];bub.append(-q);rc+=1
        us=[nV+ci*4+k for k in range(4)]
        Ar+=[rc]*len(us)+[rc];Ac+=us+[ncol-1];Av+=[1.]*len(us)+[-2.];bub.append(0.);rc+=1
    A=sparse.csr_matrix((Av,(Ar,Ac)),shape=(rc,ncol))
    eqr=[];eqb=[]
    r=np.zeros(ncol); r[:nV]=1; eqr.append(r); eqb.append(1.)
    for y,b in product(range(2),range(2)):
        r=np.zeros(ncol)
        for vi,(fB,fC) in enumerate(verts):
            if fB[y]==b: r[vi]=1
        eqr.append(r); eqb.append(sum(Q[y,0,b,c] for c in range(2)))
    for z,c in product(range(2),range(2)):
        r=np.zeros(ncol)
        for vi,(fB,fC) in enumerate(verts):
            if fC[z]==c: r[vi]=1
        eqr.append(r); eqb.append(sum(Q[0,z,b,c] for b in range(2)))
    cv=np.zeros(ncol); cv[-1]=1
    res=solve_lp(cv,A_ub=A,b_ub=np.array(bub),A_eq=np.array(eqr),b_eq=np.array(eqb),
                 bounds=[(0,None)]*ncol,context="constrained_distance")
    return res.fun
rng=np.random.default_rng(3)
def rq(r):
    n=r.normal(size=3); n/=np.linalg.norm(n)
    return projs(n[0]*X+n[1]*Y+n[2]*Z)
bad=0; maxdev=0.0
for i in range(400):
    v=rng.normal(size=4)+1j*rng.normal(size=4); v/=np.linalg.norm(v)
    PBq=[rq(rng) for _ in range(2)]; PCq=[rq(rng) for _ in range(2)]
    Q=np.zeros((2,2,2,2))
    for y,z,b,c in product(range(2),repeat=4):
        Q[y,z,b,c]=np.real(v.conj()@np.kron(PBq[y][b],PCq[z][c])@v)
    E=np.zeros((2,2))
    for y,z,b,c in product(range(2),repeat=4): E[y,z]+=((-1)**(b^c))*Q[y,z,b,c]
    S=max(abs(s[0]*E[0,0]+s[1]*E[0,1]+s[2]*E[1,0]+s[3]*E[1,1]) for s in pats[:4])
    d=constrained_distance(Q)
    dev=abs(d-max(0,(S-2)/8)); maxdev=max(maxdev,dev)
    if dev>1e-9: bad+=1
print(f"mismatches (threshold 1e-9): {bad}/400; max observed deviation {maxdev:.3e}")
assert bad==0

print("== (3) perturbation sweep near the cluster point (seeds 13/17) ==")
lp22=SigmaLP(2,2,2,PA,PB,PC,PD)
r13=np.random.default_rng(13); worst=0
for i in range(100):
    psi=cluster4()+0.02*(r13.normal(size=16)+1j*r13.normal(size=16))
    psi/=np.linalg.norm(psi)
    worst=max(worst, lp22.solve(psi,context=f"perturbed state #{i} (seed 13)"))
print(f"100 perturbed states (2% noise): max Sigma = {worst:.9f}")
assert worst<ceil+1e-9
r17=np.random.default_rng(17); worst2=0
ang=lambda t: projs(np.cos(t)*Z+np.sin(t)*X)
angA=lambda t: projs(np.cos(t)*X+np.sin(t)*Z)
for i in range(50):
    e=lambda: r17.normal(scale=0.05)
    lp=SigmaLP(2,2,2,[angA(e()),angA(np.pi/2+e())],[ang(np.pi/4+e()),ang(-np.pi/4+e())],
               [ang(e()),ang(np.pi/2+e())],[angA(e()),angA(np.pi/2+e())])
    worst2=max(worst2, lp.solve(cluster4(),context=f"perturbed config #{i} (seed 17)"))
print(f"50 perturbed measurement configs: max Sigma = {worst2:.9f}")
assert worst2<ceil+1e-9
print("== numerical diagnostics ==")
print(f"  {diagnostics_report()}")
print("all Theorem-4 / perturbation checks pass")
