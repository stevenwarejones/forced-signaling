#!/usr/bin/env python3
"""Drivers for the Table-1 rows not covered by reproduce_core.py:
setting supersets, random states + random measurements, LC5 (FULL),
locked mixed-flavor (FULL), recipient-subset pinning for Lemma 2 (FULL). Seeds fixed. Solver: HiGHS with primal/dual
feasibility tolerances set explicitly to 1e-9; every solve is checked and its
residuals recomputed from the returned solution (reported at the end)."""
import os, numpy as np
from itertools import product
from adversary import (SigmaLP, cluster4, projs, kron, PA, PB, PC, PD,
                       solve_lp, LPFailure, dependency_report, diagnostics_report)
FULL=os.environ.get("FULL")=="1"
print(f"environment: {dependency_report()}")
ceil=(np.sqrt(2)-1)/4
X=np.array([[0,1],[1,0]],dtype=complex); Z=np.diag([1,-1]).astype(complex)
Y=np.array([[0,-1j],[1j,0]])
ang=lambda t: projs(np.cos(t)*Z+np.sin(t)*X)
angA=lambda t: projs(np.cos(t)*X+np.sin(t)*Z)
psi=cluster4()

print("== setting supersets (claim: saturate ceiling to 9 digits) ==")
tests=[("B superset",SigmaLP(2,3,2,PA,[ang(np.pi/4),ang(-np.pi/4),ang(0)],PC,PD)),
       ("A superset",SigmaLP(3,2,2,[angA(0),angA(np.pi/2),angA(np.pi/4)],PB,PC,PD)),
       ("C superset",SigmaLP(2,2,3,PA,PB,[projs(Z),projs(X),projs((Z+X)/np.sqrt(2))],PD)),
       ("A+B superset",SigmaLP(3,3,2,[angA(0),angA(np.pi/2),angA(np.pi/4)],
                               [ang(np.pi/4),ang(-np.pi/4),ang(np.pi/2)],PC,PD))]
for name,lp in tests:
    s=lp.solve(psi,context=name); print(f"  {name}: {s:.9f}")
    assert abs(s-ceil)<1e-8, name

print("== random states + random measurements (claim <1e-9; seed 11) ==")
rng=np.random.default_rng(11)
worst=0
for i in range(5 if not FULL else 15):
    v=rng.normal(size=16)+1j*rng.normal(size=16); v/=np.linalg.norm(v)
    def rq(r):
        n=r.normal(size=3); n/=np.linalg.norm(n)
        return projs(n[0]*X+n[1]*Y+n[2]*Z)
    lp=SigmaLP(2,2,2,[rq(rng) for _ in range(2)],[rq(rng) for _ in range(2)],
               [rq(rng) for _ in range(2)],[rq(rng) for _ in range(2)])
    worst=max(worst, lp.solve(v,context=f"random state+measurements #{i} (seed 11)"))
print(f"  worst: {worst:.2e}")
assert worst<1e-9

if FULL:
    print("== LC5 (claim: saturates ceiling) ==")
    plus=np.array([1,1])/np.sqrt(2); p5=plus
    for _ in range(4): p5=np.kron(p5,plus)
    p5=p5.astype(complex)
    def CZ(n,i,j):
        M=np.eye(2**n)
        for s in range(2**n):
            if (s>>(n-1-i))&1 and (s>>(n-1-j))&1: M[s,s]=-1
        return M
    p5=CZ(5,0,1)@CZ(5,1,2)@CZ(5,2,3)@CZ(5,3,4)@p5
    PE5=[projs(Z),projs(X)]
    PQ={}
    for st in product(range(2),repeat=5):
        x,y,z,w,u=st
        for o in product(range(2),repeat=5):
            a,b,c,d,e=o
            M=kron(PA[x][a],PB[y][b],PC[z][c],PD[w][d],PE5[u][e])
            PQ[st+o]=np.real(p5.conj()@M@p5)
    idx={}; k=0
    for x,w,u,a,d,e,fb,fg in product(range(2),range(2),range(2),range(2),range(2),range(2),range(4),range(4)):
        idx[(x,w,u,a,d,e,fb,fg)]=k; k+=1
    N=k; fB=lambda fn,y:(fn>>y)&1; fC=lambda fn,z:(fn>>z)&1
    eq=[]; rhs=[]
    for x,w,u in product(range(2),repeat=3):
        r=np.zeros(N)
        for a,d,e,fb,fg in product(range(2),range(2),range(2),range(4),range(4)):
            r[idx[(x,w,u,a,d,e,fb,fg)]]=1
        eq.append(r); rhs.append(1.0)
    for x,y,w,u,a,b,d,e in product(range(2),repeat=8):
        r=np.zeros(N)
        for fb in range(4):
            if fB(fb,y)!=b: continue
            for fg in range(4): r[idx[(x,w,u,a,d,e,fb,fg)]]+=1
        eq.append(r); rhs.append(sum(PQ[(x,y,0,w,u,a,b,c,d,e)] for c in range(2)))
    for x,z,w,u,a,c,d,e in product(range(2),repeat=8):
        r=np.zeros(N)
        for fg in range(4):
            if fC(fg,z)!=c: continue
            for fb in range(4): r[idx[(x,w,u,a,d,e,fb,fg)]]+=1
        eq.append(r); rhs.append(sum(PQ[(x,0,z,w,u,a,b,c,d,e)] for b in range(2)))
    def mc(x,y,z,w,u,a,b,c,d,e):
        return [idx[(x,w,u,a,d,e,fb,fg)] for fb in range(4) if fB(fb,y)==b for fg in range(4) if fC(fg,z)==c]
    ctxs=[]
    for chg in 'xwu':
        for rest in product(range(2),repeat=4): ctxs.append((chg,)+rest)
    total=N; sl={}
    for ci,_ in enumerate(ctxs): sl[ci]=list(range(total,total+16)); total+=16
    delta=total; total+=1
    Aub=[]
    for ci,ctx in enumerate(ctxs):
        chg=ctx[0]
        for oi,out in enumerate(product(range(2),repeat=4)):
            uv=sl[ci][oi]; r1=np.zeros(total); r2=np.zeros(total)
            if chg=='x':
                w_,u_,y,z=ctx[1:]; b,c,d,e=out
                for a in range(2):
                    for col in mc(0,y,z,w_,u_,a,b,c,d,e): r1[col]+=1; r2[col]-=1
                    for col in mc(1,y,z,w_,u_,a,b,c,d,e): r1[col]-=1; r2[col]+=1
            elif chg=='w':
                x_,u_,y,z=ctx[1:]; a,b,c,e=out
                for d in range(2):
                    for col in mc(x_,y,z,0,u_,a,b,c,d,e): r1[col]+=1; r2[col]-=1
                    for col in mc(x_,y,z,1,u_,a,b,c,d,e): r1[col]-=1; r2[col]+=1
            else:
                x_,w_,y,z=ctx[1:]; a,b,c,d=out
                for e in range(2):
                    for col in mc(x_,y,z,w_,0,a,b,c,d,e): r1[col]+=1; r2[col]-=1
                    for col in mc(x_,y,z,w_,1,a,b,c,d,e): r1[col]-=1; r2[col]+=1
            r1[uv]-=1; r2[uv]-=1; Aub+=[r1,r2]
        r=np.zeros(total)
        for uv in sl[ci]: r[uv]=1
        r[delta]=-2; Aub.append(r)
    Aeq=np.zeros((len(eq),total))
    for i,r in enumerate(eq): Aeq[i,:N]=r
    c=np.zeros(total); c[delta]=1
    res=solve_lp(c,A_ub=np.array(Aub),b_ub=np.zeros(len(Aub)),A_eq=Aeq,b_eq=np.array(rhs),
                 bounds=[(0,None)]*total,context="LC5")
    print(f"  LC5: {res.fun:.9f}")
    assert abs(res.fun-ceil)<1e-8

    print("== locked mixed-flavor (claim: saturates ceiling to 9 digits) ==")
    PA3=[projs(X),projs(Z),projs(Y)]
    PB4=[projs((Z+X)/np.sqrt(2)),projs((Z-X)/np.sqrt(2)),projs((Z+Y)/np.sqrt(2)),projs((Z-Y)/np.sqrt(2))]
    PC3=[projs(Z),projs(X),projs(Y)]
    PD3=[projs(X),projs(Z),projs(Y)]
    nA,nB,nC,nD=3,4,3,3
    PQ2={}
    for x,y,z,w in product(range(nA),range(nB),range(nC),range(nD)):
        for a,b,c,d in product(range(2),repeat=4):
            M=kron(PA3[x][a],PB4[y][b],PC3[z][c],PD3[w][d])
            PQ2[(x,y,z,w,a,b,c,d)]=np.real(psi.conj()@M@psi)
    nFB=2**nB; nFC=2**nC
    idx2={}; k=0
    for x,w,a,d,fb,fg in product(range(nA),range(nD),range(2),range(2),range(nFB),range(nFC)):
        idx2[(x,w,a,d,fb,fg)]=k; k+=1
    N2=k
    eq=[]; rhs=[]
    for x,w in product(range(nA),range(nD)):
        r=np.zeros(N2)
        for a,d,fb,fg in product(range(2),range(2),range(nFB),range(nFC)): r[idx2[(x,w,a,d,fb,fg)]]=1
        eq.append(r); rhs.append(1.0)
    for x,y,w,a,b,d in product(range(nA),range(nB),range(nD),range(2),range(2),range(2)):
        r=np.zeros(N2)
        for fb in range(nFB):
            if fB(fb,y)!=b: continue
            for fg in range(nFC): r[idx2[(x,w,a,d,fb,fg)]]+=1
        eq.append(r); rhs.append(sum(PQ2[(x,y,0,w,a,b,c,d)] for c in range(2)))
    for x,z,w,a,c,d in product(range(nA),range(nC),range(nD),range(2),range(2),range(2)):
        r=np.zeros(N2)
        for fg in range(nFC):
            if fC(fg,z)!=c: continue
            for fb in range(nFB): r[idx2[(x,w,a,d,fb,fg)]]+=1
        eq.append(r); rhs.append(sum(PQ2[(x,0,z,w,a,b,c,d)] for b in range(2)))
    def mc2(x,y,z,w,a,b,c,d):
        return [idx2[(x,w,a,d,fb,fg)] for fb in range(nFB) if fB(fb,y)==b for fg in range(nFC) if fC(fg,z)==c]
    ctxs=[]
    for x1 in range(nA):
        for x2 in range(x1+1,nA):
            for w,y,z in product(range(nD),range(nB),range(nC)): ctxs.append(('x',x1,x2,w,y,z))
    for w1 in range(nD):
        for w2 in range(w1+1,nD):
            for x,y,z in product(range(nA),range(nB),range(nC)): ctxs.append(('w',w1,w2,x,y,z))
    total=N2; sl={}
    for ci,_ in enumerate(ctxs): sl[ci]=list(range(total,total+8)); total+=8
    delta=total; total+=1
    Aub=[]
    for ci,ctx in enumerate(ctxs):
        for oi,out in enumerate(product(range(2),repeat=3)):
            uv=sl[ci][oi]; r1=np.zeros(total); r2=np.zeros(total)
            if ctx[0]=='x':
                _,x1,x2,w,y,z=ctx; b,c,d=out
                for a in range(2):
                    for col in mc2(x1,y,z,w,a,b,c,d): r1[col]+=1; r2[col]-=1
                    for col in mc2(x2,y,z,w,a,b,c,d): r1[col]-=1; r2[col]+=1
            else:
                _,w1,w2,x,y,z=ctx; a,b,c=out
                for d in range(2):
                    for col in mc2(x,y,z,w1,a,b,c,d): r1[col]+=1; r2[col]-=1
                    for col in mc2(x,y,z,w2,a,b,c,d): r1[col]-=1; r2[col]+=1
            r1[uv]-=1; r2[uv]-=1; Aub+=[r1,r2]
        r=np.zeros(total)
        for uv in sl[ci]: r[uv]=1
        r[delta]=-2; Aub.append(r)
    Aeq=np.zeros((len(eq),total))
    for i,r in enumerate(eq): Aeq[i,:N2]=r
    c=np.zeros(total); c[delta]=1
    res=solve_lp(c,A_ub=np.array(Aub),b_ub=np.zeros(len(Aub)),A_eq=Aeq,b_eq=np.array(rhs),
                 bounds=[(0,None)]*total,context="locked mixed-flavor")
    print(f"  locked mixed-flavor: {res.fun:.9f}")
    assert abs(res.fun-ceil)<1e-8

if FULL:
    print("== recipient-subset pinning (Lemma 2 of the manuscript) ==")
    # Analytic claim: any recipient set omitting B or omitting C has its distribution
    # fixed by the reproduced ABD/ACD families, which are no-signaling, so it cannot
    # depend on the sender's setting. Only {B,C} and {B,C,D} are live for sender A.
    # This maximises the signal into each subset by LP and confirms exactly that.
    from itertools import combinations
    lpP = SigmaLP(2,2,2,PA,PB,PC,PD)
    idxP, fBp, fCp = lpP.idx, lpP.fB, lpP.fC
    def colsP(x,y,z,w,a,b,c,d):
        return [idxP[(x,w,a,d,fb,fg)] for fb in range(4) if fBp(fb,y)==b
                for fg in range(4) if fCp(fg,z)==c]
    beqP = lpP.marginals(psi)
    names = {1:'B', 2:'C', 3:'D'}
    live, pinned = {}, {}
    for size in (1,2,3):
        for R in combinations((1,2,3), size):
            best = 0.0
            for outs in product(range(2), repeat=size):
                for y,z,w in product(range(2), repeat=3):
                    obj = np.zeros(lpP.total)
                    for a,b,c,d in product(range(2), repeat=4):
                        vals = {0:a, 1:b, 2:c, 3:d}
                        if tuple(vals[i] for i in R) != outs: continue
                        for col in colsP(0,y,z,w,a,b,c,d): obj[col] += 1
                        for col in colsP(1,y,z,w,a,b,c,d): obj[col] -= 1
                    r = solve_lp(-obj, A_ub=lpP.Aub, b_ub=lpP.bub, A_eq=lpP.Aeq, b_eq=beqP,
                                 bounds=[(0,None)]*lpP.total,
                                 context=f"max signal A -> {R}")
                    best = max(best, -r.fun)
            label = ''.join(names[i] for i in R)
            (live if best > 1e-7 else pinned)[label] = best
            print(f"  sender A -> {label:4s}: max component {best:.9f}"
                  f"  {'CAN signal' if best > 1e-7 else 'pinned'}")
    assert set(pinned) == {'B','C','D','BD','CD'}, f"unexpected pinned set: {sorted(pinned)}"
    assert set(live) == {'BC','BCD'}, f"unexpected live set: {sorted(live)}"
    assert max(pinned.values()) < 1e-7, f"a pinned subset signaled: {pinned}"
    print("  only {B,C} and {B,C,D} are live, as Lemma 2 states")

print("== numerical diagnostics ==")
print(f"  {diagnostics_report()}")
print(f"  mode: {'FULL' if FULL else 'QUICK'}")
print("reproduce_extra.py: all assertions passed")
