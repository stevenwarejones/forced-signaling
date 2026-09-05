#!/usr/bin/env python3
"""Reproduce the core floating-point numerical claims of the manuscript.
Modes: QUICK (default) validates a documented subset; FULL=1 runs everything.
Solver: scipy linprog (HiGHS) with primal/dual feasibility tolerances set
explicitly to 1e-9 (HiGHS defaults are 1e-7); every solve is checked for success
and its residuals are recomputed from the returned solution and reported at the
end of the run.  Agreement in printed digits is observed agreement, not a
certified error bound -- see MANIFEST.md.  Seeds fixed below.
Claims covered here: LC4 numeric, tilt identity (with S_4 computed directly from
the state), GHZ/W/random states, chained n=3(,4 FULL), single-copy Tsirelson
distance, parallel flatness (FULL), completion spectrum (1/8 subset plus the
known optimal completion; FULL for complete), Monte-Carlo cover, budgets.
See reproduce_extra.py for: setting supersets, random states + random
measurements, LC5 (FULL), locked mixed-flavor (FULL).
NOT covered (open items, see MANIFEST): CGLMP ladder script, qutrit/Barnea evaluations."""
import os, numpy as np
from itertools import product
from scipy import sparse
from adversary import (SigmaLP, cluster4, PA, PB, PC, PD, projs, kron,
                       solve_lp, LPFailure, dependency_report, diagnostics_report)
FULL = os.environ.get("FULL") == "1"
ceil = (np.sqrt(2)-1)/4
rng = np.random.default_rng(7)
lp22 = SigmaLP(2,2,2,PA,PB,PC,PD)
print(f"environment: {dependency_report()}")

I2=np.eye(2); Xm=np.array([[0,1],[1,0]],dtype=complex); Zm=np.diag([1,-1]).astype(complex)
_A=[Xm,Zm]; _B=[(Zm+Xm)/np.sqrt(2),(Zm-Xm)/np.sqrt(2)]; _C=[Zm,Xm]; _D=[Xm,Zm]
def _ev(psi,*ops): return float(np.real(psi.conj()@kron(*ops)@psi))
def S4_of(psi):
    """The LC4 witness S_4 = R_4 + 2 L_4 evaluated directly from the state."""
    return (_ev(psi,_A[0],_B[0],I2,I2) + _ev(psi,_A[0],_B[1],I2,I2)
            + _ev(psi,_A[1],_B[0],I2,_D[0]) - _ev(psi,_A[1],_B[1],I2,_D[0])
            + 2*_ev(psi,I2,I2,_C[0],_D[0]) + 2*_ev(psi,_A[0],I2,_C[1],_D[1]))

print("== Table 1 rows (Sigma_HIC via LP) ==")
s_lc4=lp22.solve(cluster4(),context="LC4 cluster")
print(f"LC4 cluster:            {s_lc4:.9f}  (claim: saturates {ceil:.9f})")
assert abs(s_lc4-ceil)<1e-8
print("-- tilted family: Sigma_HIC vs max{0,(S4-6)/8}, S4 computed from the state --")
tilt_dev=0.0
for th in (0.85,0.90,0.95,1.00):
    psi_t=cluster4(np.pi*th)
    s = lp22.solve(psi_t,context=f"tilted theta={th}")
    s4 = S4_of(psi_t); pred = max(0.0,(s4-6)/8)
    tilt_dev=max(tilt_dev,abs(s-pred))
    print(f"  theta={th:.2f}: S4={s4:.12f}  Sigma={s:.12f}  max(0,(S4-6)/8)={pred:.12f}"
          f"  |diff|={abs(s-pred):.2e}")
assert tilt_dev<1e-8, f"tilted identity deviates by {tilt_dev:.3e}"
print(f"  identity holds at the 4 tested points; max deviation {tilt_dev:.2e}"
      "  (tested points only -- not a claim for all theta)")

print("-- fixed cluster-point dual applied to the tilted states' marginals --")
# The dual of the marginal-reproduction equalities is feasible independently of
# the right-hand side, so the dual obtained AT THE CLUSTER POINT is an affine
# functional of any other state's marginals, and a valid lower bound on its
# Sigma. The claim checked here is the stronger one: on this family that fixed
# functional equals the raw witness expression (S4-6)/8, including where it goes
# negative -- which is exactly why Sigma = max{0,(S4-6)/8} rather than (S4-6)/8.
mu_cluster,_ = lp22.equality_duals(cluster4(),context="cluster-point dual")
dual_dev=0.0
for th in (0.85,0.90,0.95,1.00):
    psi_t=cluster4(np.pi*th)
    val=float(mu_cluster@lp22.marginals(psi_t)); raw=(S4_of(psi_t)-6)/8
    dual_dev=max(dual_dev,abs(val-raw))
    print(f"  theta={th:.2f}: fixed dual={val:+.12f}  (S4-6)/8={raw:+.12f}"
          f"  |diff|={abs(val-raw):.2e}")
assert dual_dev<1e-8, f"fixed-dual functional deviates by {dual_dev:.3e}"
print(f"  fixed dual reproduces (S4-6)/8 at the tested points; max deviation {dual_dev:.2e}")
ghz=np.zeros(16,dtype=complex); ghz[0]=ghz[15]=1/np.sqrt(2)
w4=np.zeros(16,dtype=complex)
for i in (1,2,4,8): w4[i]=0.5
s_ghz=lp22.solve(ghz,context="GHZ4"); s_w4=lp22.solve(w4,context="W4")
print(f"GHZ4: {s_ghz:.2e}   W4: {s_w4:.2e}  (claim: <1e-9)")
assert s_ghz<1e-9 and s_w4<1e-9, f"GHZ4={s_ghz:.3e}, W4={s_w4:.3e}"
worst=0.0
for i in range(10 if not FULL else 40):
    v=rng.normal(size=16)+1j*rng.normal(size=16); v/=np.linalg.norm(v)
    worst=max(worst, lp22.solve(v,context=f"random state #{i} (seed 7)"))
print(f"random states, worst: {worst:.2e}  (claim: <1e-9; seed 7)")
assert worst<1e-9

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
    r=solve_lp(c,A_ub=A,b_ub=np.array(bub),A_eq=Ae,b_eq=np.array([1.]),
               bounds=[(0,None)]*ncol,context=f"local_distance(nset={nset},nout={nout})")
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
    d2=local_distance(4,4,Q2)
    print(f"parallel two-copy: {d2:.12f}  vs single {d1:.12f}  |diff|={abs(d2-d1):.2e}")
    print("  (observed agreement of two floating-point LPs; no exact certificate)")
    assert abs(d2-d1)<1e-9, f"parallel flatness deviates by {abs(d2-d1):.3e}"

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
    r=solve_lp(-obj,A_ub=Aub,b_ub=bub,A_eq=AeqT,b_eq=np.ones(4),
               bounds=[(0,None)]*total,context=f"completion {comp}")
    return round((-r.fun-6)/eps,6)
comps=list(product(product(range(2),range(2)),product(range(2),range(2)),range(2),range(2),product(range(2),range(2)),range(2)))
OPT_COMP=((0,1),(0,1),0,0,(1,0),0)   # the K=8 completion certified in Theorem 1
from collections import Counter
if FULL:
    subset=comps
else:
    # The plain 1/8 stride omits both K=8 and K=10, so the optimal completion --
    # the one the manuscript's sensitivity claim rests on -- is added explicitly.
    subset=comps[::8]+[OPT_COMP]
cnt=Counter(K_of(c) for c in subset)
print(f"completions evaluated: {len(subset)} "
      f"({'FULL (all 512)' if FULL else 'QUICK: 1/8 stride + the optimal completion'}); "
      f"spectrum: {dict(sorted(cnt.items()))}")
k_opt=K_of(OPT_COMP)
print(f"optimal completion {OPT_COMP}: K={k_opt} (exact optimality is Corollary 1, not this LP)")
assert k_opt==8, f"optimal completion gave K={k_opt}, expected 8"
assert max(cnt)<=16 and min(cnt)>=8, f"K outside the claimed range: {dict(cnt)}"
if FULL:
    expected={8:64,10:128,12:128,14:128,16:64}
    assert dict(sorted(cnt.items()))==expected, f"spectrum {dict(sorted(cnt.items()))} != {expected}"
    print("  FULL spectrum matches the claimed {8:64, 10:128, 12:128, 14:128, 16:64}")
else:
    print("  (QUICK mode: multiplicities are NOT exercised; run FULL=1 for the full spectrum)")

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

print("== numerical diagnostics ==")
print(f"  {diagnostics_report()}")
print(f"  mode: {'FULL' if FULL else 'QUICK'}; {dependency_report()}")
print("reproduce_core.py: all assertions passed")
