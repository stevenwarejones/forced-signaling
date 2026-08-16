#!/usr/bin/env python3
"""Independent minimal verifier for the K=8 certificate (Theorem 1).
Rebuilds the LP constraint matrix from first principles (integers only),
loads K8_certificate.json, and verifies dual feasibility in exact fraction
arithmetic. Depends only on the Python standard library + numpy for indexing.
Run: python3 verify_K8.py"""
import json, numpy as np
from fractions import Fraction
from itertools import product

idx={}; n=0
for x,w,a,d,fb,fg in product(range(2),range(2),range(2),range(2),range(4),range(4)):
    idx[(x,w,a,d,fb,fg)]=n; n+=1
N=n
f=lambda fn,y:(fn>>y)&1
def P_cols(x,y,z,w,a,b,c,d):
    return [idx[(x,w,a,d,fb,fg)] for fb in range(4) if f(fb,y)==b for fg in range(4) if f(fg,z)==c]
# objective: the K=8-optimal completion ((0,1),(0,1),0,0,(1,0),0)
obj=np.zeros(N,dtype=int)
def add(coef,x,y,z,w,sgn):
    for a,b,c,d in product(range(2),repeat=4):
        for col in P_cols(x,y,z,w,a,b,c,d): obj[col]+=coef*sgn(a,b,c,d)
add(1,0,0,0,1,lambda a,b,c,d:(-1)**(a^b)); add(1,0,1,0,1,lambda a,b,c,d:(-1)**(a^b))
add(1,1,0,0,0,lambda a,b,c,d:(-1)**(a^b^d)); add(-1,1,1,0,0,lambda a,b,c,d:(-1)**(a^b^d))
add(2,1,0,0,0,lambda a,b,c,d:(-1)**(c^d)); add(2,0,0,1,1,lambda a,b,c,d:(-1)**(a^c^d))
A_eq=np.zeros((4,N),dtype=int)
for i,(x,w) in enumerate(product(range(2),range(2))):
    for a,d,fb,fg in product(range(2),range(2),range(4),range(4)):
        A_eq[i, idx[(x,w,a,d,fb,fg)]]=1
contexts=[('x',)+t for t in product(range(2),repeat=3)]+[('w',)+t for t in product(range(2),repeat=3)]
total=N; ctx_slack={}
for ci,_ in enumerate(contexts):
    ctx_slack[ci]=list(range(total,total+8)); total+=8
rows=[]
for ci,ctx in enumerate(contexts):
    for oi,out in enumerate(product(range(2),repeat=3)):
        u=ctx_slack[ci][oi]; r1=np.zeros(total,dtype=int); r2=np.zeros(total,dtype=int)
        if ctx[0]=='x':
            _,w,y,z=ctx; b,c,d=out
            for a in range(2):
                for col in P_cols(0,y,z,w,a,b,c,d): r1[col]+=1; r2[col]-=1
                for col in P_cols(1,y,z,w,a,b,c,d): r1[col]-=1; r2[col]+=1
        else:
            _,x,y,z=ctx; a,b,c=out
            for d in range(2):
                for col in P_cols(x,y,z,0,a,b,c,d): r1[col]+=1; r2[col]-=1
                for col in P_cols(x,y,z,1,a,b,c,d): r1[col]-=1; r2[col]+=1
        r1[u]-=1; r2[u]-=1; rows.append(r1); rows.append(r2)
tv_start=len(rows)
for ci,_ in enumerate(contexts):
    r=np.zeros(total,dtype=int)
    for u in ctx_slack[ci]: r[u]=1
    rows.append(r)
Aub=np.array(rows)
import os
cert=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'K8_certificate.json')))
y=[Fraction(0)]*Aub.shape[0]
for i,v in cert['nonzero_dual_entries'].items(): y[int(i)]=Fraction(v)
mu=[Fraction(v) for v in cert['normalization_duals']]
assert all(v>=0 for v in y), "dual sign violation"
ok=True; minslack=None
for j in range(total):
    s=sum(Fraction(int(Aub[i,j]))*y[i] for i in np.nonzero(Aub[:,j])[0])
    if j<N: s+=sum(Fraction(int(A_eq[i,j]))*mu[i] for i in range(4) if A_eq[i,j])
    target=Fraction(int(obj[j])) if j<N else Fraction(0)
    sl=s-target
    if sl<0: ok=False
    minslack=sl if minslack is None else min(minslack,sl)
tvsum=sum(y[tv_start+i] for i in range(len(contexts)))
musum=sum(mu)
print("dual feasible:",ok,"| min slack:",minslack)
print("sum TV duals:",tvsum,"(need 4) | sum normalization duals:",musum,"(need 6)")
print("VERDICT:", "CERTIFICATE VALID: S4^op <= 6 + 8*Delta_sig for all Delta_sig >= 0"
      if ok and tvsum==4 and musum==6 else "FAILED")
