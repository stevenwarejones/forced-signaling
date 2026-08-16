#!/usr/bin/env python3
"""Standalone verifier for the Theorem 2 certificates (Sigma_LC4_certificates.json).
Rebuilds the exact LC4 quantum marginals and the LP constraint system from first
principles (sympy exact arithmetic over Q(sqrt2)); does NOT import the generating code.
Checks: primal feasibility (nonneg, 132 equalities, 272 inequalities) with objective
(sqrt2-1)/4, and dual feasibility with the same objective. Run: python3 verify_Sigma.py"""
import json, numpy as np, sympy as sp
from itertools import product
R2 = sp.sqrt(2)

# exact quantum marginals of |LC4> with the Li et al. measurements
def skron(*Ms):
    out = sp.Matrix([[1]])
    for M in Ms:
        out = sp.Matrix(np.kron(np.array(out.tolist(),dtype=object), np.array(M.tolist(),dtype=object)))
    return out
X=sp.Matrix([[0,1],[1,0]]); Z=sp.Matrix([[1,0],[0,-1]]); I=sp.eye(2)
plus=sp.Matrix([1,1])/R2
psi=plus
for _ in range(3):
    psi=sp.Matrix(np.kron(np.array(psi.tolist(),dtype=object),np.array(plus.tolist(),dtype=object)))
def CZ(n,i,j):
    M=sp.zeros(2**n,2**n)
    for s in range(2**n): M[s,s]=-1 if ((s>>(n-1-i))&1 and (s>>(n-1-j))&1) else 1
    return M
psi=CZ(4,0,1)*CZ(4,1,2)*CZ(4,2,3)*psi
pp=lambda O:{0:(I+O)/2,1:(I-O)/2}
PA=[pp(X),pp(Z)]; PB=[pp((Z+X)/R2),pp((Z-X)/R2)]; PC=[pp(Z),pp(X)]; PD=[pp(X),pp(Z)]
def qp(x,y,z,w,a,b,c,d):
    return sp.nsimplify(sp.expand((psi.T*skron(PA[x][a],PB[y][b],PC[z][c],PD[w][d])*psi)[0,0]),[R2])

idx={}; k=0
for x,w,a,d,fb,fg in product(range(2),range(2),range(2),range(2),range(4),range(4)):
    idx[(x,w,a,d,fb,fg)]=k; k+=1
N=k; fB=lambda fn,y:(fn>>y)&1; fC=lambda fn,z:(fn>>z)&1
eq=[]; rhs=[]
for x,w in product(range(2),range(2)):
    r=np.zeros(N,dtype=int)
    for a,d,fb,fg in product(range(2),range(2),range(4),range(4)): r[idx[(x,w,a,d,fb,fg)]]=1
    eq.append(r); rhs.append(sp.Integer(1))
for x,y,w,a,b,d in product(range(2),repeat=6):
    r=np.zeros(N,dtype=int)
    for fb in range(4):
        if fB(fb,y)!=b: continue
        for fg in range(4): r[idx[(x,w,a,d,fb,fg)]]+=1
    eq.append(r); rhs.append(sp.simplify(sum(qp(x,y,0,w,a,b,c,d) for c in range(2))))
for x,z,w,a,c,d in product(range(2),repeat=6):
    r=np.zeros(N,dtype=int)
    for fg in range(4):
        if fC(fg,z)!=c: continue
        for fb in range(4): r[idx[(x,w,a,d,fb,fg)]]+=1
    eq.append(r); rhs.append(sp.simplify(sum(qp(x,0,z,w,a,b,c,d) for b in range(2))))
contexts=[('x',)+t for t in product(range(2),repeat=3)]+[('w',)+t for t in product(range(2),repeat=3)]
total=N; cs={}
for ci,_ in enumerate(contexts): cs[ci]=list(range(total,total+8)); total+=8
delta=total; total+=1
def mcols(x,y,z,w,a,b,c,d):
    return [idx[(x,w,a,d,fb,fg)] for fb in range(4) if fB(fb,y)==b for fg in range(4) if fC(fg,z)==c]
Aub=[]
for ci,ctx in enumerate(contexts):
    for oi,out in enumerate(product(range(2),repeat=3)):
        uv=cs[ci][oi]; r1=np.zeros(total,dtype=int); r2=np.zeros(total,dtype=int)
        if ctx[0]=='x':
            _,w,y,z=ctx; b,c,d=out
            for a in range(2):
                for col in mcols(0,y,z,w,a,b,c,d): r1[col]+=1; r2[col]-=1
                for col in mcols(1,y,z,w,a,b,c,d): r1[col]-=1; r2[col]+=1
        else:
            _,x,y,z=ctx; a,b,c=out
            for d in range(2):
                for col in mcols(x,y,z,0,a,b,c,d): r1[col]+=1; r2[col]-=1
                for col in mcols(x,y,z,1,a,b,c,d): r1[col]-=1; r2[col]+=1
        r1[uv]-=1; r2[uv]-=1; Aub+=[r1,r2]
    r=np.zeros(total,dtype=int)
    for uv in cs[ci]: r[uv]=1
    r[delta]=-2; Aub.append(r)
Aub=np.array(Aub)   # 272 rows
assert Aub.shape[0]==272 and len(eq)==132

import os
cert=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'Sigma_LC4_certificates.json')))
parse=lambda s: sp.nsimplify(sp.sympify(s),[R2])
tF=[sp.Integer(0)]*total
for i,v in cert['primal_t_nonzero'].items(): tF[int(i)]=parse(v)
lam=[sp.Integer(0)]*272
for i,v in cert['dual_lambda_nonzero'].items(): lam[int(i)]=parse(v)
mu=[parse(v) for v in cert['dual_mu']]
target=(R2-1)/4

okPn=all(e>=0 for e in tF)  # full primal vector incl. slacks and Delta
okPe=all(sp.simplify(sum(sp.Integer(int(eq[i][j]))*tF[j] for j in np.nonzero(eq[i])[0])-rhs[i])==0 for i in range(132))
okPu=all(sp.simplify(sum(sp.Integer(int(Aub[i,j]))*tF[j] for j in np.nonzero(Aub[i])[0]))<=0 for i in range(272))
okPv=sp.simplify(tF[delta]-target)==0
print(f"PRIMAL: nonneg {okPn}, equalities {okPe}, inequalities {okPu}, value=(sqrt2-1)/4 {okPv}")
okDl=all(l>=0 for l in lam)
okDf=True
for j in range(total):
    s=sum(-sp.Integer(int(Aub[i,j]))*lam[i] for i in np.nonzero(Aub[:,j])[0])
    for i in range(132):
        if eq[i][j] if j<N else 0: s+=sp.Integer(int(eq[i][j]))*mu[i]
    if sp.simplify(s-(sp.Integer(1) if j==delta else sp.Integer(0)))>0: okDf=False; break
okDv=sp.simplify(sum(mu[i]*rhs[i] for i in range(132))-target)==0
print(f"DUAL: sign {okDl}, feasibility {okDf}, value=(sqrt2-1)/4 {okDv}")
ok=okPn and okPe and okPu and okPv and okDl and okDf and okDv
print("VERDICT:", "CERTIFICATES VALID: Sigma_HIC(Q_LC4) = (sqrt2-1)/4 exactly" if ok else "FAILED")
