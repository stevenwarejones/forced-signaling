import numpy as np
from adversary import solve_lp, LPFailure
from itertools import product

I2=np.eye(2); X=np.array([[0,1],[1,0]]); Z=np.array([[1,0],[0,-1]])
phip=np.array([1,0,0,1],dtype=complex)/np.sqrt(2)
PSI=np.kron(phip,phip)  # qubit order A, b1, b2, C

def qubit_proj(theta):
    O=np.cos(theta)*Z+np.sin(theta)*X
    vals,vecs=np.linalg.eigh(O)
    return {0:np.outer(vecs[:,np.argmax(vals)],vecs[:,np.argmax(vals)].conj()),
            1:np.outer(vecs[:,np.argmin(vals)],vecs[:,np.argmin(vals)].conj())}

class BarneaLP:
    def __init__(self):
        idx={}; k=0
        for y,b,al,ga in product(range(3),range(3),range(4),range(4)):
            idx[(y,b,al,ga)]=k; k+=1
        self.idx=idx; self.N=k
        fA=lambda fn,x:(fn>>x)&1; fC=lambda fn,z:(fn>>z)&1
        eq=[]; self.specs=[]
        for y in range(3):
            r=np.zeros(self.N)
            for b,al,ga in product(range(3),range(4),range(4)): r[idx[(y,b,al,ga)]]=1
            eq.append(r); self.specs.append(('norm',))
        for x,y,a,b in product(range(2),range(3),range(2),range(3)):
            r=np.zeros(self.N)
            for al in range(4):
                if fA(al,x)!=a: continue
                for ga in range(4): r[idx[(y,b,al,ga)]]+=1
            eq.append(r); self.specs.append(('AB',x,y,a,b))
        for y,z,b,c in product(range(3),range(2),range(3),range(2)):
            r=np.zeros(self.N)
            for ga in range(4):
                if fC(ga,z)!=c: continue
                for al in range(4): r[idx[(y,b,al,ga)]]+=1
            eq.append(r); self.specs.append(('BC',y,z,b,c))
        contexts=[]
        for y1 in range(3):
            for y2 in range(y1+1,3):
                for x,z in product(range(2),range(2)):
                    contexts.append((y1,y2,x,z))
        total=self.N; slack={}
        for ci,_ in enumerate(contexts):
            slack[ci]=list(range(total,total+4)); total+=4
        delta=total; total+=1
        Aub=[]
        for ci,(y1,y2,x,z) in enumerate(contexts):
            for oi,(a,c) in enumerate(product(range(2),range(2))):
                uv=slack[ci][oi]; r1=np.zeros(total); r2=np.zeros(total)
                for b in range(3):
                    for al in range(4):
                        if fA(al,x)!=a: continue
                        for ga in range(4):
                            if fC(ga,z)!=c: continue
                            r1[idx[(y1,b,al,ga)]]+=1; r1[idx[(y2,b,al,ga)]]-=1
                            r2[idx[(y1,b,al,ga)]]-=1; r2[idx[(y2,b,al,ga)]]+=1
                r1[uv]-=1; r2[uv]-=1; Aub+=[r1,r2]
            r=np.zeros(total)
            for uv in slack[ci]: r[uv]=1
            r[delta]=-2; Aub.append(r)
        self.Aub=np.array(Aub); self.bub=np.zeros(len(Aub))
        A=np.zeros((len(eq),total))
        for i,r in enumerate(eq): A[i,:self.N]=r
        self.Aeq=A; self.total=total
        self.c=np.zeros(total); self.c[delta]=1
    def solve(self, PA, PB, PC, psi=PSI):
        PQ={}
        for x,y,z in product(range(2),range(3),range(2)):
            for a,b,c in product(range(2),range(3),range(2)):
                M=np.kron(np.kron(PA[x][a],PB[y][b]),PC[z][c])
                PQ[(x,y,z,a,b,c)]=np.real(psi.conj()@M@psi)
        beq=[]
        for spec in self.specs:
            if spec[0]=='norm': beq.append(1.0)
            elif spec[0]=='AB':
                _,x,y,a,b=spec
                beq.append(sum(PQ[(x,y,0,a,b,c)] for c in range(2)))
            else:
                _,y,z,b,c=spec
                beq.append(sum(PQ[(0,y,z,a,b,c)] for a in range(2)))
        res=solve_lp(self.c,A_ub=self.Aub,b_ub=self.bub,A_eq=self.Aeq,b_eq=np.array(beq),
                     bounds=[(0,None)]*self.total,context="barnea tripartite")
        return res.fun

def bell_measurement_coarse(u_local=None):
    """3-outcome coarse-grained Bell measurement on b1,b2: {Phi+}, {Phi-}, {Psi+ or Psi-}.
       Optionally pre-rotate b1 by 2x2 unitary u_local."""
    ph_p=np.array([1,0,0,1],dtype=complex)/np.sqrt(2)
    ph_m=np.array([1,0,0,-1],dtype=complex)/np.sqrt(2)
    ps_p=np.array([0,1,1,0],dtype=complex)/np.sqrt(2)
    ps_m=np.array([0,1,-1,0],dtype=complex)/np.sqrt(2)
    P=[np.outer(ph_p,ph_p.conj()),np.outer(ph_m,ph_m.conj()),
       np.outer(ps_p,ps_p.conj())+np.outer(ps_m,ps_m.conj())]
    if u_local is not None:
        U=np.kron(u_local,I2)
        P=[U@p@U.conj().T for p in P]
    return P

def rot(theta):
    return np.array([[np.cos(theta/2),-np.sin(theta/2)],[np.sin(theta/2),np.cos(theta/2)]],dtype=complex)
