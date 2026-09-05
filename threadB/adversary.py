import sys
import numpy as np
from scipy.optimize import linprog
from itertools import product

# --------------------------------------------------------------- solver policy
# HiGHS defaults are 1e-7 primal/dual feasibility; the manuscript's numerical
# record is stated at ~1e-9, so the tolerances are set explicitly rather than
# inherited.  FEAS_TOL is the acceptance threshold applied to the residuals
# actually computed from the returned solution -- a success flag alone is not
# treated as evidence.
HIGHS_OPTIONS = {'primal_feasibility_tolerance': 1e-9,
                 'dual_feasibility_tolerance': 1e-9}
FEAS_TOL = 1e-8


class LPFailure(RuntimeError):
    """Raised when an LP solve fails or its solution violates the stated tolerance.

    Never caught-and-zeroed: a failed solve must not be reported as Sigma = 0,
    which would silently manufacture support for the ceiling conjecture."""


DIAG = {'solves': 0, 'max_eq_resid': 0.0, 'max_ub_viol': 0.0, 'max_bound_viol': 0.0}


def _mv(A, x):
    return np.asarray(A @ x).ravel()


def solve_lp(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None,
             context='', tol=FEAS_TOL):
    """linprog(method='highs') with explicit tolerances, hard failure on an
    unsuccessful or nonfinite solve, and residual diagnostics computed from the
    returned solution.  Returns the OptimizeResult."""
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=bounds, method='highs', options=HIGHS_OPTIONS)
    where = f" [{context}]" if context else ""
    if not res.success:
        raise LPFailure(f"LP solve failed{where}: status={res.status} "
                        f"({getattr(res, 'message', '')!r})")
    if res.fun is None or not np.isfinite(res.fun):
        raise LPFailure(f"LP returned a missing/nonfinite objective{where}: {res.fun!r}")
    if res.x is None:
        raise LPFailure(f"LP returned no solution vector{where}")
    x = np.asarray(res.x, dtype=float)
    if not np.all(np.isfinite(x)):
        raise LPFailure(f"LP returned a nonfinite solution vector{where}")

    eq_r = ub_v = bd_v = 0.0
    if A_eq is not None and b_eq is not None:
        eq_r = float(np.max(np.abs(_mv(A_eq, x) - np.asarray(b_eq, float)))) if x.size else 0.0
    if A_ub is not None and b_ub is not None:
        ub_v = float(np.max(np.maximum(_mv(A_ub, x) - np.asarray(b_ub, float), 0.0)))
    if bounds is not None:
        lo = np.array([-np.inf if b[0] is None else b[0] for b in bounds], dtype=float)
        hi = np.array([np.inf if b[1] is None else b[1] for b in bounds], dtype=float)
        bd_v = float(max(np.max(np.maximum(lo - x, 0.0)), np.max(np.maximum(x - hi, 0.0))))
    worst = max(eq_r, ub_v, bd_v)
    if worst > tol:
        raise LPFailure(f"LP solution violates tolerance {tol:g}{where}: "
                        f"eq residual {eq_r:.3e}, ub violation {ub_v:.3e}, "
                        f"bound violation {bd_v:.3e}")
    DIAG['solves'] += 1
    DIAG['max_eq_resid'] = max(DIAG['max_eq_resid'], eq_r)
    DIAG['max_ub_viol'] = max(DIAG['max_ub_viol'], ub_v)
    DIAG['max_bound_viol'] = max(DIAG['max_bound_viol'], bd_v)
    return res


def dependency_report():
    import scipy, platform
    try:
        import sympy; sv = sympy.__version__
    except ImportError:
        sv = "not installed"
    return (f"python {platform.python_version()}, numpy {np.__version__}, "
            f"scipy {scipy.__version__}, sympy {sv}; "
            f"HiGHS options {HIGHS_OPTIONS}")


def diagnostics_report():
    return (f"{DIAG['solves']} LP solves | max equality residual "
            f"{DIAG['max_eq_resid']:.3e} | max inequality violation "
            f"{DIAG['max_ub_viol']:.3e} | max bound violation "
            f"{DIAG['max_bound_viol']:.3e} | acceptance tolerance {FEAS_TOL:g}")

I2=np.eye(2); X=np.array([[0,1],[1,0]]); Z=np.array([[1,0],[0,-1]])
def kron(*ops):
    out=np.array([[1.0+0j]])
    for o in ops: out=np.kron(out,o)
    return out
def projs(O):
    vals,vecs=np.linalg.eigh(O); P={}
    for o in (0,1):
        t=1 if o==0 else -1
        P[o]=sum(np.outer(vecs[:,k],vecs[:,k].conj()) for k in range(2) if abs(vals[k]-t)<1e-9)
    return P

class SigmaLP:
    """Fixed scenario: A nA settings, B nB, C nC, D 2. Build structure once; per-state solve."""
    def __init__(self, nA, nB, nC, PA, PB, PC, PD):
        self.nA,self.nB,self.nC=nA,nB,nC
        self.PA,self.PB,self.PC,self.PD=PA,PB,PC,PD
        nFB=2**nB; nFC=2**nC
        idx={}; k=0
        for x,w,a,d,fb,fg in product(range(nA),range(2),range(2),range(2),range(nFB),range(nFC)):
            idx[(x,w,a,d,fb,fg)]=k; k+=1
        self.idx=idx; self.N=k
        fB=lambda fn,y:(fn>>y)&1; fC=lambda fn,z:(fn>>z)&1
        self.fB,self.fC=fB,fC
        eq_rows=[]
        self.marg_specs=[]
        for x,w in product(range(nA),range(2)):
            r=np.zeros(self.N)
            for a,d,fb,fg in product(range(2),range(2),range(nFB),range(nFC)):
                r[idx[(x,w,a,d,fb,fg)]]=1
            eq_rows.append(r); self.marg_specs.append(('norm',))
        for x,y,w,a,b,d in product(range(nA),range(nB),range(2),range(2),range(2),range(2)):
            r=np.zeros(self.N)
            for fb in range(nFB):
                if fB(fb,y)!=b: continue
                for fg in range(nFC): r[idx[(x,w,a,d,fb,fg)]]+=1
            eq_rows.append(r); self.marg_specs.append(('ABD',x,y,w,a,b,d))
        for x,z,w,a,c,d in product(range(nA),range(nC),range(2),range(2),range(2),range(2)):
            r=np.zeros(self.N)
            for fg in range(nFC):
                if fC(fg,z)!=c: continue
                for fb in range(nFB): r[idx[(x,w,a,d,fb,fg)]]+=1
            eq_rows.append(r); self.marg_specs.append(('ACD',x,z,w,a,c,d))
        contexts=[]
        for x1 in range(nA):
            for x2 in range(x1+1,nA):
                for w,y,z in product(range(2),range(nB),range(nC)):
                    contexts.append(('x',x1,x2,w,y,z))
        for x,y,z in product(range(nA),range(nB),range(nC)):
            contexts.append(('w',x,y,z))
        total=self.N; ctx_slack={}
        for ci,_ in enumerate(contexts):
            ctx_slack[ci]=list(range(total,total+8)); total+=8
        delta=total; total+=1
        def mcols(x,y,z,w,a,b,c,d):
            return [idx[(x,w,a,d,fb,fg)] for fb in range(nFB) if fB(fb,y)==b for fg in range(nFC) if fC(fg,z)==c]
        Aub=[]; 
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
                r1[uv]-=1; r2[uv]-=1; Aub+=[r1,r2]
            r=np.zeros(total)
            for uv in ctx_slack[ci]: r[uv]=1
            r[delta]=-2; Aub.append(r)
        self.Aub=np.array(Aub); self.bub=np.zeros(len(Aub))
        Aeq=np.zeros((len(eq_rows),total))
        for i,r in enumerate(eq_rows): Aeq[i,:self.N]=r
        self.Aeq=Aeq; self.total=total
        self.cobj=np.zeros(total); self.cobj[delta]=1
    def solve(self, psi, context=''):
        """Returns Sigma_HIC for the state psi.  Raises LPFailure if the solve
        does not succeed -- callers must not substitute 0 for a failure."""
        PQ={}
        for x in range(self.nA):
            for y in range(self.nB):
                for z in range(self.nC):
                    for w in range(2):
                        for a,b,c,d in product(range(2),repeat=4):
                            M=kron(self.PA[x][a],self.PB[y][b],self.PC[z][c],self.PD[w][d])
                            PQ[(x,y,z,w,a,b,c,d)]=np.real(psi.conj()@M@psi)
        beq=[]
        for spec in self.marg_specs:
            if spec[0]=='norm': beq.append(1.0)
            elif spec[0]=='ABD':
                _,x,y,w,a,b,d=spec
                beq.append(sum(PQ[(x,y,0,w,a,b,c,d)] for c in range(2)))
            else:
                _,x,z,w,a,c,d=spec
                beq.append(sum(PQ[(x,0,z,w,a,b,c,d)] for b in range(2)))
        res=solve_lp(self.cobj,A_ub=self.Aub,b_ub=self.bub,A_eq=self.Aeq,b_eq=np.array(beq),
                     bounds=[(0,None)]*self.total,
                     context=context or f"SigmaLP({self.nA},{self.nB},{self.nC})")
        return res.fun

def cluster4(phase=np.pi):
    plus=np.array([1,1])/np.sqrt(2); psi=plus
    for _ in range(3): psi=np.kron(psi,plus)
    psi=psi.astype(complex)
    def CP(n,i,j,th):
        M=np.eye(2**n,dtype=complex)
        for s in range(2**n):
            if (s>>(n-1-i))&1 and (s>>(n-1-j))&1: M[s,s]=np.exp(1j*th)
        return M
    return CP(4,0,1,phase)@CP(4,1,2,phase)@CP(4,2,3,phase)@psi

PA=[projs(X),projs(Z)]; PB=[projs((Z+X)/np.sqrt(2)),projs((Z-X)/np.sqrt(2))]
PC=[projs(Z),projs(X)]; PD=[projs(X),projs(Z)]
