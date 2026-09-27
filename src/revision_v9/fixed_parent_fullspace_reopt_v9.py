#!/usr/bin/env python3
"""Fixed-parent full-space reoptimization for the canonical projected Cu sequence.

Keeps the selected parent generators and their order fixed, frees all amplitudes,
and minimizes the full fixed-Ms energy subject to a hard target-doublet-weight
constraint p_D >= target using SLSQP with analytic gradients.

Default campaign:
  * p_D >= 0.999 from zero + three N(0,0.05) starts and the projected amplitudes;
  * p_D >= 0.9999 from zero and the projected amplitudes.

The determinant-excitation exponentials are evaluated exactly as disjoint
two-level rotations; this is algebraically equivalent to expm_multiply for
these parent generators.
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import expm_multiply

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/"src"))
from spin_fci import (read_fcidump,pure_spin_hamiltonian,hf_determinant,
                      excitation_specs_from_hf,excitation_generator,
                      spin_square_expectation)

def setup(fcidump, checkpoint):
    ck=json.load(open(checkpoint,"r",encoding="utf-8"))
    D=read_fcidump(str(fcidump))
    na=(D.nelec+D.ms2)//2; nb=D.nelec-na
    basis,H,U,Hd,Sp=pure_spin_hamiltonian(D,na,nb,.5)
    specs=excitation_specs_from_hf(D.norb,na,nb,max_rank=3)
    mats=[excitation_generator(basis,*specs[k]).astype(complex)
          for k in ck["selected_indices"]]
    rots=[]
    for A in mats:
        c=A.tocoo()
        q=[(r,col,v) for r,col,v in zip(c.row,c.col,c.data)
           if r<col and abs(v)>0]
        rots.append((np.array([x[0] for x in q],int),
                     np.array([x[1] for x in q],int),
                     np.array([float(np.real(x[2])) for x in q])))
    pos={int(d):i for i,d in enumerate(basis)}
    psi0=np.zeros(len(basis),complex)
    psi0[pos[hf_determinant(D.norb,na,nb)]]=1
    th0=np.asarray(ck["theta"],float)
    ed,vd=np.linalg.eigh(Hd); ED=float(ed[0])
    exactD=U@vd[:,0]; exactD/=np.linalg.norm(exactD)
    phi=U.conj().T@psi0
    for t,A in zip(th0,mats):
        Ad=U.conj().T@(A@U); Ad=.5*(Ad-Ad.conj().T)
        phi=expm_multiply(float(t)*Ad,phi)
    proj=U@phi; proj/=np.linalg.norm(proj)
    return basis,H,U,Sp,mats,rots,psi0,th0,ED,exactD,proj

def apply_rot(v,rot,t):
    i,j,s=rot; y=v.copy(); c=math.cos(t); z=math.sin(t)
    xi=v[i]; xj=v[j]
    y[i]=c*xi+s*z*xj; y[j]=-s*z*xi+c*xj
    return y

def apply_A(v,rot):
    i,j,s=rot; y=np.zeros_like(v)
    y[i]=s*v[j]; y[j]=-s*v[i]
    return y

def campaign(fcidump, checkpoint, seed=9272026, sigma=.05):
    basis,H,U,Sp,mats,rots,psi0,th0,ED,exactD,proj=setup(fcidump,checkpoint)
    def Q(v): return v-U@(U.conj().T@v)
    def opfg(th,which):
        states=[psi0]; v=psi0
        for t,r in zip(th,rots):
            v=apply_rot(v,r,float(t)); states.append(v)
        psi=states[-1]; w=H@psi if which=="E" else Q(psi)
        f=float(np.vdot(psi,w).real); g=np.empty(len(th))
        for k in range(len(th)-1,-1,-1):
            g[k]=2*np.real(np.vdot(w,apply_A(states[k+1],rots[k])))
            w=apply_rot(w,rots[k],-float(th[k]))
        return f,g
    def evaluate(th):
        p=psi0.copy()
        for t,r in zip(th,rots): p=apply_rot(p,r,float(t))
        p/=np.linalg.norm(p); E=float(np.vdot(p,H@p).real)
        q=float(np.linalg.norm(Q(p))**2); c=U.conj().T@p
        pd=U@(c/np.linalg.norm(c)); Epd=float(np.vdot(pd,H@pd).real)
        return dict(error_mEh=(E-ED)*1e3,doublet_weight=1-q,
                    S2=float(spin_square_expectation(p,Sp,ms=.5)),
                    fidelity_to_projected=float(abs(np.vdot(proj,p))**2),
                    postprojected_error_mEh=(Epd-ED)*1e3,
                    fidelity_projected_vs_postprojected=float(abs(np.vdot(proj,pd))**2),
                    fidelity_exact_doublet_vs_postprojected=float(abs(np.vdot(exactD,pd))**2))
    rng=np.random.default_rng(seed)
    starts={"zero":np.zeros_like(th0),
            "zero_perturbation_1":rng.normal(0,sigma,len(th0)),
            "zero_perturbation_2":rng.normal(0,sigma,len(th0)),
            "zero_perturbation_3":rng.normal(0,sigma,len(th0)),
            "projected_amplitudes":th0.copy()}
    results=[]
    for target,names in [(0.999,list(starts)),(0.9999,["zero","projected_amplitudes"])]:
        qmax=1-target
        for name in names:
            def cfun(x):
                q,g=opfg(x,"Q"); return qmax-q,-g
            res=minimize(lambda x:opfg(x,"E"),starts[name],jac=True,method="SLSQP",
                         constraints=[{"type":"ineq",
                                       "fun":lambda x:cfun(x)[0],
                                       "jac":lambda x:cfun(x)[1]}],
                         options={"maxiter":1000,"ftol":1e-12,"disp":False})
            results.append(dict(target_pD=target,initialization=name,
                                success=bool(res.success),status=int(res.status),
                                nit=int(res.nit),nfev=int(res.nfev),
                                message=str(res.message),**evaluate(np.asarray(res.x,float)),
                                theta=list(map(float,res.x))))
    return results

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fcidump",default=str(ROOT/"data/18q/active.FCIDUMP"))
    ap.add_argument("--checkpoint",default=str(ROOT/"results/revision_v6/cu/canonical_projected_121_checkpoint_compact.json"))
    ap.add_argument("--output",default=str(ROOT/"results/revision_v9/cu/fixed_parent_reoptimization_v9.json"))
    a=ap.parse_args()
    out=campaign(Path(a.fcidump),Path(a.checkpoint))
    p=Path(a.output); p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(out,indent=2),encoding="utf-8")
    summary=[{k:v for k,v in x.items() if k!="theta"} for x in out]
    print(json.dumps(summary,indent=2))

if __name__=="__main__":
    main()
