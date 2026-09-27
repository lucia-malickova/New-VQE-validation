#!/usr/bin/env python3
"""Direct representability test for the canonical 121-generator parent sequence.

Keeps the selected full-space parent generators and their order fixed, then
maximizes overlap with the lifted projected state. Reports unconstrained
fidelity optimization and a hard target-doublet-purity version p_D >= 0.999.
This is a local optimization diagnostic, not a proof of the global maximum.
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

def setup(fcidump,checkpoint):
    ck=json.load(open(checkpoint,"r",encoding="utf-8"))
    D=read_fcidump(str(fcidump)); na=(D.nelec+D.ms2)//2; nb=D.nelec-na
    basis,H,U,Hd,Sp=pure_spin_hamiltonian(D,na,nb,.5)
    specs=excitation_specs_from_hf(D.norb,na,nb,max_rank=3)
    mats=[excitation_generator(basis,*specs[k]).astype(complex)
          for k in ck["selected_indices"]]
    rots=[]
    for A in mats:
        c=A.tocoo(); q=[(r,col,v) for r,col,v in zip(c.row,c.col,c.data)
                         if r<col and abs(v)>0]
        rots.append((np.array([x[0] for x in q],int),
                     np.array([x[1] for x in q],int),
                     np.array([float(np.real(x[2])) for x in q])))
    pos={int(d):i for i,d in enumerate(basis)}
    psi0=np.zeros(len(basis),complex)
    psi0[pos[hf_determinant(D.norb,na,nb)]]=1
    th0=np.asarray(ck["theta"],float)
    ed,vd=np.linalg.eigh(Hd); ED=float(ed[0]); exactD=U@vd[:,0]; exactD/=np.linalg.norm(exactD)
    phi=U.conj().T@psi0
    for t,A in zip(th0,mats):
        Ad=U.conj().T@(A@U); Ad=.5*(Ad-Ad.conj().T)
        phi=expm_multiply(float(t)*Ad,phi)
    proj=U@phi; proj/=np.linalg.norm(proj)
    return H,U,Sp,rots,psi0,th0,ED,exactD,proj

def rotate(v,rot,t):
    i,j,s=rot; y=v.copy(); c=math.cos(t); z=math.sin(t)
    xi=v[i]; xj=v[j]; y[i]=c*xi+s*z*xj; y[j]=-s*z*xi+c*xj
    return y

def apply_A(v,rot):
    i,j,s=rot; y=np.zeros_like(v); y[i]=s*v[j]; y[j]=-s*v[i]; return y

def campaign(fcidump,checkpoint,seed=1977):
    H,U,Sp,rots,psi0,th0,ED,exactD,proj=setup(fcidump,checkpoint)
    def Q(v): return v-U@(U.conj().T@v)
    def expfg(th,op):
        states=[psi0]; v=psi0
        for t,r in zip(th,rots):
            v=rotate(v,r,float(t)); states.append(v)
        psi=states[-1]; w=op(psi); val=float(np.vdot(psi,w).real); grad=np.empty(len(th))
        for k in range(len(th)-1,-1,-1):
            grad[k]=2*np.real(np.vdot(w,apply_A(states[k+1],rots[k])))
            w=rotate(w,rots[k],-float(th[k]))
        return val,grad
    def fidelity_fg(th):
        val,grad=expfg(th,lambda v:proj*np.vdot(proj,v)); return -val,-grad
    def energy_fg(th): return expfg(th,lambda v:H@v)
    def q_fg(th): return expfg(th,Q)
    def evaluate(th):
        p=psi0.copy()
        for t,r in zip(th,rots): p=rotate(p,r,float(t))
        p/=np.linalg.norm(p); E=float(np.vdot(p,H@p).real); q=float(np.linalg.norm(Q(p))**2)
        c=U.conj().T@p; pd=U@(c/np.linalg.norm(c)); Epd=float(np.vdot(pd,H@pd).real)
        return dict(fidelity_to_projected=float(abs(np.vdot(proj,p))**2),
                    energy_error_mEh=(E-ED)*1e3,doublet_weight=1-q,
                    S2=float(spin_square_expectation(p,Sp,ms=.5)),
                    postprojected_error_mEh=(Epd-ED)*1e3,
                    fidelity_exact_doublet_vs_postprojected=float(abs(np.vdot(exactD,pd))**2))
    qmax=.001
    def constraint_fun(th):
        q,g=q_fg(th); return qmax-q,-g
    er=minimize(energy_fg,np.zeros_like(th0),jac=True,method="SLSQP",
                constraints=[{"type":"ineq","fun":lambda x:constraint_fun(x)[0],
                              "jac":lambda x:constraint_fun(x)[1]}],
                options={"maxiter":1000,"ftol":1e-12,"disp":False})
    if not er.success: raise RuntimeError("energy-constrained start failed: "+str(er.message))
    base=np.asarray(er.x,float); rng=np.random.default_rng(seed)
    starts={"energy_constrained_best":base,
            "energy_constrained_perturbation_1":base+rng.normal(0,.02,len(base)),
            "energy_constrained_perturbation_2":base+rng.normal(0,.02,len(base)),
            "energy_constrained_perturbation_3":base+rng.normal(0,.02,len(base))}
    results=[]; uncon=[]
    for name,x0 in starts.items():
        r=minimize(fidelity_fg,x0,jac=True,method="L-BFGS-B",
                   options={"maxiter":2000,"ftol":1e-15,"gtol":1e-10,"maxls":60,"maxcor":40})
        rec=dict(mode="unconstrained_fidelity",initialization=name,success=bool(r.success),
                 status=int(r.status),nit=int(r.nit),nfev=int(r.nfev),message=str(r.message),
                 **evaluate(np.asarray(r.x,float)),theta=list(map(float,r.x)))
        results.append(rec); uncon.append(rec)
    best=max(uncon,key=lambda x:x["fidelity_to_projected"])
    cstarts={"energy_constrained_best":base,
             "best_unconstrained_fidelity":np.asarray(best["theta"],float),
             "energy_constrained_perturbation":starts["energy_constrained_perturbation_1"]}
    for name,x0 in cstarts.items():
        def cfun(x):
            q,g=q_fg(x); return qmax-q,-g
        r=minimize(fidelity_fg,x0,jac=True,method="SLSQP",
                   constraints=[{"type":"ineq","fun":lambda x:cfun(x)[0],
                                 "jac":lambda x:cfun(x)[1]}],
                   options={"maxiter":1500,"ftol":1e-13,"disp":False})
        results.append(dict(mode="fidelity_with_pD_ge_0.999",initialization=name,
                            success=bool(r.success),status=int(r.status),nit=int(r.nit),
                            nfev=int(r.nfev),message=str(r.message),
                            **evaluate(np.asarray(r.x,float)),theta=list(map(float,r.x))))
    return results

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fcidump",default=str(ROOT/"data/18q/active.FCIDUMP"))
    ap.add_argument("--checkpoint",default=str(ROOT/"results/revision_v6/cu/canonical_projected_121_checkpoint_compact.json"))
    ap.add_argument("--output",default=str(ROOT/"results/revision_v10/cu/fixed_parent_fidelity_max_v10.json"))
    a=ap.parse_args(); out=campaign(Path(a.fcidump),Path(a.checkpoint))
    p=Path(a.output); p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps([{k:v for k,v in x.items() if k!="theta"} for x in out],indent=2))
if __name__=="__main__": main()
