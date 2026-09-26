#!/usr/bin/env python3
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import expm_multiply

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
sys.path.insert(0,str(ROOT/"src"))
from independent_radical_benchmark import generate
from spin_fci import read_fcidump,pure_spin_hamiltonian,hf_determinant,spin_square_expectation
from physical_pool import build_pool

OUT=ROOT/"analysis"/"oh_escape_audit";OUT.mkdir(parents=True,exist_ok=True)

def main():
    fcidump,cfg,meta=generate("OH",OUT)
    d=read_fcidump(str(fcidump));na=(d.nelec+d.ms2)//2;nb=d.nelec-na;ms=.5*(na-nb)
    basis,H,U,Hd,Sp=pure_spin_hamiltonian(d,na,nb,target_s=.5)
    E0=float(np.linalg.eigvalsh(Hd)[0])
    pool=build_pool(basis,U,Sp,d.norb,na,nb,max_rank=2)
    pos={int(x):i for i,x in enumerate(basis)}[hf_determinant(d.norb,na,nb)]
    psi0=np.zeros(len(basis),complex);psi0[pos]=1

    selected=[];theta=np.zeros(0)
    def B(i):return pool[i].matrix
    def build(th,sel=None):
        if sel is None:sel=selected
        v=psi0.copy()
        for t,i in zip(th,sel):v=expm_multiply(float(t)*B(i),v)
        return v/np.linalg.norm(v)
    def objgrad(th,sel=None):
        if sel is None:sel=selected
        states=[psi0]
        for t,i in zip(th,sel):states.append(expm_multiply(float(t)*B(i),states[-1]))
        psi=states[-1];w=H@psi;E=float(np.vdot(psi,w).real);g=np.zeros(len(th))
        for k in range(len(th)-1,-1,-1):
            Bi=B(sel[k]);g[k]=2*np.real(np.vdot(w,Bi@states[k+1]));w=expm_multiply(-float(th[k])*Bi,w)
        return E,g

    # Reproduce the six meaningful first-order ADAPT selections before the gradient stall.
    trace=[]
    for it in range(6):
        psi=build(theta);w=H@psi
        grads=np.array([abs(2*np.real(np.vdot(w,e.matrix@psi))) for e in pool])
        j=int(np.argmax(grads));selected.append(j);theta=np.append(theta,0.)
        res=minimize(objgrad,theta,method="L-BFGS-B",jac=True,
                     options={"gtol":1e-10,"ftol":1e-14,"maxiter":1200,"maxls":80})
        theta=np.asarray(res.x,float);E=float(res.fun)
        trace.append({"iter":it+1,"selected":j,"max_gradient_before":float(grads[j]),
                      "E":E,"dev_mEh":(E-E0)*1e3,"optimizer_success":bool(res.success)})

    psi=build(theta);E=float(np.vdot(psi,H@psi).real);w=H@psi
    rawgr=np.array([2*np.real(np.vdot(w,e.matrix@psi)) for e in pool])
    gmax=float(np.max(np.abs(rawgr)))
    stalled={"selected":selected,"theta":theta.tolist(),"E":E,"dev_mEh":(E-E0)*1e3,
             "max_abs_gradient":gmax,"S2":spin_square_expectation(psi,Sp,ms=ms),
             "doublet_weight":float(np.linalg.norm(U.conj().T@psi)**2)}

    # Two-operator finite-step lookahead. This is a diagnostic, not a claimed algorithm.
    epsilons=[0.01,0.02,0.05,0.10]
    best=None
    for eps in epsilons:
        for i in range(len(pool)):
            for si in (-1.0,1.0):
                v1=expm_multiply(si*eps*B(i),psi)
                for j in range(len(pool)):
                    for sj in (-1.0,1.0):
                        v2=expm_multiply(sj*eps*B(j),v1);v2/=np.linalg.norm(v2)
                        e=float(np.vdot(v2,H@v2).real)
                        rec=(e,i,j,si*eps,sj*eps)
                        if best is None or e<best[0]:best=rec
        print("eps",eps,"best dev",(best[0]-E0)*1e3,"pair",best[1:3],flush=True)

    bestE,i,j,ti,tj=best
    sel2=selected+[int(i),int(j)]
    th2=np.concatenate([theta,[ti,tj]])
    res=minimize(lambda x:objgrad(x,sel2),th2,method="L-BFGS-B",jac=True,
                 options={"gtol":1e-11,"ftol":1e-14,"maxiter":3000,"maxls":100})
    th2=np.asarray(res.x,float)
    psi2=build(th2,sel2);E2=float(np.vdot(psi2,H@psi2).real)
    escape={"lookahead_best_preopt_E":float(bestE),
            "lookahead_best_preopt_dev_mEh":float((bestE-E0)*1e3),
            "pair":[int(i),int(j)],"initial_pair_angles":[float(ti),float(tj)],
            "post_global_opt_E":E2,"post_global_opt_dev_mEh":(E2-E0)*1e3,
            "optimizer_success":bool(res.success),"optimizer_message":str(res.message),
            "theta":th2.tolist(),"selected":sel2,
            "S2":spin_square_expectation(psi2,Sp,ms=ms),
            "doublet_weight":float(np.linalg.norm(U.conj().T@psi2)**2)}

    out={"meta":meta,"trace":trace,"stalled":stalled,"escape":escape}
    (OUT/"SUMMARY.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__":main()
