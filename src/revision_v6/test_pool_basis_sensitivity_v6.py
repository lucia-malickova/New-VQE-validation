#!/usr/bin/env python3
"""Controlled null-space-basis sensitivity test for generalized physical ADAPT."""
from __future__ import annotations
import argparse, json
from collections import defaultdict
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize
from scipy.sparse.linalg import expm_multiply

from spin_fci import read_fcidump, pure_spin_hamiltonian, hf_determinant
from physical_pool import build_pool, Entry

def rotated_pool(pool, seed):
    if seed is None:
        return list(pool), 0
    rng=np.random.default_rng(seed); groups=defaultdict(list)
    for i,e in enumerate(pool): groups[e.group_key].append(i)
    out=list(pool); nrot=0
    for _,idx in groups.items():
        if len(idx)!=2: continue
        i,j=idx; e1,e2=pool[i],pool[j]
        a=float(rng.uniform(0,2*np.pi)); c,s=np.cos(a),np.sin(a)
        for target,x,y in ((i,c,s),(j,-s,c)):
            B=(x*e1.matrix+y*e2.matrix).tocsr(); B.eliminate_zeros()
            out[target]=Entry(target,e1.rank,e1.group_key,(),np.array([],float),(),B,
                              max(e1.leakage_ratio,e2.leakage_ratio),
                              max(e1.commutator_ratio,e2.commutator_ratio),
                              float(sp.linalg.norm(B)))
        nrot+=1
    return out,nrot

def run_one(fcidump,seed,max_ops=50):
    D=read_fcidump(str(fcidump)); na=(D.nelec+D.ms2)//2; nb=D.nelec-na
    basis,H,U,Hd,Sp=pure_spin_hamiltonian(D,na,nb,0.5)
    evals,evecs=np.linalg.eigh(Hd); E0=float(evals[0]); exact=U@evecs[:,0]; exact/=np.linalg.norm(exact)
    pos={int(x):i for i,x in enumerate(basis)}
    psi0=np.zeros(len(basis),complex); psi0[pos[hf_determinant(D.norb,na,nb)]]=1.
    pool0=build_pool(basis,U,Sp,D.norb,na,nb,max_rank=2); pool,nrot=rotated_pool(pool0,seed)
    selected=[]; theta=np.zeros(0); history=[]
    def build(th=None):
        th=theta if th is None else th; v=psi0.copy()
        for t,k in zip(th,selected): v=expm_multiply(float(t)*pool[k].matrix,v)
        return v
    def full_fg(th):
        states=[psi0]
        for t,k in zip(th,selected): states.append(expm_multiply(float(t)*pool[k].matrix,states[-1]))
        psi=states[-1]; w=H@psi; E=float(np.vdot(psi,w).real); g=np.zeros(len(th))
        for i in range(len(th)-1,-1,-1):
            B=pool[selected[i]].matrix; g[i]=2*np.real(np.vdot(w,B@states[i+1]))
            w=expm_multiply(-float(th[i])*B,w)
        return E,g
    def optimize_window(window=12,maxiter=80):
        nonlocal theta
        start=max(0,len(theta)-window); fixed=theta[:start].copy(); p=psi0.copy()
        for t,k in zip(fixed,selected[:start]): p=expm_multiply(float(t)*pool[k].matrix,p)
        sel=selected[start:]
        def fg(x):
            states=[p]
            for t,k in zip(x,sel): states.append(expm_multiply(float(t)*pool[k].matrix,states[-1]))
            psi=states[-1]; w=H@psi; E=float(np.vdot(psi,w).real); g=np.zeros(len(x))
            for i in range(len(x)-1,-1,-1):
                B=pool[sel[i]].matrix; g[i]=2*np.real(np.vdot(w,B@states[i+1]))
                w=expm_multiply(-float(x[i])*B,w)
            return E,g
        r=minimize(fg,theta[start:].copy(),jac=True,method="L-BFGS-B",
                   options={"maxiter":maxiter,"ftol":1e-12,"gtol":1e-8,"maxls":30})
        theta=np.concatenate([fixed,np.asarray(r.x,float)])
    def optimize_global(maxiter):
        nonlocal theta
        r=minimize(full_fg,theta,jac=True,method="L-BFGS-B",
                   options={"maxiter":maxiter,"ftol":1e-13,"gtol":1e-8,"maxls":40})
        theta=np.asarray(r.x,float)
    for step in range(1,max_ops+1):
        psi=build(); w=H@psi
        grads=np.array([abs(2*np.real(np.vdot(w,e.matrix@psi))) for e in pool])
        k=int(np.argmax(grads)); selected.append(k); theta=np.append(theta,0.0)
        optimize_window()
        if step%20==0: optimize_global(160)
        psi=build(); psi/=np.linalg.norm(psi); E=float(np.vdot(psi,H@psi).real)
        history.append({"step":step,"error_mEh":(E-E0)*1e3,"fidelity":float(abs(np.vdot(exact,psi))**2)})
    optimize_global(300); psi=build(); psi/=np.linalg.norm(psi); E=float(np.vdot(psi,H@psi).real)
    return {"rotation_seed":seed,"two_dimensional_groups_rotated":nrot,"n_ops":len(selected),
            "step30_error_mEh":next(x["error_mEh"] for x in history if x["step"]==30),
            "step50_error_before_final_global_mEh":history[-1]["error_mEh"],
            "step50_error_after_final_global_mEh":(E-E0)*1e3,
            "step50_fidelity_exact_doublet":float(abs(np.vdot(exact,psi))**2),
            "step50_residual_Eh":float(np.linalg.norm(H@psi-E*psi)),
            "step50_unique_generators":len(set(selected))}
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--fcidump",required=True)
    ap.add_argument("--seeds",nargs="*",type=int,default=[1,2,3]); ap.add_argument("--output",default="pool_basis_sensitivity_v6.json")
    a=ap.parse_args(); rows=[run_one(a.fcidump,None)]+[run_one(a.fcidump,s) for s in a.seeds]
    Path(a.output).write_text(json.dumps(rows,indent=2)); print(json.dumps(rows,indent=2))
if __name__=="__main__": main()
