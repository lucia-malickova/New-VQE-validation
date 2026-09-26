#!/usr/bin/env python3
from __future__ import annotations
import json, math, subprocess, sys
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
    evals,evecs=np.linalg.eigh(Hd);E0=float(evals[0])
    ng=int(np.sum(evals-evals[0] <= 1e-8))
    ground_sub=U@evecs[:,:ng]
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
    def diagnostics(psi):
        E=float(np.vdot(psi,H@psi).real)
        return {"E":E,"dev_mEh":(E-E0)*1e3,"S2":spin_square_expectation(psi,Sp,ms=ms),
                "doublet_weight":float(np.linalg.norm(U.conj().T@psi)**2),
                "weight_lowest_doublet_eigenspace":float(np.linalg.norm(ground_sub.conj().T@psi)**2)}

    trace=[]
    for it in range(6):
        psi=build(theta);w=H@psi
        grads=np.array([abs(2*np.real(np.vdot(w,e.matrix@psi))) for e in pool])
        j=int(np.argmax(grads));selected.append(j);theta=np.append(theta,0.)
        res=minimize(objgrad,theta,method="L-BFGS-B",jac=True,
                     options={"gtol":1e-10,"ftol":1e-14,"maxiter":1200,"maxls":80})
        theta=np.asarray(res.x,float);psi=build(theta)
        dg=diagnostics(psi)
        trace.append({"iter":it+1,"selected":j,"max_gradient_before":float(grads[j]),
                      **dg,"optimizer_success":bool(res.success)})

    psi=build(theta);w=H@psi
    rawgr=np.array([2*np.real(np.vdot(w,e.matrix@psi)) for e in pool])
    stalled={"selected":selected,"theta":theta.tolist(),"max_abs_gradient":float(np.max(np.abs(rawgr))),
             **diagnostics(psi)}

    # Exhaustive two-operator finite-step lookahead. Diagnostic only.
    epsilons=[0.01,0.02,0.05,0.10]
    per_eps=[]
    global_best=None
    for eps in epsilons:
        best=None
        for i in range(len(pool)):
            for si in (-1.0,1.0):
                v1=expm_multiply(si*eps*B(i),psi)
                for j in range(len(pool)):
                    for sj in (-1.0,1.0):
                        v2=expm_multiply(sj*eps*B(j),v1);v2/=np.linalg.norm(v2)
                        e=float(np.vdot(v2,H@v2).real)
                        rec=(e,i,j,si*eps,sj*eps)
                        if best is None or e<best[0]:best=rec
        e,i,j,ti,tj=best
        sel2=selected+[int(i),int(j)]; th2=np.concatenate([theta,[ti,tj]])
        res=minimize(lambda x:objgrad(x,sel2),th2,method="L-BFGS-B",jac=True,
                     options={"gtol":1e-11,"ftol":1e-14,"maxiter":3000,"maxls":100})
        x=np.asarray(res.x,float);v=build(x,sel2);dg=diagnostics(v)
        row={"epsilon":eps,"preopt_E":float(e),"preopt_dev_mEh":float((e-E0)*1e3),
             "pair":[int(i),int(j)],"initial_pair_angles":[float(ti),float(tj)],
             "postopt_theta":x.tolist(),"postopt_selected":sel2,
             "optimizer_success":bool(res.success),"optimizer_message":str(res.message),**dg}
        per_eps.append(row)
        if global_best is None or dg["E"]<global_best["E"]:global_best=row
        print("eps",eps,"pair",(i,j),"predev",(e-E0)*1e3,"postdev",dg["dev_mEh"],flush=True)

    # Two-dimensional Hessian at the stalled state for the best escape pair,
    # with all prior amplitudes fixed. This tests whether a first-gradient stall
    # still has a second-order descent direction.
    bi,bj=global_best["pair"]
    def pairE(x,y):
        v=expm_multiply(float(x)*B(bi),psi)
        v=expm_multiply(float(y)*B(bj),v);v/=np.linalg.norm(v)
        return float(np.vdot(v,H@v).real)
    h=1e-3;e00=pairE(0,0)
    Hii=(pairE(h,0)-2*e00+pairE(-h,0))/h**2
    Hjj=(pairE(0,h)-2*e00+pairE(0,-h))/h**2
    Hij=(pairE(h,h)-pairE(h,-h)-pairE(-h,h)+pairE(-h,-h))/(4*h**2)
    Hess=np.array([[Hii,Hij],[Hij,Hjj]])
    hess_eigs=np.linalg.eigvalsh(Hess)
    hessian={"pair":[int(bi),int(bj)],"step":h,"matrix_Eh_per_rad2":Hess.tolist(),
             "eigenvalues_Eh_per_rad2":hess_eigs.tolist(),
             "negative_direction_exists":bool(np.min(hess_eigs)<0)}

    # Serialize the best eight-operator corrected physical ansatz and validate circuits.
    best=global_best
    ck={"algorithm":"OH two-operator lookahead diagnostic after first-gradient ADAPT stall",
        "selected_generators":[pool[i].as_dict() for i in best["postopt_selected"]],
        "theta":best["postopt_theta"]}
    ckpath=OUT/"oh_escape_physical_ansatz.json";ckpath.write_text(json.dumps(ck,indent=2),encoding="utf-8")
    circuits=[]
    for reps in (1,2,4):
        op=OUT/f"oh_escape_circuit_suzuki2_r{reps}.json"
        cmd=[sys.executable,str(ROOT/"src"/"circuit_validate_physical.py"),str(ckpath),
             "--fcidump",str(fcidump),"--method","suzuki2","--reps",str(reps),
             "--optimization-level","1","--output",str(op)]
        subprocess.run(cmd,check=True,cwd=str(ROOT/"src"))
        x=json.loads(op.read_text())
        circuits.append({"reps":reps,"CX":x["transpiled_gate_counts"].get("cx",0),
                         "depth":x["transpiled_depth"],
                         "synthesis_error_mEh":x["actual_transpiled_circuit"]["energy_error_vs_exact_prefix_mEh"],
                         "fidelity_full":x["actual_transpiled_circuit"]["fidelity_full_with_exact_physical_prefix"],
                         "S2":x["actual_transpiled_circuit"]["S2_conditioned_on_target_Ms"],
                         "PASS_strict":x["PASS"]})

    out={"meta":meta,"lowest_doublet_degeneracy":ng,"trace":trace,"stalled":stalled,
         "lookahead_by_epsilon":per_eps,"best_escape":global_best,"hessian":hessian,"circuits":circuits}
    (OUT/"SUMMARY.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__":main()
