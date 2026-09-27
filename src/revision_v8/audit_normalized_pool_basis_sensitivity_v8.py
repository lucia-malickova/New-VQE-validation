#!/usr/bin/env python3
"""Normalization-controlled null-space-basis sensitivity audit."""
from __future__ import annotations
import argparse,csv,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import expm_multiply
from spin_fci import read_fcidump,pure_spin_hamiltonian,hf_determinant
from physical_pool import build_pool

def write_csv(path,rows):
    with open(path,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--fcidump",required=True);ap.add_argument("--checkpoint",required=True)
    ap.add_argument("--n-random",type=int,default=100);ap.add_argument("--seed-start",type=int,default=1)
    ap.add_argument("--output-prefix",default="pool_basis_normalized_v8");a=ap.parse_args()
    ck=json.load(open(a.checkpoint,encoding="utf-8"));D=read_fcidump(a.fcidump);na=(D.nelec+D.ms2)//2;nb=D.nelec-na
    basis,H,U,Hd,Sp=pure_spin_hamiltonian(D,na,nb,0.5);pool=build_pool(basis,U,Sp,D.norb,na,nb,max_rank=2)
    selected=list(map(int,ck["selected_pool_indices"]));theta=np.asarray(ck["theta"],float)
    groups=defaultdict(list)
    for i,e in enumerate(pool):groups[e.group_key].append(i)
    gl=list(groups.values());two=[x for x in gl if len(x)==2];geom={};ca=[];fa=[];fr=[];hr=[]
    pos={int(d):i for i,d in enumerate(basis)};hf=np.zeros(len(basis));hf[pos[hf_determinant(D.norb,na,nb)]]=1.
    for inds in two:
        i,j=inds;e1,e2=pool[i],pool[j];c1=np.asarray(e1.coeffs);c2=np.asarray(e2.coeffs)
        C=np.array([[c1@c1,c1@c2],[c2@c1,c2@c2]],float);B1,B2=e1.matrix,e2.matrix
        F=np.array([[float(B1.multiply(B1).sum()),float(B1.multiply(B2).sum())],[float(B2.multiply(B1).sum()),float(B2.multiply(B2).sum())]],float)
        geom[(i,j)]=(C,F);ca.append(abs(C[0,1])/math.sqrt(C[0,0]*C[1,1]));fa.append(abs(F[0,1])/math.sqrt(F[0,0]*F[1,1]))
        ev=np.linalg.eigvalsh(F);fr.append(math.sqrt(ev[-1]/ev[0]));v1=B1@hf;v2=B2@hf
        A=np.array([[np.vdot(v1,v1).real,np.vdot(v1,v2).real],[np.vdot(v2,v1).real,np.vdot(v2,v2).real]])
        ev=np.linalg.eigvalsh(A)
        if ev[0]>1e-30:hr.append(math.sqrt(ev[-1]/ev[0]))
    seeds=range(a.seed_start,a.seed_start+a.n_random);angles={}
    for seed in seeds:
        rng=np.random.default_rng(seed);angles[seed]={tuple(x):float(rng.uniform(0,2*np.pi)) for x in two}
    psi=hf.astype(complex);detail=[]
    for pref in range(len(selected)+1):
        if pref:psi=expm_multiply(float(theta[pref-1])*pool[selected[pref-1]].matrix,psi)
        p=psi/np.linalg.norm(psi);w=H@p;signed=np.array([2*np.real(np.vdot(w,e.matrix@p)) for e in pool]);ag=np.abs(signed)
        canon=int(np.argmax(ag));cg=str(pool[canon].group_key)
        for metric in ("coeff","fro"):
            inv=-1.;ig=None
            for inds in gl:
                if len(inds)==1:val=float(abs(signed[inds[0]]))
                else:
                    i,j=inds;M=geom[(i,j)][0 if metric=="coeff" else 1];g=np.array([signed[i],signed[j]])
                    val=float(np.sqrt(max(0.,g@np.linalg.solve(M,g))))
                    if metric=="fro":val*=math.sqrt(.5*(M[0,0]+M[1,1]))
                if val>inv:inv=val;ig=str(pool[inds[0]].group_key)
            changed=0;tops=set();vals=[]
            for seed in seeds:
                gn=ag.copy()
                for inds in two:
                    i,j=inds;M=geom[(i,j)][0 if metric=="coeff" else 1];rho=M[0,1]/math.sqrt(M[0,0]*M[1,1])
                    alpha=angles[seed][tuple(inds)];c=math.cos(alpha);s=math.sin(alpha)
                    n1=math.sqrt(max(1e-30,1+2*c*s*rho));n2=math.sqrt(max(1e-30,1-2*c*s*rho))
                    gn[i]=abs((c*signed[i]+s*signed[j])/n1);gn[j]=abs((-s*signed[i]+c*signed[j])/n2)
                top=int(np.argmax(gn));tg=str(pool[top].group_key);changed+=int(tg!=cg);tops.add(tg);vals.append(float(gn[top]))
            detail.append({"prefix":pref,"normalization":metric,"canonical_max_gradient":float(ag[canon]),"invariant_group_opt_gradient":inv,
                "canonical_top_is_invariant_opt":cg==ig,"different_top_group_fraction":changed/float(a.n_random),"distinct_top_groups":len(tops),
                "random_top_gradient_min":min(vals),"random_top_gradient_median":float(np.median(vals)),"random_top_gradient_max":max(vals)})
    summary=[]
    for metric in ("coeff","fro"):
        x=[r for r in detail if r["normalization"]==metric]
        summary.append({"normalization_metric":metric,"prefix_states_tested":len(x),"random_bases_per_prefix":a.n_random,
        "prefixes_with_any_top_group_change":sum(r["different_top_group_fraction"]>0 for r in x),
        "max_fraction_random_bases_changing_top_group":max(r["different_top_group_fraction"] for r in x),
        "prefixes_where_canonical_top_differs_from_invariant_group_optimum":sum(not r["canonical_top_is_invariant_opt"] for r in x)})
    norms=[{"n_two_dimensional_groups":len(two),"abs_coefficient_overlap_min":min(ca),"abs_coefficient_overlap_max":max(ca),
    "abs_frobenius_overlap_min":min(fa),"abs_frobenius_overlap_max":max(fa),"max_raw_O2_frobenius_norm_ratio":max(fr),
    "max_raw_O2_HF_action_norm_ratio_nonzero":max(hr)}]
    p=Path(a.output_prefix);write_csv(p.with_name(p.name+"_allprefix.csv"),detail);write_csv(p.with_name(p.name+"_summary.csv"),summary);write_csv(p.with_name(p.name+"_normalization.csv"),norms)
    print(json.dumps({"summary":summary,"normalization":norms[0]},indent=2))
if __name__=="__main__":main()
