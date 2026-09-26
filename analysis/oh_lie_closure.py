#!/usr/bin/env python3
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
sys.path.insert(0,str(ROOT/"src"))
from independent_radical_benchmark import generate
from spin_fci import read_fcidump,pure_spin_hamiltonian
from physical_pool import build_pool

OUT=ROOT/"analysis"/"oh_lie_closure"
OUT.mkdir(parents=True,exist_ok=True)

def vec_skew(A):
    A=np.asarray(np.real_if_close(A),float)
    iu=np.triu_indices(A.shape[0],1)
    return A[iu]

def add_orth(A,Q,tol=1e-10):
    v=vec_skew(A).copy()
    # reorthogonalized MGS
    for _ in range(2):
        for q in Q:
            v-=np.dot(q,v)*q
    n=np.linalg.norm(v)
    if n<=tol:
        return False
    Q.append(v/n)
    return True

def main():
    fcidump,cfg,meta=generate("OH",OUT)
    d=read_fcidump(str(fcidump));na=(d.nelec+d.ms2)//2;nb=d.nelec-na
    basis,H,U,Hd,Sp=pure_spin_hamiltonian(d,na,nb,target_s=.5)
    pool=build_pool(basis,U,Sp,d.norb,na,nb,max_rank=2)
    mats=[]
    for e in pool:
        X=U.conj().T@(e.matrix@U)
        X=.5*(X-X.conj().T)
        mats.append(np.asarray(np.real_if_close(X),float))

    Q=[]
    basis_mats=[]
    for A in mats:
        if add_orth(A,Q):
            basis_mats.append(A/np.linalg.norm(vec_skew(A)))
    initial=len(Q)
    target=U.shape[1]*(U.shape[1]-1)//2
    rounds=[{"round":0,"dimension":initial}]
    frontier=list(basis_mats)

    # Repeated commutators with the original generator set span the generated Lie algebra.
    for rnd in range(1,30):
        new=[]
        current_basis=list(basis_mats)
        for A in frontier:
            for G in mats:
                C=A@G-G@A
                if np.linalg.norm(C)<1e-12:
                    continue
                if add_orth(C,Q):
                    C=C/np.linalg.norm(vec_skew(C))
                    basis_mats.append(C);new.append(C)
                    if len(Q)>=target:
                        break
            if len(Q)>=target:
                break
        rounds.append({"round":rnd,"dimension":len(Q),"new":len(new)})
        print(f"round {rnd} dim {len(Q)} new {len(new)} target {target}",flush=True)
        if not new or len(Q)>=target:
            break
        frontier=new

    out={
        "system":"OH STO-3G CAS(7e,5o)",
        "doublet_dimension":int(U.shape[1]),
        "pool_size":len(pool),
        "initial_linear_span_dimension":initial,
        "max_real_skew_dimension":target,
        "lie_closure_dimension":len(Q),
        "full_so_dimension_reached":bool(len(Q)==target),
        "rounds":rounds,
        "interpretation_if_full":"Within real state space, the generated Lie algebra is so(d); the pool is therefore state-universal on the connected real unit sphere. A stalled first-gradient ADAPT search is then a search/optimization roadblock rather than a linear expressivity proof.",
    }
    (OUT/"SUMMARY.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__":main()
