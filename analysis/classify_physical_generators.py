#!/usr/bin/env python3
from __future__ import annotations
import json, sys
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from spin_fci import read_fcidump,pure_spin_hamiltonian
from physical_pool import build_pool

OUT=ROOT/"analysis"/"generator_classification"
OUT.mkdir(parents=True,exist_ok=True)

def canon_pattern(coeffs,tol=1e-10):
    a=np.asarray(coeffs,float)
    nz=a[np.abs(a)>tol]
    return tuple(round(float(x),10) for x in nz)

def one(label,fcidump,checkpoint):
    d=read_fcidump(str(fcidump));na=(d.nelec+d.ms2)//2;nb=d.nelec-na
    basis,H,U,Hd,Sp=pure_spin_hamiltonian(d,na,nb,target_s=.5)
    pool=build_pool(basis,U,Sp,d.norb,na,nb,max_rank=2)
    ck=json.loads(Path(checkpoint).read_text())
    sel=list(map(int,ck["selected_pool_indices"]))
    entries=[pool[i] for i in sel]
    unique={i:pool[i] for i in sorted(set(sel))}
    rows=[]
    for pos,i in enumerate(sel,1):
        e=pool[i]
        rows.append({
            "position":pos,"pool_index":i,"rank":e.rank,
            "support_size":int(np.sum(np.abs(e.coeffs)>1e-12)),
            "coeffs":[float(x) for x in e.coeffs],
            "specs":[{"remove":list(r),"add":list(a)} for r,a in e.specs],
            "group_key":repr(e.group_key),
            "leakage_ratio":e.leakage_ratio,
            "commutator_ratio":e.commutator_ratio,
        })
    def dist(vals):return {str(k):int(v) for k,v in sorted(Counter(vals).items())}
    summary={
        "label":label,"pool_size":len(pool),"applications":len(sel),"unique":len(unique),
        "applications_rank":dist(e.rank for e in entries),
        "applications_support_size":dist(int(np.sum(np.abs(e.coeffs)>1e-12)) for e in entries),
        "unique_rank":dist(e.rank for e in unique.values()),
        "unique_support_size":dist(int(np.sum(np.abs(e.coeffs)>1e-12)) for e in unique.values()),
        "unique_coefficient_patterns":dist(canon_pattern(e.coeffs) for e in unique.values()),
        "max_selected_leakage":max(e.leakage_ratio for e in entries),
        "max_selected_commutator":max(e.commutator_ratio for e in entries),
        "selected":rows,
    }
    (OUT/f"{label}.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return summary

def main():
    res={}
    res["18q"]=one("18q",ROOT/"data/18q/active.FCIDUMP",ROOT/"data/18q/best50_compact.json")
    res["24q"]=one("24q",ROOT/"data/24q/active.FCIDUMP",ROOT/"data/24q/best61_compact.json")
    compact={k:{kk:vv for kk,vv in v.items() if kk!="selected"} for k,v in res.items()}
    (OUT/"SUMMARY.json").write_text(json.dumps(compact,indent=2),encoding="utf-8")
    print(json.dumps(compact,indent=2))

if __name__=="__main__":main()
