#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,hashlib
from spin_fci import read_fcidump,pure_spin_hamiltonian
from physical_pool import build_pool

def sha256(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('compact_json');ap.add_argument('--fcidump',required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
 r=json.load(open(args.compact_json,encoding='utf-8')); d=read_fcidump(args.fcidump); na=(d.nelec+d.ms2)//2; nb=d.nelec-na
 basis,H,U,Hd,Sp=pure_spin_hamiltonian(d,na,nb,target_s=.5); pool=build_pool(basis,U,Sp,d.norb,na,nb,max_rank=2)
 idx=[int(x) for x in r['selected_pool_indices']]
 if max(idx)>=len(pool):raise SystemExit('pool index out of range')
 out=dict(r); out['selected_generators']=[pool[i].as_dict() for i in idx]
 out['metadata']=dict(out.get('metadata',{}));out['metadata']['fcidump']=args.fcidump;out['metadata']['fcidump_sha256']=sha256(args.fcidump)
 json.dump(out,open(args.output,'w',encoding='utf-8'),indent=2)
 print(f'materialized {len(idx)} generators from pool size {len(pool)} -> {args.output}')
if __name__=='__main__':main()
