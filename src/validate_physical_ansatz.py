#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,hashlib,numpy as np
from scipy.sparse.linalg import expm_multiply
from spin_fci import read_fcidump,pure_spin_hamiltonian,hf_determinant,spin_square_expectation
from physical_pool import entry_from_dict

def sha256(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('result_json');ap.add_argument('--fcidump',required=True);ap.add_argument('--output',default='physical_ansatz_validation.json');args=ap.parse_args()
 r=json.load(open(args.result_json,encoding='utf-8'));d=read_fcidump(args.fcidump);na=(d.nelec+d.ms2)//2;nb=d.nelec-na;ms=.5*(na-nb)
 basis,H,U,Hd,Sp=pure_spin_hamiltonian(d,na,nb,target_s=.5); ew,ev=np.linalg.eigh(Hd);E0=float(ew[0]);psi0=U@ev[:,0]
 pos={int(x):i for i,x in enumerate(basis)}[hf_determinant(d.norb,na,nb)];psi=np.zeros(len(basis),complex);psi[pos]=1
 entries=[entry_from_dict(x,basis) for x in r['selected_generators']];theta=np.asarray(r['theta'],float)
 for t,e in zip(theta,entries):psi=expm_multiply(float(t)*e.matrix,psi)
 psi/=np.linalg.norm(psi);E=float(np.vdot(psi,H@psi).real);s2=spin_square_expectation(psi,Sp,ms=ms);dc=U.conj().T@psi;pD=float(np.vdot(dc,dc).real);fid=float(abs(np.vdot(psi0,psi))**2);res=float(np.linalg.norm(H@psi-E*psi));dev=(E-E0)*1e3
 saved=(r.get('summary') or {}).get('energy_Eh'); match=True if saved is None else abs(E-float(saved))<1e-9; passed=bool(match and abs(s2-.75)<2e-8 and abs(pD-1)<2e-8)
 out={'result_json':args.result_json,'result_json_sha256':sha256(args.result_json),'fcidump':args.fcidump,'fcidump_sha256':sha256(args.fcidump),'n_selected':len(entries),'reconstructed_energy_Eh':E,'exact_doublet_Eh':E0,'dev_mEh':dev,'S2':s2,'doublet_weight':pD,'fidelity_exact_doublet':fid,'residual_norm_Eh':res,'saved_final_energy_Eh':saved,'PASS':passed}
 json.dump(out,open(args.output,'w'),indent=2);print(json.dumps(out,indent=2))
 if not passed:raise SystemExit(2)
if __name__=='__main__':main()
