#!/usr/bin/env python3
"""State-specific representation-equivalence certificate for projected ADAPT."""
import argparse,json,math,sys
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import expm_multiply
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("checkpoint_json"); ap.add_argument("fcidump")
    ap.add_argument("--src",default="src"); ap.add_argument("--output",default="state_specific_representation_certificate.json")
    a=ap.parse_args(); sys.path.insert(0,a.src)
    from spin_fci import read_fcidump,pure_spin_hamiltonian,hf_determinant,excitation_specs_from_hf,excitation_generator
    ck=json.load(open(a.checkpoint_json)); D=read_fcidump(a.fcidump); na=(D.nelec+D.ms2)//2; nb=D.nelec-na
    basis,H,U,Hd,Sp=pure_spin_hamiltonian(D,na,nb,.5); pos={int(d):i for i,d in enumerate(basis)}
    psi0=np.zeros(len(basis),complex); psi0[pos[hf_determinant(D.norb,na,nb)]]=1.
    phi=U.conj().T@psi0; phi/=np.linalg.norm(phi); psi_proj=U@phi; psi_phys=psi0.copy()
    specs=excitation_specs_from_hf(D.norb,na,nb,max_rank=3); cache={}
    def A(k):
        if k not in cache: cache[k]=excitation_generator(basis,*specs[k]).astype(complex)
        return cache[k]
    rows=[]; Dseq=0.
    for i,(k,t) in enumerate(zip(ck["selected_indices"],ck["theta"]),1):
        Ak=A(int(k)); Ad=U.conj().T@(Ak@U); Ad=.5*(Ad-Ad.conj().T)
        local_full=expm_multiply(float(t)*Ak,psi_proj); phi_next=expm_multiply(float(t)*Ad,phi); proj_next=U@phi_next
        delta=float(np.linalg.norm(local_full-proj_next)); Dseq+=delta
        psi_phys=expm_multiply(float(t)*Ak,psi_phys); psi_proj=proj_next; phi=phi_next
        p=psi_phys/np.linalg.norm(psi_phys); q=psi_proj/np.linalg.norm(psi_proj)
        F=float(abs(np.vdot(q,p))**2); dproj=math.sqrt(max(0.,2.-2.*abs(np.vdot(q,p))))
        rows.append({"k":i,"pool_index":int(k),"theta":float(t),"local_delta":delta,"cumulative_Dseq":Dseq,
                     "prefix_fidelity":F,"projective_state_distance":dproj})
    Ffinal=rows[-1]["prefix_fidelity"]; dactual=rows[-1]["projective_state_distance"]
    Flow=max(0.,1.-Dseq*Dseq/2.)**2 if Dseq<=math.sqrt(2.) else 0.; Ptarget=max(0.,1.-Dseq*Dseq)
    out={"n_ops":len(rows),"D_seq":Dseq,"actual_projective_state_distance":dactual,"actual_fidelity":Ffinal,
         "certified_fidelity_lower_bound":Flow,"certified_target_weight_lower_bound":Ptarget,"rows":rows}
    Path(a.output).write_text(json.dumps(out,indent=2)); print(json.dumps({k:v for k,v in out.items() if k!="rows"},indent=2))
if __name__=="__main__": main()
