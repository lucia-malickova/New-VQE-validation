#!/usr/bin/env python3
"""First-order Hartree-Fock target-space tangent-accessibility audit."""
import argparse,json,sys
from pathlib import Path
import numpy as np
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("fcidump"); ap.add_argument("--src",default="src")
    ap.add_argument("--output",default="HF_tangent_accessibility.json"); a=ap.parse_args(); sys.path.insert(0,a.src)
    from spin_fci import read_fcidump,pure_spin_hamiltonian,hf_determinant,excitation_specs_from_hf,excitation_generator
    from physical_pool import build_pool
    D=read_fcidump(a.fcidump); na=(D.nelec+D.ms2)//2; nb=D.nelec-na
    basis,H,U,Hd,Sp=pure_spin_hamiltonian(D,na,nb,.5); pos={int(d):i for i,d in enumerate(basis)}
    hf=np.zeros(len(basis)); hf[pos[hf_determinant(D.norb,na,nb)]]=1.
    rank=lambda v:0 if not v else int(np.linalg.matrix_rank(np.column_stack(v),tol=1e-10))
    specs=excitation_specs_from_hf(D.norb,na,nb,max_rank=3); proj=[]; inv=[]; n_inv=0
    for s in specs:
        A=excitation_generator(basis,*s); AU=A@U; inside=U@(U.T@AU)
        ell=float(np.linalg.norm(AU-inside))/max(float(np.linalg.norm(AU)),1e-300)
        v=A@hf; vp=U@(U.T@v)
        if np.linalg.norm(vp)>1e-12: proj.append(vp)
        if ell<1e-10:
            n_inv+=1
            if np.linalg.norm(v)>1e-12: inv.append(v)
    pool=build_pool(basis,U,Sp,D.norb,na,nb,max_rank=2); gen=[]
    for e in pool:
        v=e.matrix@hf
        if np.linalg.norm(v)>1e-12: gen.append(v)
    out={"target_dimension":int(U.shape[1]),"maximum_tangent_dimension":int(U.shape[1]-1),
         "projected_bare_pool_size":len(specs),"projected_bare_HF_tangent_rank":rank(proj),
         "invariant_bare_pool_size":n_inv,"invariant_bare_HF_tangent_rank":rank(inv),
         "generalized_physical_pool_size":len(pool),"generalized_physical_nonzero_on_HF":len(gen),
         "generalized_physical_HF_tangent_rank":rank(gen)}
    Path(a.output).write_text(json.dumps(out,indent=2)); print(json.dumps(out,indent=2))
if __name__=="__main__": main()
