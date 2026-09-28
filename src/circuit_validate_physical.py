#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,hashlib,numpy as np
from scipy.sparse.linalg import expm_multiply
from spin_fci import read_fcidump,pure_spin_hamiltonian,hf_determinant,spin_square_expectation
from physical_pool import entry_from_dict,fermionic_terms_for_entry

def sha256(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('result_json');ap.add_argument('--fcidump',required=True);ap.add_argument('--prefix',type=int,default=0);ap.add_argument('--method',choices=('lie','suzuki2'),default='suzuki2');ap.add_argument('--reps',type=int,default=1);ap.add_argument('--optimization-level',type=int,default=1,choices=(0,1,2,3));ap.add_argument('--basis-gates',default='rz,sx,x,cx');ap.add_argument('--seed-transpiler',type=int,default=None);ap.add_argument('--output',required=True);ap.add_argument('--strict',action='store_true');args=ap.parse_args()
 import qiskit,qiskit_aer,qiskit_nature
 from qiskit import QuantumCircuit,transpile
 from qiskit.circuit.library import PauliEvolutionGate
 from qiskit.synthesis import LieTrotter,SuzukiTrotter
 from qiskit_nature.second_q.mappers import JordanWignerMapper
 from qiskit_nature.second_q.operators import FermionicOp
 from qiskit_aer import AerSimulator
 r=json.load(open(args.result_json,encoding='utf-8'));d=read_fcidump(args.fcidump);na=(d.nelec+d.ms2)//2;nb=d.nelec-na;ms=.5*(na-nb)
 basis,H,U,Hd,Sp=pure_spin_hamiltonian(d,na,nb,target_s=.5);E0=float(np.linalg.eigvalsh(Hd)[0]);defs=r['selected_generators'];theta=np.asarray(r['theta'],float);nuse=len(defs) if args.prefix<=0 else min(args.prefix,len(defs));defs=defs[:nuse];theta=theta[:nuse];entries=[entry_from_dict(x,basis) for x in defs]
 pos={int(x):i for i,x in enumerate(basis)}[hf_determinant(d.norb,na,nb)];psi=np.zeros(len(basis),complex);psi[pos]=1
 for t,e in zip(theta,entries):psi=expm_multiply(float(t)*e.matrix,psi)
 psi/=np.linalg.norm(psi);Eexact=float(np.vdot(psi,H@psi).real);S2exact=spin_square_expectation(psi,Sp,ms=ms)
 nso=2*d.norb;hfd=hf_determinant(d.norb,na,nb);qc=QuantumCircuit(nso)
 for q in range(nso):
  if (hfd>>q)&1:qc.x(q)
 mapper=JordanWignerMapper();synth=LieTrotter(reps=args.reps) if args.method=='lie' else SuzukiTrotter(order=2,reps=args.reps);paulis=0
 for gd,t in zip(defs,theta):
  fop=FermionicOp(fermionic_terms_for_entry(gd),num_spin_orbitals=nso);qop=mapper.map(fop);Hh=((-1j)*qop).simplify(atol=1e-12);paulis+=len(Hh);qc.append(PauliEvolutionGate(Hh,time=-float(t),synthesis=synth),range(nso))
 bg=[x.strip() for x in args.basis_gates.split(',') if x.strip()];tqc=transpile(qc,basis_gates=bg,optimization_level=args.optimization_level,seed_transpiler=args.seed_transpiler);depth=int(tqc.depth());counts={str(k):int(v) for k,v in tqc.count_ops().items()};twoq=sum(v for k,v in counts.items() if k in ('cx','cz','ecr','iswap','rzz','rxx','ryy'))
 sim=AerSimulator(method='statevector');run=tqc.copy();run.save_statevector();sv=np.asarray(sim.run(run).result().get_statevector(run),complex);idx=np.asarray([int(x) for x in basis],dtype=np.int64);qms=sv[idx];pms=float(np.vdot(qms,qms).real)
 if pms>0:
  qcond=qms/np.sqrt(pms);Eq=float(np.vdot(qcond,H@qcond).real);S2q=spin_square_expectation(qcond,Sp,ms=ms);dc=U.conj().T@qcond;pD=float(np.vdot(dc,dc).real);Fcond=float(abs(np.vdot(psi,qcond))**2)
 else:Eq=S2q=pD=Fcond=float('nan')
 Ffull=float(abs(np.vdot(psi,qms))**2);probs=np.abs(sv)**2;pN=float(sum(probs[i] for i in range(len(probs)) if bin(int(i)).count('1')==d.nelec));de=(Eq-Eexact)*1e3
 passed=bool(Ffull>=.999 and abs(de)<=.1 and abs(S2q-.75)<=1e-3 and pms>=.999)
 out={'result_json':args.result_json,'result_json_sha256':sha256(args.result_json),'fcidump':args.fcidump,'fcidump_sha256':sha256(args.fcidump),'versions':{'qiskit':qiskit.__version__,'qiskit_aer':qiskit_aer.__version__,'qiskit_nature':qiskit_nature.__version__},'n_qubits':nso,'n_operators_used':nuse,'product_formula':args.method,'product_formula_reps':args.reps,'optimization_level':args.optimization_level,'seed_transpiler':args.seed_transpiler,'basis_gates':bg,'sum_pauli_terms_before_synthesis':paulis,'transpiled_depth':depth,'transpiled_gate_counts':counts,'two_qubit_gate_count_proxy':twoq,'exact_physical_prefix':{'energy_Eh':Eexact,'dev_vs_doublet_mEh':(Eexact-E0)*1e3,'S2':S2exact},'actual_transpiled_circuit':{'total_N_weight':pN,'target_Ms_weight':pms,'energy_Eh_conditioned_on_target_Ms':Eq,'energy_error_vs_exact_prefix_mEh':de,'S2_conditioned_on_target_Ms':S2q,'doublet_weight_within_target_Ms':pD,'fidelity_full_with_exact_physical_prefix':Ffull,'fidelity_conditioned_on_target_Ms':Fcond},'validation_criteria':{'fidelity_full_min':.999,'abs_energy_error_mEh_max':.1,'abs_S2_minus_0p75_max':1e-3,'target_Ms_weight_min':.999},'PASS':passed,'resource_note':'Abstract gate-set resources only; not device-mapped hardware resources.'}
 json.dump(out,open(args.output,'w'),indent=2);print(json.dumps(out,indent=2))
 if args.strict and not passed:raise SystemExit(2)
if __name__=='__main__':main()
