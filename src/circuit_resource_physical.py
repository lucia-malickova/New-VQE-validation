#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,hashlib,numpy as np
from spin_fci import read_fcidump,hf_determinant
from physical_pool import fermionic_terms_for_entry

def sha256(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('result_json');ap.add_argument('--fcidump',required=True);ap.add_argument('--prefix',type=int,default=0);ap.add_argument('--method',choices=('lie','suzuki2'),default='suzuki2');ap.add_argument('--reps',type=int,default=1);ap.add_argument('--optimization-level',type=int,default=1);ap.add_argument('--basis-gates',default='rz,sx,x,cx');ap.add_argument('--output',required=True);args=ap.parse_args()
 import qiskit,qiskit_nature
 from qiskit import QuantumCircuit,transpile
 from qiskit.circuit.library import PauliEvolutionGate
 from qiskit.synthesis import LieTrotter,SuzukiTrotter
 from qiskit_nature.second_q.mappers import JordanWignerMapper
 from qiskit_nature.second_q.operators import FermionicOp
 r=json.load(open(args.result_json));d=read_fcidump(args.fcidump);na=(d.nelec+d.ms2)//2;nb=d.nelec-na;nso=2*d.norb;hfd=hf_determinant(d.norb,na,nb);defs=r['selected_generators'];theta=np.asarray(r['theta'],float);nuse=len(defs) if args.prefix<=0 else min(args.prefix,len(defs));defs=defs[:nuse];theta=theta[:nuse]
 qc=QuantumCircuit(nso)
 for q in range(nso):
  if (hfd>>q)&1:qc.x(q)
 mapper=JordanWignerMapper();synth=LieTrotter(reps=args.reps) if args.method=='lie' else SuzukiTrotter(order=2,reps=args.reps);paulis=0
 for gd,t in zip(defs,theta):
  qop=mapper.map(FermionicOp(fermionic_terms_for_entry(gd),num_spin_orbitals=nso));Hh=((-1j)*qop).simplify(atol=1e-12);paulis+=len(Hh);qc.append(PauliEvolutionGate(Hh,time=-float(t),synthesis=synth),range(nso))
 bg=[x.strip() for x in args.basis_gates.split(',') if x.strip()];tqc=transpile(qc,basis_gates=bg,optimization_level=args.optimization_level);counts={str(k):int(v) for k,v in tqc.count_ops().items()};twoq=sum(v for k,v in counts.items() if k in ('cx','cz','ecr','iswap','rzz','rxx','ryy'))
 out={'result_json':args.result_json,'result_json_sha256':sha256(args.result_json),'fcidump':args.fcidump,'fcidump_sha256':sha256(args.fcidump),'qiskit_version':qiskit.__version__,'qiskit_nature_version':qiskit_nature.__version__,'n_qubits':nso,'n_operators_used':nuse,'product_formula':args.method,'product_formula_reps':args.reps,'optimization_level':args.optimization_level,'basis_gates':bg,'sum_pauli_terms_before_synthesis':paulis,'transpiled_depth':int(tqc.depth()),'gate_counts':counts,'two_qubit_gate_count_proxy':int(twoq),'warning':'Abstract gate-set transpilation only; not device-mapped hardware resources.'};json.dump(out,open(args.output,'w'),indent=2);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
