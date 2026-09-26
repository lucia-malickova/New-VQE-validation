#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import sys
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from novelty_upgrade_audit import (
    OUT as DEFAULT_OUT,
    projected_adapt,
    audit_projected_sequence,
    physical_adapt,
    plot_prefix,
    jdump,
)

def generate(system: str, outdir: Path):
    from pyscf import gto, scf, mcscf, ao2mo
    from pyscf.tools import fcidump

    if system == "NO":
        cfg = dict(
            atom="N 0 0 0; O 0 0 1.1508",
            geometry={"N":[0.0,0.0,0.0],"O":[0.0,0.0,1.1508]},
            charge=0, spin=1, ncas=6, nelecas=7, ncore=4,
            label="NO_STO3G_CAS7e6o", max_iters=120,
        )
    elif system == "OH":
        cfg = dict(
            atom="O 0 0 0; H 0 0 0.9700",
            geometry={"O":[0.0,0.0,0.0],"H":[0.0,0.0,0.9700]},
            charge=0, spin=1, ncas=5, nelecas=7, ncore=1,
            label="OH_STO3G_CAS7e5o", max_iters=100,
        )
    else:
        raise ValueError(system)

    mol = gto.M(
        atom=cfg["atom"], unit="Angstrom", basis="sto-3g",
        charge=cfg["charge"], spin=cfg["spin"], symmetry=False, verbose=0
    )
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-12
    mf.max_cycle = 200
    ehf = mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"{system} ROHF did not converge")

    expected_ncore = (mol.nelectron - cfg["nelecas"]) // 2
    if expected_ncore != cfg["ncore"]:
        raise RuntimeError((expected_ncore,cfg["ncore"]))

    mc = mcscf.CASCI(mf, cfg["ncas"], cfg["nelecas"])
    h1eff, ecore = mc.get_h1eff(mo_coeff=mf.mo_coeff)
    h2eff = mc.get_h2eff(mo_coeff=mf.mo_coeff)
    h2eff = ao2mo.restore(8, h2eff, cfg["ncas"])
    fpath = outdir / f"{cfg['label']}.FCIDUMP"
    fcidump.from_integrals(
        str(fpath), h1eff, h2eff, cfg["ncas"], cfg["nelecas"],
        nuc=float(ecore), ms=1, tol=1e-14
    )
    meta = {
        "system":system,
        "geometry_A":cfg["geometry"],
        "basis":"STO-3G",
        "reference":"ROHF doublet",
        "ROHF_energy_Eh":float(ehf),
        "full_electrons":int(mol.nelectron),
        "full_spatial_orbitals":int(mf.mo_coeff.shape[1]),
        "frozen_lowest_spatial_orbitals":int(cfg["ncore"]),
        "active_space":f"CAS({cfg['nelecas']}e,{cfg['ncas']}o)",
        "active_spin_orbitals":2*cfg["ncas"],
        "note":"Independent algorithmic open-shell benchmark; no converged chemical prediction is claimed.",
    }
    jdump(outdir/f"{cfg['label']}_metadata.json",meta)
    return fpath,cfg,meta

def add_controls(fcidump, projected, audit, tag, outdir):
    import numpy as np
    import scipy.linalg
    from scipy.sparse.linalg import expm_multiply

    sys.path.insert(0, str(ROOT/"src"))
    from spin_fci import (
        read_fcidump,pure_spin_hamiltonian,hf_determinant,
        excitation_specs_from_hf,excitation_generator,spin_square_expectation
    )

    d=read_fcidump(str(fcidump))
    na=(d.nelec+d.ms2)//2; nb=d.nelec-na; ms=.5*(na-nb)
    basis,H,U,Hd,Sp=pure_spin_hamiltonian(d,na,nb,target_s=.5)
    E0=float(np.linalg.eigvalsh(Hd)[0])
    hfd=hf_determinant(d.norb,na,nb); pos={int(x):i for i,x in enumerate(basis)}[hfd]
    psi0=np.zeros(len(basis),complex);psi0[pos]=1
    phi0=U.conj().T@psi0;phi0/=np.linalg.norm(phi0)
    specs=excitation_specs_from_hf(d.norb,na,nb,max_rank=int(projected["max_rank"]))
    sel=list(map(int,projected["selected_indices"])); th=np.asarray(projected["theta"],float)
    ratios=np.asarray(audit["leakage_ratio_F"],float)

    cache={}
    def A(k):
        if k not in cache: cache[k]=excitation_generator(basis,*specs[k]).tocsr()
        return cache[k]
    def Ad(k):
        x=U.conj().T@(A(k)@U)
        return .5*(x-x.conj().T)

    def states(scale=1.0, mask=None):
        pp=phi0.copy();bb=psi0.copy()
        for i,(t,k) in enumerate(zip(th,sel)):
            if mask is not None and not bool(mask[i]): continue
            pp=expm_multiply(scale*float(t)*Ad(k),pp)
            bb=expm_multiply(scale*float(t)*A(k),bb)
        pp/=np.linalg.norm(pp); bb/=np.linalg.norm(bb)
        proj=U@pp;proj/=np.linalg.norm(proj)
        c=U.conj().T@bb;pD=float(np.vdot(c,c).real)
        F=float(abs(np.vdot(proj,bb))**2)
        s2=spin_square_expectation(bb,Sp,ms=ms)
        Eb=float(np.vdot(bb,H@bb).real);Ep=float(np.vdot(proj,H@proj).real)
        postF=None
        if pD>1e-14:
            post=U@(c/math.sqrt(pD));post/=np.linalg.norm(post)
            postF=float(abs(np.vdot(proj,post))**2)
        return dict(scale=float(scale),fidelity=F,doublet_weight=pD,S2=s2,
                    projected_Eh=Ep,bare_Eh=Eb,bare_minus_projected_mEh=(Eb-Ep)*1e3,
                    postprojected_fidelity=postF)

    # Scale all amplitudes to test the small-amplitude behavior.
    scales=[0.0,0.01,0.02,0.05,0.1,0.2,0.4,0.6,0.8,1.0]
    scaled=[states(x) for x in scales]

    # Negative/positive controls on the same selected sequence.
    invariant=ratios<=1e-10
    leaking=~invariant
    controls={
        "n_invariant":int(np.sum(invariant)),
        "n_leaking":int(np.sum(leaking)),
        "invariant_only":states(1.0,invariant),
        "leaking_only":states(1.0,leaking),
    }

    # Estimate small-scale power laws from alpha <= 0.1, excluding zero.
    xs=[]; inf=[]; leakp=[]
    for r in scaled:
        if 0 < r["scale"] <= .1:
            xs.append(math.log(r["scale"]))
            inf.append(math.log(max(1e-300,1-r["fidelity"])))
            leakp.append(math.log(max(1e-300,1-r["doublet_weight"])))
    if len(xs)>=2:
        controls["small_scale_infidelity_power"]=float(np.polyfit(xs,inf,1)[0])
        controls["small_scale_leakage_probability_power"]=float(np.polyfit(xs,leakp,1)[0])

    out={"scaled_amplitudes":scaled,"controls":controls}
    jdump(outdir/f"{tag}_scaled_and_controls.json",out)

    import csv
    with open(outdir/f"{tag}_scaled_amplitudes.csv","w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=list(scaled[0].keys()));w.writeheader();w.writerows(scaled)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--system",choices=["NO","OH"],required=True)
    args=ap.parse_args()
    outdir=ROOT/"analysis"/"independent_benchmark_results"/args.system
    outdir.mkdir(parents=True,exist_ok=True)

    # Redirect helper outputs for this process.
    import novelty_upgrade_audit as nua
    nua.OUT=outdir

    fcidump,cfg,meta=generate(args.system,outdir)
    tag=args.system.lower()
    proj=projected_adapt(fcidump,tag,max_rank=3,max_iters=cfg["max_iters"],tol_meh=1.6)
    audit=audit_projected_sequence(fcidump,proj,tag)
    plot_prefix(tag,audit)
    controls=add_controls(fcidump,proj,audit,tag,outdir)
    phys=physical_adapt(fcidump,tag,max_iters=100,tol_meh=1.6)

    # Validate the independent physical ansatz at the qubit-circuit layer.
    circuit_results=[]
    result_json=outdir/f"{tag}_physical_adapt.json"
    for reps in (1,2,4):
        opath=outdir/f"{tag}_circuit_suzuki2_r{reps}.json"
        cmd=[
            sys.executable, str(ROOT/"src"/"circuit_validate_physical.py"),
            str(result_json), "--fcidump", str(fcidump),
            "--method","suzuki2","--reps",str(reps),
            "--optimization-level","1","--output",str(opath)
        ]
        subprocess.run(cmd,check=True,cwd=str(ROOT/"src"))
        circuit_results.append(json.loads(opath.read_text()))

    summary={
        "metadata":meta,
        "projected":{k:v for k,v in proj.items() if k not in ("selected_indices","theta","history")},
        "audit":{k:v for k,v in audit.items() if k not in ("prefix","lambda_op2","leakage_ratio_F")},
        "controls":controls["controls"],
        "physical":{k:v for k,v in phys.items() if k not in ("selected_pool_indices","selected_generators","theta","history")},
        "circuits":[{
            "reps":x["product_formula_reps"],
            "CX":x["transpiled_gate_counts"].get("cx",0),
            "depth":x["transpiled_depth"],
            "synthesis_error_mEh":x["actual_transpiled_circuit"]["energy_error_vs_exact_prefix_mEh"],
            "fidelity_full":x["actual_transpiled_circuit"]["fidelity_full_with_exact_physical_prefix"],
            "S2":x["actual_transpiled_circuit"]["S2_conditioned_on_target_Ms"],
            "Ms_weight":x["actual_transpiled_circuit"]["target_Ms_weight"],
            "PASS_strict":x["PASS"],
        } for x in circuit_results],
    }
    jdump(outdir/"SUMMARY.json",summary)
    print(json.dumps(summary,indent=2))

if __name__=="__main__":
    main()
