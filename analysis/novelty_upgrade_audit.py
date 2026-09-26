#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg
import scipy.sparse as sp
from scipy.optimize import minimize
from scipy.sparse.linalg import expm_multiply

ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "legacy" / "initial_repository_snapshot"
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(LEGACY))

from spin_fci import (
    read_fcidump,
    pure_spin_hamiltonian,
    sector_hamiltonian,
    hf_determinant,
    excitation_specs_from_hf,
    excitation_generator,
    spin_square_expectation,
)
from physical_pool import build_pool

OUT = ROOT / "analysis" / "novelty_upgrade_results"
OUT.mkdir(parents=True, exist_ok=True)


def jdump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2), encoding="utf-8")


def state_setup(data):
    nalpha = (data.nelec + data.ms2) // 2
    nbeta = data.nelec - nalpha
    ms = 0.5 * (nalpha - nbeta)
    if abs(ms - 0.5) > 1e-12:
        raise RuntimeError(f"Expected M_S=1/2, got {ms}")
    basis, Hms, U, Hd, Splus = pure_spin_hamiltonian(data, nalpha, nbeta, target_s=0.5)
    hfd = hf_determinant(data.norb, nalpha, nbeta)
    pos = {int(d): i for i, d in enumerate(basis)}[hfd]
    psi0 = np.zeros(len(basis), dtype=complex)
    psi0[pos] = 1.0
    phi0 = U.conj().T @ psi0
    pn = float(np.vdot(phi0, phi0).real)
    if abs(pn - 1.0) > 1e-9:
        raise RuntimeError(f"HF reference not pure doublet: projected norm={pn}")
    phi0 /= math.sqrt(pn)
    return nalpha, nbeta, ms, basis, Hms, U, Hd, Splus, psi0, phi0


def projected_adapt(fcidump, tag, max_rank=3, max_iters=160, tol_meh=1.6):
    data = read_fcidump(str(fcidump))
    nalpha, nbeta, ms, basis, Hms, U, Hd, Splus, psi0, phi0 = state_setup(data)
    E0 = float(np.linalg.eigvalsh(Hd)[0])
    Eq = float(np.linalg.eigvalsh(sector_hamiltonian(data, nalpha + 1, nbeta - 1)[1])[0])
    specs = excitation_specs_from_hf(data.norb, nalpha, nbeta, max_rank=max_rank)

    Ams_cache = {}
    Ad_cache = {}
    def Ams(k):
        if k not in Ams_cache:
            Ams_cache[k] = excitation_generator(basis, *specs[k])
        return Ams_cache[k]
    def Ad(k):
        if k not in Ad_cache:
            M = U.conj().T @ (Ams(k) @ U)
            Ad_cache[k] = 0.5 * (M - M.conj().T)
        return Ad_cache[k]

    selected = []
    theta = np.zeros(0, dtype=float)
    hist = []

    def build(th):
        v = phi0.copy()
        for t, k in zip(th, selected):
            v = expm_multiply(float(t) * Ad(k), v)
        return v

    def obj_grad(th):
        states = [phi0]
        for t, k in zip(th, selected):
            states.append(expm_multiply(float(t) * Ad(k), states[-1]))
        psi = states[-1]
        w = Hd @ psi
        E = float(np.vdot(psi, w).real)
        g = np.zeros(len(th), dtype=float)
        for i in range(len(th)-1, -1, -1):
            Ai = Ad(selected[i])
            g[i] = 2.0 * np.real(np.vdot(w, Ai @ states[i+1]))
            w = expm_multiply(-float(th[i]) * Ai, w)
        return E, g

    t0 = time.time()
    for it in range(1, max_iters + 1):
        phi = build(theta)
        psi_ms = U @ phi
        u_ms = Hms @ psi_ms
        grads = np.full(len(specs), -np.inf, dtype=float)
        used = set(selected)
        for k in range(len(specs)):
            if k in used:
                continue
            grads[k] = abs(2.0 * np.real(np.vdot(u_ms, Ams(k) @ psi_ms)))
        kbest = int(np.argmax(grads))
        gmax = float(grads[kbest])
        if not np.isfinite(gmax) or gmax < 1e-8:
            break
        selected.append(kbest)
        theta = np.append(theta, 0.0)
        res = minimize(
            obj_grad, theta, method="L-BFGS-B", jac=True,
            options={"gtol":1e-9, "ftol":1e-13, "maxiter":700, "maxls":50}
        )
        theta = np.asarray(res.x, dtype=float)
        E = float(res.fun)
        dev = (E - E0) * 1e3
        phi = build(theta)
        psi_ms = U @ phi
        s2 = spin_square_expectation(psi_ms, Splus, ms=ms)
        resid = float(np.linalg.norm(Hd @ phi - E * phi))
        hist.append([len(selected), E, dev, gmax, s2, resid, bool(res.success)])
        print(f"[{tag}] projected iter={it:3d} ops={len(selected):3d} dev={dev:10.6f} mEh g={gmax:.3e} s2={s2:.12f}", flush=True)
        if dev <= tol_meh:
            break

    result = {
        "tag": tag,
        "fcidump": str(fcidump),
        "norb": data.norb,
        "nelec": data.nelec,
        "nalpha": nalpha,
        "nbeta": nbeta,
        "Ms_dim": len(basis),
        "doublet_dim": U.shape[1],
        "pool_size": len(specs),
        "max_rank": max_rank,
        "E_doublet": E0,
        "E_quartet": Eq,
        "selected_indices": [int(x) for x in selected],
        "theta": [float(x) for x in theta],
        "history": hist,
        "elapsed_s": time.time()-t0,
    }
    jdump(OUT / f"{tag}_projected_adapt.json", result)
    return result


def audit_projected_sequence(fcidump, projected, tag):
    data = read_fcidump(str(fcidump))
    nalpha, nbeta, ms, basis, Hms, U, Hd, Splus, psi0, phi0 = state_setup(data)
    E0 = float(np.linalg.eigvalsh(Hd)[0])
    specs = excitation_specs_from_hf(data.norb, nalpha, nbeta, max_rank=int(projected["max_rank"]))
    selected = list(map(int, projected["selected_indices"]))
    theta = np.asarray(projected["theta"], dtype=float)

    A_cache = {}
    Ad_cache = {}
    def A(k):
        if k not in A_cache:
            A_cache[k] = excitation_generator(basis, *specs[k]).tocsr()
        return A_cache[k]
    def Ad(k):
        if k not in Ad_cache:
            X = U.conj().T @ (A(k) @ U)
            Ad_cache[k] = 0.5*(X-X.conj().T)
        return Ad_cache[k]

    lambdas = []
    fro_ratios = []
    for k in selected:
        AU = A(k) @ U
        inside = U @ (U.conj().T @ AU)
        leak = AU - inside
        denomF = float(np.linalg.norm(AU))
        fro_ratios.append(0.0 if denomF == 0 else float(np.linalg.norm(leak))/denomF)
        # spectral norm of Q A P represented as leak acting on the orthonormal U coordinates
        lambdas.append(float(scipy.linalg.svdvals(np.asarray(leak))[0]) if leak.size else 0.0)

    phi = phi0.copy()
    bare = psi0.copy()
    L = 0.0
    rows = []
    for i, (t, k, lam, fr) in enumerate(zip(theta, selected, lambdas, fro_ratios), start=1):
        phi = expm_multiply(float(t)*Ad(k), phi)
        bare = expm_multiply(float(t)*A(k), bare)
        phi /= np.linalg.norm(phi)
        bare /= np.linalg.norm(bare)
        proj = U @ phi
        proj /= np.linalg.norm(proj)
        cD = U.conj().T @ bare
        pD = float(np.vdot(cD,cD).real)
        pD = min(1.0,max(0.0,pD))
        F = float(abs(np.vdot(proj,bare))**2)
        Fpost = None
        if pD > 1e-14:
            post = U @ (cD/math.sqrt(pD))
            post /= np.linalg.norm(post)
            Fpost = float(abs(np.vdot(proj,post))**2)
        Eproj = float(np.vdot(proj,Hms@proj).real)
        Ebare = float(np.vdot(bare,Hms@bare).real)
        s2b = spin_square_expectation(bare,Splus,ms=ms)
        L += abs(float(t))*float(lam)
        rows.append({
            "prefix":i, "pool_index":int(k), "theta":float(t),
            "lambda_op2":float(lam), "leakage_ratio_F":float(fr),
            "L_cumulative":float(L), "doublet_weight_bare":pD,
            "state_leakage_norm":math.sqrt(max(0.0,1.0-pD)),
            "fidelity_projected_bare":F,
            "fidelity_projected_postprojected_bare":Fpost,
            "E_projected_Eh":Eproj, "E_bare_Eh":Ebare,
            "bare_minus_projected_mEh":(Ebare-Eproj)*1e3,
            "S2_bare":s2b,
        })

    proj = U @ phi
    proj /= np.linalg.norm(proj)
    Eproj = float(np.vdot(proj,Hms@proj).real)
    Ebare = float(np.vdot(bare,Hms@bare).real)
    cD = U.conj().T @ bare
    pD = float(np.vdot(cD,cD).real)
    postE = None
    postF = None
    if pD > 1e-14:
        post = U @ (cD/math.sqrt(pD))
        post /= np.linalg.norm(post)
        postE = float(np.vdot(post,Hms@post).real)
        postF = float(abs(np.vdot(proj,post))**2)

    out = {
        "tag":tag,
        "n_ops":len(selected),
        "E_projected_Eh":Eproj,
        "projected_dev_mEh":(Eproj-E0)*1e3,
        "E_bare_Eh":Ebare,
        "bare_dev_mEh":(Ebare-E0)*1e3,
        "fidelity_projected_bare":float(abs(np.vdot(proj,bare))**2),
        "doublet_weight_bare":pD,
        "S2_bare":spin_square_expectation(bare,Splus,ms=ms),
        "postprojected_bare_Eh":postE,
        "postprojected_bare_dev_mEh":None if postE is None else (postE-E0)*1e3,
        "fidelity_projected_postprojected_bare":postF,
        "lambda_op2":lambdas,
        "leakage_ratio_F":fro_ratios,
        "L_total":float(L),
        "n_noninvariant_gt_1e10":int(np.sum(np.asarray(fro_ratios)>1e-10)),
        "prefix":rows,
    }
    jdump(OUT/f"{tag}_representation_audit.json",out)

    # CSV for figure/source data
    import csv
    with open(OUT/f"{tag}_prefix_audit.csv","w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    return out


def generate_no_active_fcidump():
    from pyscf import gto, scf, mcscf, ao2mo
    from pyscf.tools import fcidump

    # Independent open-shell molecular benchmark. Geometry is used only to define
    # a reproducible algorithmic Hamiltonian; no spectroscopic claim is made.
    mol = gto.M(
        atom="N 0 0 0; O 0 0 1.1508",
        unit="Angstrom",
        basis="sto-3g",
        charge=0,
        spin=1,  # 2S = 1, doublet
        symmetry=False,
        verbose=0,
    )
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-12
    mf.max_cycle = 200
    ehf = mf.kernel()
    if not mf.converged:
        raise RuntimeError("NO ROHF did not converge")

    ncas = 6
    nelecas = 7
    mc = mcscf.CASCI(mf, ncas, nelecas)
    # Default CASCI orbital ordering retains the last ncas orbitals active,
    # freezing four lowest doubly occupied orbitals for 15e/10o STO-3G NO.
    h1eff, ecore = mc.get_h1eff(mo_coeff=mf.mo_coeff)
    h2eff = mc.get_h2eff(mo_coeff=mf.mo_coeff)
    h2eff = ao2mo.restore(8, h2eff, ncas)
    path = OUT/"NO_STO3G_CAS7e6o.FCIDUMP"
    fcidump.from_integrals(
        str(path), h1eff, h2eff, ncas, nelecas, nuc=float(ecore), ms=1, tol=1e-14
    )
    meta = {
        "molecule":"NO",
        "geometry_A":{"N":[0,0,0],"O":[0,0,1.1508]},
        "basis":"STO-3G",
        "reference":"ROHF doublet",
        "ROHF_energy_Eh":float(ehf),
        "full_spatial_orbitals":int(mf.mo_coeff.shape[1]),
        "frozen_lowest_spatial_orbitals":4,
        "active_space":"CAS(7e,6o)",
        "active_spin_orbitals":12,
        "note":"Algorithmic independent open-shell benchmark; not a converged chemistry claim."
    }
    jdump(OUT/"NO_STO3G_CAS7e6o_metadata.json",meta)
    return path, meta


def physical_adapt(fcidump, tag, max_iters=80, tol_meh=1.6):
    data = read_fcidump(str(fcidump))
    nalpha, nbeta, ms, basis, Hms, U, Hd, Splus, psi0, phi0 = state_setup(data)
    E0 = float(np.linalg.eigvalsh(Hd)[0])
    pool = build_pool(basis,U,Splus,data.norb,nalpha,nbeta,max_rank=2)
    print(f"[{tag}] physical pool size={len(pool)}",flush=True)

    selected = []
    theta = np.zeros(0)
    hist = []

    def B(idx):
        return pool[idx].matrix

    def build(th):
        v = psi0.copy()
        for t,idx in zip(th,selected):
            v = expm_multiply(float(t)*B(idx),v)
        return v

    def obj_grad(th):
        states=[psi0]
        for t,idx in zip(th,selected):
            states.append(expm_multiply(float(t)*B(idx),states[-1]))
        psi=states[-1]
        w=Hms@psi
        E=float(np.vdot(psi,w).real)
        g=np.zeros(len(th))
        for i in range(len(th)-1,-1,-1):
            Bi=B(selected[i])
            g[i]=2.0*np.real(np.vdot(w,Bi@states[i+1]))
            w=expm_multiply(-float(th[i])*Bi,w)
        return E,g

    t0=time.time()
    for it in range(1,max_iters+1):
        psi=build(theta)
        w=Hms@psi
        grads=np.empty(len(pool))
        for j,e in enumerate(pool):
            grads[j]=abs(2.0*np.real(np.vdot(w,e.matrix@psi)))
        jbest=int(np.argmax(grads))
        gmax=float(grads[jbest])
        if gmax<1e-8:
            break
        selected.append(jbest)
        theta=np.append(theta,0.0)
        res=minimize(obj_grad,theta,method="L-BFGS-B",jac=True,
                     options={"gtol":1e-9,"ftol":1e-13,"maxiter":700,"maxls":50})
        theta=np.asarray(res.x,float)
        E=float(res.fun)
        dev=(E-E0)*1e3
        psi=build(theta); psi/=np.linalg.norm(psi)
        s2=spin_square_expectation(psi,Splus,ms=ms)
        pD=float(np.linalg.norm(U.conj().T@psi)**2)
        resid=float(np.linalg.norm(Hms@psi-E*psi))
        hist.append([len(selected),E,dev,gmax,s2,pD,resid,bool(res.success)])
        print(f"[{tag}] physical iter={it:3d} ops={len(selected):3d} dev={dev:10.6f} mEh g={gmax:.3e} s2={s2:.12f}",flush=True)
        if dev<=tol_meh:
            break

    psi=build(theta); psi/=np.linalg.norm(psi)
    out={
        "tag":tag,
        "pool_size":len(pool),
        "n_operator_applications":len(selected),
        "n_unique_generators":len(set(selected)),
        "selected_pool_indices":[int(x) for x in selected],
        "selected_generators":[pool[i].as_dict() for i in selected],
        "theta":[float(x) for x in theta],
        "E_doublet":E0,
        "energy_Eh":float(np.vdot(psi,Hms@psi).real),
        "dev_mEh":float((np.vdot(psi,Hms@psi).real-E0)*1e3),
        "S2":spin_square_expectation(psi,Splus,ms=ms),
        "doublet_weight":float(np.linalg.norm(U.conj().T@psi)**2),
        "fidelity_exact_doublet":None,
        "residual_Eh":float(np.linalg.norm(Hms@psi-float(np.vdot(psi,Hms@psi).real)*psi)),
        "max_selected_leakage_ratio":max((pool[i].leakage_ratio for i in selected),default=0.0),
        "max_selected_commutator_ratio":max((pool[i].commutator_ratio for i in selected),default=0.0),
        "history":hist,
        "elapsed_s":time.time()-t0,
    }
    # exact doublet fidelity
    evals,evecs=np.linalg.eigh(Hd)
    exact=U@evecs[:,0]
    exact/=np.linalg.norm(exact)
    out["fidelity_exact_doublet"]=float(abs(np.vdot(exact,psi))**2)
    jdump(OUT/f"{tag}_physical_adapt.json",out)
    return out


def plot_prefix(tag,audit):
    import matplotlib.pyplot as plt
    pref=np.array([r["prefix"] for r in audit["prefix"]])
    F=np.array([r["fidelity_projected_bare"] for r in audit["prefix"]])
    pD=np.array([r["doublet_weight_bare"] for r in audit["prefix"]])
    L=np.array([r["L_cumulative"] for r in audit["prefix"]])
    fig,ax=plt.subplots(figsize=(7.2,4.8))
    ax.plot(pref,F,label="Projected--bare fidelity")
    ax.plot(pref,pD,label="Bare doublet weight")
    ax.set_xlabel("Prefix length")
    ax.set_ylabel("State metric")
    ax.set_ylim(-0.02,1.02)
    ax2=ax.twinx()
    ax2.plot(pref,L,linestyle="--",label="Cumulative bound budget L")
    ax2.set_ylabel("L")
    lines,labels=ax.get_legend_handles_labels()
    lines2,labels2=ax2.get_legend_handles_labels()
    ax.legend(lines+lines2,labels+labels2,loc="best")
    fig.tight_layout()
    fig.savefig(OUT/f"{tag}_prefix_audit.pdf")
    fig.savefig(OUT/f"{tag}_prefix_audit.png",dpi=180)
    plt.close(fig)


def main():
    summary={"started":time.time()}

    # 1. Recover the canonical 123-operator projected Cu calculation deterministically.
    cu_json=OUT/"cu18_projected_rerun.json"
    cu_npz=OUT/"cu18_projected_rerun.npz"
    cmd=[
        sys.executable, str(LEGACY/"run_adapt_doublet.py"),
        str(ROOT/"data/18q/active.FCIDUMP"),
        "--max-iters","130","--max-rank","3","--energy-tol-meh","1.6",
        "--output",str(cu_npz),"--json",str(cu_json)
    ]
    print("RUNNING CU PROJECTED RERUN",flush=True)
    subprocess.run(cmd,check=True,cwd=str(ROOT))
    cu_proj=json.loads(cu_json.read_text())
    cu_projected={
        "max_rank":3,
        "selected_indices":cu_proj["selected_indices"],
        "theta":cu_proj["theta"],
    }
    cu_audit=audit_projected_sequence(ROOT/"data/18q/active.FCIDUMP",cu_projected,"cu18_rerun")
    plot_prefix("cu18_rerun",cu_audit)

    canonical=json.loads((ROOT/"results/final/projected_vs_physical_summary.json").read_text())
    checks={
        "n_ops_123":len(cu_proj["selected_indices"])==123,
        "projected_energy_match":abs(cu_audit["E_projected_Eh"]-canonical["projected_ADAPT_state"]["energy_Eh"])<1e-8,
        "bare_energy_match":abs(cu_audit["E_bare_Eh"]-canonical["same_parameters_unprojected_generators"]["energy_Eh"])<1e-8,
        "bare_fidelity_match":abs(cu_audit["fidelity_projected_bare"]-canonical["same_parameters_unprojected_generators"]["fidelity_with_projected_state"])<1e-8,
        "doublet_weight_match":abs(cu_audit["doublet_weight_bare"]-canonical["same_parameters_unprojected_generators"]["doublet_weight"])<1e-8,
    }
    checks["canonical_reproduction_pass"]=all(checks.values())
    summary["cu18_reproduction_checks"]=checks
    summary["cu18_audit"]= {k:v for k,v in cu_audit.items() if k not in ("prefix","lambda_op2","leakage_ratio_F")}

    # 2. Independent NO/STO-3G CAS(7e,6o) benchmark.
    no_fc,no_meta=generate_no_active_fcidump()
    no_proj=projected_adapt(no_fc,"no12",max_rank=3,max_iters=120,tol_meh=1.6)
    no_audit=audit_projected_sequence(no_fc,no_proj,"no12")
    plot_prefix("no12",no_audit)
    no_phys=physical_adapt(no_fc,"no12",max_iters=80,tol_meh=1.6)
    summary["NO_metadata"]=no_meta
    summary["NO_projected"]={k:v for k,v in no_proj.items() if k not in ("theta","selected_indices","history")}
    summary["NO_audit"]={k:v for k,v in no_audit.items() if k not in ("prefix","lambda_op2","leakage_ratio_F")}
    summary["NO_physical"]={k:v for k,v in no_phys.items() if k not in ("theta","selected_pool_indices","history")}

    summary["elapsed_s"]=time.time()-summary["started"]
    jdump(OUT/"SUMMARY.json",summary)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":
    main()

# CI trigger: 2026-09-26 novelty audit
