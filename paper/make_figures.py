#!/usr/bin/env python3
"""Regenerate manuscript and SI figures from committed CSV source data."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent

def save(fig,name):
    fig.tight_layout();fig.savefig(HERE/name,bbox_inches='tight');plt.close(fig)

def fig1():
 d=pd.read_csv(HERE/'source_data_method_ablation_v5.csv')
 labels=['Projected\nsubspace','Bare same\nangles','Invariant bare\n+ repeats','Restricted\nS+D','Restricted\nS+D+T','Restricted\nS+D+T + repeats','Generalized\n18q','Generalized\n24q']
 fig,ax=plt.subplots(figsize=(9.4,5.5));x=np.arange(len(d));ax.bar(x,d.energy_error_mEh);ax.axhline(1.6,linestyle='--',linewidth=1.2,label='1.6 mEh model-space benchmark');ax.axvline(6.5,linewidth=1.0);ax.text(7,max(d.energy_error_mEh)*.62,'expanded\n24q benchmark',ha='center',va='center',fontsize=9);ax.set_yscale('log');ax.set_ylabel('Energy error relative to target doublet (mEh)');ax.set_xticks(x);ax.set_xticklabels(labels);ax.set_title('Representation fidelity, spin symmetry, and expressivity are distinct requirements');ax.legend(fontsize=8);save(fig,'fig1_method_ablation_v5.pdf')

def fig2():
 d=pd.read_csv(HERE/'source_data_projected_negative_control_prefix_v5.csv')
 fig,ax=plt.subplots(figsize=(8.2,5.2));ax.plot(d.ops,d.projected_error_mEh,label='Projected ADAPT');ax.plot(d.ops,d.bare_error_vs_doublet_mEh,label='Same angles, bare sequence');ax.plot(d.ops,d.postprojected_doublet_error_mEh,label='Bare state projected back to doublet');ax.axhline(1.6,linestyle='--',linewidth=1.2,label='1.6 mEh model-space benchmark');ax.set_yscale('log');ax.set_xlabel('Number of selected operator applications');ax.set_ylabel('Energy error relative to exact target doublet (mEh)');ax.set_title('Projected optimization and physical evolution diverge along the ADAPT sequence');ax.legend(fontsize=8);save(fig,'fig2_projected_prefix_divergence.pdf')

def fig3():
 d=pd.read_csv(HERE/'source_data_circuit_synthesis_v5.csv');d=d[d.qubits==18]
 fig,ax=plt.subplots(figsize=(7.6,5.2));ax.plot(d.cx,abs(d.synthesis_error_mEh),marker='o',label='|Circuit - fermionic ansatz|');ax.plot(d.cx,d.total_error_mEh,marker='s',label='Circuit - exact doublet');ax.axhline(1.6,linestyle='--',linewidth=1.2,label='1.6 mEh model-space benchmark');ax.axhline(.1,linestyle=':',linewidth=1.2,label='0.1 mEh strict synthesis criterion');ax.set_yscale('log');ax.set_xlabel('CX count after abstract-basis transpilation');ax.set_ylabel('Energy error (mEh)');
 for _,r in d.iterrows():ax.annotate(f'r={int(r.suzuki_reps)}',(r.cx,r.total_error_mEh),xytext=(4,5),textcoords='offset points')
 ax.set_title('18-qubit Suzuki synthesis: fidelity is bought with gates');ax.legend(fontsize=8);save(fig,'fig3_synthesis_tradeoff_18q.pdf')

def figS1():
 d=pd.read_csv(HERE/'source_data_natural_occupations_v5.csv');fig,ax=plt.subplots(figsize=(7.5,5.2))
 for case,marker in [('18q','o'),('24q','s')]:
  x=d[d.case==case];ax.plot(x.orbital_index,x.absolute_error,marker=marker,label=case.replace('q',' qubits'))
 ax.set_yscale('log');ax.set_xlabel('Natural-orbital index');ax.set_ylabel('Absolute occupation-number error');ax.set_title('One-particle occupations remain close to exact doublet values');ax.legend();save(fig,'figS1_natural_occupation_errors.pdf')

def figS2():
 d=pd.read_csv(HERE/'source_data_multistart_v5.csv');fig,ax=plt.subplots(figsize=(7.5,5.2))
 for case,marker in [('18q','o'),('24q','s')]:
  x=d[d.case==case];ax.scatter(np.arange(1,len(x)+1),x.dev_mEh,marker=marker,label=case.replace('q',' qubits'))
 ax.axhline(1.6,linestyle='--',linewidth=1.2,label='1.6 mEh model-space benchmark');ax.set_xlabel('Perturbed-start reoptimization');ax.set_ylabel('Final energy error (mEh)');ax.set_title('Compact ansatz solutions are locally robust');ax.legend();save(fig,'figS2_multistart.pdf')

def figS3():
 d=pd.read_csv(HERE/'source_data_representation_scaling.csv');fig,ax=plt.subplots(figsize=(7.2,5.2));ax.loglog(d.theta,d.leak_amp,marker='o',label=r'Leakage amplitude $\|Qe^{\theta A}|\Phi\rangle\|$');ax.loglog(d.theta,d.compressed_mismatch_norm,marker='s',label='Compressed-evolution mismatch');m=d.postprojected_infidelity>0;ax.loglog(d.loc[m,'theta'],d.loc[m,'postprojected_infidelity'],marker='^',label='Post-projection infidelity');ax.set_xlabel(r'Amplitude $\theta$');ax.set_ylabel('Norm / infidelity');ax.set_title('Single-generator representation mismatch in the Cu active space');ax.legend(fontsize=8);save(fig,'figS3_representation_scaling.pdf')

def figS4():
 d=pd.read_csv(HERE/'source_data_projected_negative_control_prefix_v5.csv');fig,ax=plt.subplots(figsize=(8.0,5.0));ax.plot(d.ops,d.full_fidelity,label='Fidelity: projected vs bare');ax.plot(d.ops,d.doublet_weight,label='Bare-state doublet weight');ax.plot(d.ops,d.postprojected_doublet_fidelity,label='Fidelity after post-projection');ax.set_xlabel('Number of selected operator applications');ax.set_ylabel('State overlap / weight');ax.set_ylim(0,1.03);ax.set_title('State diagnostics along the canonical projected-ADAPT negative control');ax.legend(fontsize=8);save(fig,'figS4_projected_prefix_state_diagnostics.pdf')

def figS5():
 d=pd.read_csv(HERE/'source_data_hubbard4_representation_v5.csv');fig,ax=plt.subplots(figsize=(7.4,5.0));ax.plot(d.U_over_t,d.projected_bare_fidelity,marker='o',label='Projected vs bare fidelity');ax.plot(d.U_over_t,d.projected_postprojected_fidelity,marker='s',label='Projected vs post-projected bare fidelity');ax.set_xlabel(r'$U/t$');ax.set_ylabel('Fidelity');ax.set_title('Independent Hubbard-model stress test of representation equivalence');ax.legend();save(fig,'figS5_hubbard_representation_audit.pdf')

def figS6():
    d=pd.read_csv(HERE/"source_data_local_suzuki_predictor_v5.csv")
    fig,ax=plt.subplots(figsize=(7.4,5.2))
    for case,marker in [("18q","o"),("24q","s"),("NO","^")]:
        q=d[(d["case"]==case)&(d["comm_relF"]>0)&(d["local_infidelity"]>0)]
        ax.loglog(q["theta_abs"],q["local_infidelity"],linestyle="None",
                  marker=marker,label=case)
    ax.set_xlabel(r"Optimized amplitude magnitude $|\\theta|$")
    ax.set_ylabel("Single-application Suzuki-2 infidelity")
    ax.set_title("Among noncommuting factors, amplitude dominates local Suzuki difficulty")
    ax.legend()
    save(fig,"figS6_local_suzuki_amplitude_v6.pdf")

def figS7():
 d=pd.read_csv(HERE/'source_data_measurement_shots_v5.csv');fig,ax=plt.subplots(figsize=(7.4,5.2));ax.plot(d.sigma_mEh,d.shots_18q,marker='o',label='18 qubits');ax.plot(d.sigma_mEh,d.shots_24q,marker='s',label='24 qubits');ax.set_yscale('log');ax.invert_xaxis();ax.set_xlabel('Target statistical uncertainty (mEh)');ax.set_ylabel('Optimally allocated shots for greedy QWC groups');ax.set_title('Measurement burden at the mEh scale');ax.legend();save(fig,'figS7_measurement_shots.pdf')

if __name__=='__main__':
 for f in (fig1,fig2,fig3,figS1,figS2,figS3,figS4,figS5,figS6,figS7):f()
 print('Wrote manuscript/SI figures to',HERE)
