#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent

def save(fig,name):
 fig.tight_layout();fig.savefig(HERE/name,bbox_inches='tight');plt.close(fig)

def fig1():
 d=pd.read_csv(HERE/'source_data_method_ablation.csv');labels=['Projected\nsubspace','Bare same\nangles','Restricted\nSD','Restricted\nSDT','Restricted\nSDT + repeats','Generalized\n18q','Generalized\n24q'];fig,ax=plt.subplots(figsize=(9.2,5.4));ax.bar(np.arange(len(d)),d['energy_error_mEh']);ax.axhline(1.6,ls='--',label='1.6 mEh model-space benchmark');ax.set_yscale('log');ax.set_ylabel('Energy error relative to target doublet (mEh)');ax.set_xticks(np.arange(len(d)));ax.set_xticklabels(labels);ax.legend();save(fig,'fig1_method_ablation.pdf')
def fig2():
 d=pd.read_csv(HERE/'source_data_circuit_synthesis.csv');d=d[d.qubits==18];fig,ax=plt.subplots(figsize=(7.6,5.2));ax.plot(d.cx,abs(d.synthesis_error_mEh),marker='o',label='|Circuit - fermionic ansatz|');ax.plot(d.cx,d.total_error_mEh,marker='s',label='Circuit - exact doublet');ax.axhline(1.6,ls='--',label='1.6 mEh benchmark');ax.axhline(.1,ls=':',label='0.1 mEh synthesis criterion');ax.set_yscale('log');ax.set_xlabel('CX count');ax.set_ylabel('Energy error (mEh)');ax.legend();save(fig,'fig2_synthesis_tradeoff_18q.pdf')
def fig3():
 d=pd.read_csv(HERE/'source_data_measurement_shots.csv');fig,ax=plt.subplots(figsize=(7.4,5.2));ax.plot(d.sigma_mEh,d.shots_18q,marker='o',label='18 qubits');ax.plot(d.sigma_mEh,d.shots_24q,marker='s',label='24 qubits');ax.set_yscale('log');ax.invert_xaxis();ax.set_xlabel('Target statistical uncertainty (mEh)');ax.set_ylabel('Optimally allocated shots');ax.legend();save(fig,'fig3_measurement_shots.pdf')
def fig4():
 e18=np.array([0.00048357573,0.00047039804,0.00059436820,0.00026595395,0.00155171307,0.00033863674,0.00210936298,0.00085949634,0.00012553096]);e24=np.array([0.00019340111,0.00025390935,0.00048622471,0.00060305063,0.00160472119,0.00119589413,0.00086318520,0.00142126393,0.00068346131,0.00050670502,0.00080445950,0.00017160307]);fig,ax=plt.subplots(figsize=(7.5,5.2));ax.plot(np.arange(1,len(e18)+1),e18,marker='o',label='18 qubits');ax.plot(np.arange(1,len(e24)+1),e24,marker='s',label='24 qubits');ax.set_yscale('log');ax.set_xlabel('Natural-orbital index');ax.set_ylabel('Absolute occupation-number error');ax.legend();save(fig,'fig4_natural_occupation_errors.pdf')
def fig5():
 e18=[1.49266014,1.49265740,1.49269760,1.49266191,1.49266224,1.49269023,1.49305900];e24=[1.18070434,1.18075629,1.18095827,1.18119887,1.18124885];fig,ax=plt.subplots(figsize=(7.5,5.2));ax.scatter(np.arange(1,len(e18)+1),e18,label='18 qubits');ax.scatter(np.arange(1,len(e24)+1),e24,label='24 qubits');ax.axhline(1.6,ls='--',label='1.6 mEh benchmark');ax.set_xlabel('Perturbed-start reoptimization');ax.set_ylabel('Final energy error (mEh)');ax.legend();save(fig,'fig5_multistart.pdf')
if __name__=='__main__':
 fig1();fig2();fig3();fig4();fig5()
