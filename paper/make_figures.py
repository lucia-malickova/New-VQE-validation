#!/usr/bin/env python3
"""Regenerate manuscript and SI figures from committed CSV source data."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = Path(__file__).resolve().parent

mpl.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'legend.fontsize': 10,
    'xtick.labelsize': 10.5,
    'ytick.labelsize': 10.5,
    'figure.dpi': 200,
    'savefig.dpi': 300,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.9,
})

PALETTE = {
    'blue': '#1f77b4',
    'orange': '#d95f02',
    'green': '#2ca25f',
    'purple': '#756bb1',
    'gray': '#6b7280',
    'lightblue': '#E8F1FA',
    'lightred': '#FBE9E0',
    'lightgreen': '#E8F6EE',
    'lightpurple': '#EEEAF9',
}

LW = 2.4
MS = 4.5


def stylize(ax, ygrid=True):
    if ygrid:
        ax.grid(True, axis='y', which='major', color='#D8DDE6', linewidth=0.8)
        ax.grid(True, axis='y', which='minor', color='#EEF1F4', linewidth=0.5)
    ax.tick_params(direction='out', length=3.2, width=0.8)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(HERE / name, bbox_inches='tight')
    plt.close(fig)


def fig1():
    fig, ax = plt.subplots(figsize=(11.0, 3.25))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    boxes = [
        (0.02, 0.28, 0.205, 0.50, PALETTE['lightblue'], PALETTE['blue'],
         'Projected / reduced-space\noptimization', r'$\{U^\dagger A_kU,\,\theta_k\}$'),
        (0.285, 0.28, 0.205, 0.50, PALETTE['lightred'], PALETTE['orange'],
         'Parent full-space\nfermionic ansatz', r'$\prod_k e^{\theta_k A_k}$'),
        (0.55, 0.28, 0.18, 0.50, PALETTE['lightgreen'], PALETTE['green'],
         'Synthesized\nqubit circuit', 'Suzuki / factorization'),
        (0.79, 0.28, 0.19, 0.50, PALETTE['lightpurple'], PALETTE['purple'],
         'Measurement /\nestimator', 'Sampling / statistics'),
    ]
    for x, y, w, h, fc, ec, title, subtitle in boxes:
        patch = FancyBboxPatch((x, y), w, h,
                               boxstyle='round,pad=0.012,rounding_size=0.03',
                               linewidth=1.5, edgecolor=ec, facecolor=fc)
        ax.add_patch(patch)
        ax.text(x + w/2, y + h*0.65, title, ha='center', va='center',
                fontsize=10.8, weight='bold')
        ax.text(x + w/2, y + h*0.28, subtitle, ha='center', va='center', fontsize=9.8)

    links = [
        ((0.228, 0.53), (0.282, 0.53), PALETTE['orange'],
         r'Representation check\n$QA_kP=0$'),
        ((0.493, 0.53), (0.547, 0.53), PALETTE['gray'],
         'Synthesis check'),
        ((0.733, 0.53), (0.787, 0.53), PALETTE['gray'],
         'Estimation layer'),
    ]
    for (x1, y1), (x2, y2), c, label in links:
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',
                                    mutation_scale=14, lw=1.7, color=c))
        ax.text((x1+x2)/2, 0.84, label, ha='center', va='bottom',
                fontsize=9.1, color=c)

    ax.text(0.39, 0.08,
            'If the representation check fails, a faithful downstream circuit can implement a different variational state from the projected optimum.',
            ha='center', va='center', fontsize=9.2, color='#333333', wrap=True)
    save(fig, 'fig1_validation_hierarchy.pdf')


def fig2():
    d = pd.read_csv(HERE / 'source_data_projected_negative_control_prefix_v5.csv')
    fig, ax = plt.subplots(figsize=(8.2, 5.3))
    ax.plot(d.ops, d.projected_error_mEh, label='Projected final-angle prefix', color=PALETTE['blue'], lw=LW)
    ax.plot(d.ops, d.bare_error_vs_doublet_mEh, label='Parent full-space final-angle prefix', color=PALETTE['orange'], lw=LW)
    ax.plot(d.ops, d.postprojected_doublet_error_mEh, label='Post-projected parent prefix', color=PALETTE['green'], lw=LW)
    ax.axhline(1.6, linestyle='--', linewidth=1.4, color=PALETTE['gray'], label='1.6 mEh model-space benchmark')
    ax.set_yscale('log')
    ax.set_xlabel('Final-sequence prefix length')
    ax.set_ylabel('Energy error relative to exact target doublet (mEh)')
    ax.set_xlim(d.ops.min(), d.ops.max())
    stylize(ax)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=2, frameon=False)
    save(fig, 'fig3_projected_prefix_divergence.pdf')


def figS1():
    d = pd.read_csv(HERE / 'source_data_natural_occupations_v5.csv')
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    for case, marker, color in [('18q', 'o', PALETTE['blue']), ('24q', 's', PALETTE['orange'])]:
        x = d[d.case == case]
        ax.plot(x.orbital_index, x.absolute_error, marker=marker, ms=MS, lw=2.0, color=color, label=case.replace('q', ' qubits'))
    ax.set_yscale('log')
    ax.set_xlabel('Natural-orbital index')
    ax.set_ylabel('Absolute occupation-number error')
    ax.set_title('One-particle occupations remain close to exact doublet values')
    stylize(ax)
    ax.legend(frameon=False)
    save(fig, 'figS5_natural_occupation_errors.pdf')


def figS2():
    d = pd.read_csv(HERE / 'source_data_multistart_v5.csv')
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    styles = [('18q', 'o', PALETTE['blue']), ('24q', 's', PALETTE['orange'])]
    for case, marker, color in styles:
        x = d[d.case == case].reset_index(drop=True)
        xpos = np.arange(1, len(x) + 1)
        ok = x.success.astype(str).str.lower().eq('true')
        ax.scatter(xpos[ok], x.loc[ok, 'dev_mEh'], marker=marker, s=38,
                   color=color, label=case.replace('q', ' qubits') + ' (converged)')
        ax.scatter(xpos[~ok], x.loc[~ok, 'dev_mEh'], marker=marker, s=38,
                   facecolors='none', edgecolors=color, linewidths=1.4,
                   label=case.replace('q', ' qubits') + ' (iteration cap)')
    ax.axhline(1.6, linestyle='--', linewidth=1.3, color=PALETTE['gray'],
               label='1.6 mEh model-space benchmark')
    ax.set_xlabel('Perturbed-start reoptimization')
    ax.set_ylabel('Final energy error (mEh)')
    ax.set_title('Compact ansatz solutions are locally robust')
    stylize(ax)
    ax.legend(frameon=False, fontsize=9)
    save(fig, 'figS6_multistart_robustness.pdf')


def figS3():
    d = pd.read_csv(HERE / 'source_data_representation_scaling.csv')
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.loglog(d.theta, d.leak_amp, marker='o', ms=MS, lw=2.0, color=PALETTE['blue'], label=r'Leakage amplitude $\|Qe^{\theta A}|\Phi\rangle\|$')
    ax.loglog(d.theta, d.compressed_mismatch_norm, marker='s', ms=MS, lw=2.0, color=PALETTE['orange'], label='Compressed-evolution mismatch')
    m = d.postprojected_infidelity > 0
    ax.loglog(d.loc[m, 'theta'], d.loc[m, 'postprojected_infidelity'], marker='^', ms=MS, lw=2.0, color=PALETTE['green'], label='Post-projection infidelity')
    ax.set_xlabel(r'Amplitude $\theta$')
    ax.set_ylabel('Norm / infidelity')
    ax.set_title('Single-generator representation mismatch in the Cu active space')
    stylize(ax)
    ax.legend(frameon=False, fontsize=9)
    save(fig, 'fig2_single_generator_scaling.pdf')


def figS4():
    d = pd.read_csv(HERE / 'source_data_projected_negative_control_prefix_v5.csv')
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.plot(d.ops, d.full_fidelity, lw=LW, color=PALETTE['blue'], label='Fidelity: projected vs bare')
    ax.plot(d.ops, d.doublet_weight, lw=LW, color=PALETTE['orange'], label='Bare-state doublet weight')
    ax.plot(d.ops, d.postprojected_doublet_fidelity, lw=LW, color=PALETTE['green'], label='Fidelity after post-projection')
    ax.set_xlabel('Final-sequence prefix length')
    ax.set_ylabel('State overlap / weight')
    ax.set_ylim(0, 1.03)
    ax.set_title('State diagnostics for final-angle prefixes of the Cu negative control')
    stylize(ax)
    ax.legend(frameon=False)
    save(fig, 'figS1_prefix_state_diagnostics.pdf')


def figS5():
    d = pd.read_csv(HERE / 'source_data_hubbard4_representation_v5.csv')
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.plot(d.U_over_t, d.projected_bare_fidelity, marker='o', ms=MS, lw=2.0, color=PALETTE['blue'], label='Projected vs bare fidelity')
    ax.plot(d.U_over_t, d.projected_postprojected_fidelity, marker='s', ms=MS, lw=2.0, color=PALETTE['orange'], label='Projected vs post-projected bare fidelity')
    ax.set_xlabel(r'$U/t$')
    ax.set_ylabel('Fidelity')
    ax.set_title('Independent Hubbard-model stress test of representation equivalence')
    stylize(ax)
    ax.legend(frameon=False)
    save(fig, 'figS2_hubbard_representation_audit.pdf')


def figS6():
    d = pd.read_csv(HERE / 'source_data_local_suzuki_predictor_v5.csv')
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    for case, marker, color in [('18q', 'o', PALETTE['blue']), ('24q', 's', PALETTE['orange']), ('NO', '^', PALETTE['green'])]:
        q = d[(d['case'] == case) & (d['comm_relF'] > 0) & (d['local_infidelity'] > 0)]
        ax.loglog(q['theta_abs'], q['local_infidelity'], linestyle='None', marker=marker, markersize=5.5, color=color, label=case)
    ax.set_xlabel(r'Optimized amplitude magnitude $|\theta|$')
    ax.set_ylabel('Single-application Suzuki-2 infidelity')
    ax.set_title('Among noncommuting factors, amplitude dominates local Suzuki difficulty')
    stylize(ax)
    ax.legend(frameon=False)
    save(fig, 'figS3_local_suzuki_amplitude.pdf')


def figS7():
    d = pd.read_csv(HERE / 'source_data_measurement_shots_v5.csv')
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.plot(d.sigma_mEh, d.shots_18q, marker='o', ms=MS, lw=2.0, color=PALETTE['blue'], label='18 qubits')
    ax.plot(d.sigma_mEh, d.shots_24q, marker='s', ms=MS, lw=2.0, color=PALETTE['orange'], label='24 qubits')
    ax.set_yscale('log')
    ax.invert_xaxis()
    ax.set_xlabel('Target statistical uncertainty (mEh)')
    ax.set_ylabel('Optimally allocated shots for greedy QWC groups')
    ax.set_title('Measurement burden at the mEh scale')
    stylize(ax)
    ax.legend(frameon=False)
    save(fig, 'figS4_measurement_shots.pdf')


if __name__ == '__main__':
    for f in (fig1, fig2, figS1, figS2, figS3, figS4, figS5, figS6, figS7):
        f()
    print('Wrote manuscript/SI figures to', HERE)
