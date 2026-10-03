# Representation Inequivalence in Open-Shell ADAPT-VQE

This repository contains the code, fixed Hamiltonians, checkpoints, validation records, and figure source data for

**Representation Inequivalence in Open-Shell ADAPT-VQE: Projected Optimization versus Full-Space Implementation**

by Lucia Malíčková and Petr Klenovský.

**Submission state (3 October 2026):** the manuscript is prepared as a Research Paper for *Electronic Structure*. The repository is the public reproducibility archive cited by the submission. The earlier JCP submission files remain under paper/jcp_submission/ for provenance only and are not the current manuscript.

## Scientific result

A projected/reduced-space ADAPT-VQE optimization and the corresponding unprojected parent-generator implementation are not automatically the same variational construction. For a determinant excitation T, the anti-Hermitian parent generator is A = T - T†. For a target-space projector P and Q = I - P, generator-by-generator equivalence at arbitrary amplitudes requires QAP = 0.

The invariant-subspace algebra is standard. The contribution of the work is the explicit representation-consistency framework, finite-sequence/state-specific diagnostics, and controlled numerical demonstration of the resulting mismatch in open-shell ADAPT-VQE workflows.

The exact projector-based diagnostics are small-system validation tools, not a universal scalable certification protocol. For the spin-resolved workflows studied here, preservation of M_S together with [A,S^2]=0 provides a projector-free sufficient generator-level test of target-spin invariance.

## Key validated results

- Canonical Cu CAS(15e,9o) projected-doublet ADAPT: **1.5719 mEh** with 121 operators.
- Same parent sequence and final optimized amplitudes in full space: **464.730 mEh**, projected-state fidelity **0.07343**, target-doublet weight **0.45549**.
- Best tested target-spin-constrained full-space reoptimization of the same fixed sequence: **8.4815 mEh** at p_D=0.999, where p_D=||U†Psi||^2; no global optimum is claimed.
- Independent NO benchmark: projected **1.4286 mEh**, same-angle parent **24.537 mEh**, fidelity **0.97205**.
- OH near-equivalent counterexample: projected **1.0551 mEh**, same-angle parent **1.9328 mEh**, fidelity **0.999229**.
- Four-site Hubbard stress tests show that the mechanism is not molecule specific.
- Representation-faithful positive controls: **1.4927 mEh** (18 qubits) and **1.1813 mEh** (24 qubits), with reduced/full-space agreement to numerical precision.
- Fixed-seed selective 18q Qiskit validation: synthesis error **0.0950208 mEh**, circuit-to-fermionic fidelity **0.99958246**, **19,660 CX** gates, depth **23,016**; all four predeclared circuit-equivalence criteria pass.

The Cu prefix plots are **final-angle prefix reconstructions**: prefix k uses the first k factors and first k entries of the final globally optimized 121-parameter vector. They are diagnostics of the final ansatz, not the per-iteration ADAPT convergence history, and therefore need not be monotonic.

## Review entry points

- paper/CURRENT_MANUSCRIPT.md: concise manuscript-level summary and scope.
- paper/ELECTRONIC_STRUCTURE_SUBMISSION.md: submission-state abstract, article type, and manuscript/repository mapping.
- REPRODUCIBILITY.md: rerun commands and canonical input hashes.
- paper/make_all_figures.py: regenerates all current manuscript/SI figures, including the benchmark-system overview.
- results/revision_v9/ and results/revision_v10/: fixed-parent energy/fidelity reoptimization records.
- results/revision_v12/selective_suzuki_qiskit/: fixed-seed selective-Qiskit validation.

## Scientific scope

The work does **not** claim a new invariant-subspace theorem, a new preferred spin-adapted ADAPT algorithm, a chemically converged blue-copper prediction, quantum advantage, hardware-level resource bounds, or a globally optimal circuit/fixed-sequence parameter solution.

The 24-qubit Hamiltonian is retained as a fixed algorithmic FCIDUMP benchmark because the dedicated 12-orbital selection driver was not preserved; no geometry-to-Hamiltonian reproducibility or active-space-convergence claim is made for that case.

## Repository layout

    src/                         implementation and validation scripts
    data/18q/                    canonical CAS(15e,9o) Hamiltonian and checkpoints
    data/24q/                    fixed CAS(21e,12o) algorithmic benchmark
    benchmarks/                  Cu protocol and independent NO/OH benchmarks
    results/                     validated outputs and audit records
    paper/                       figure source data and deterministic figure generators
    paper/jcp_submission/        archived previous JCP submission materials
    docs/                        validation and data-audit notes
    legacy/                      superseded workflow retained for provenance

Canonical 18q FCIDUMP SHA256:

    4fe2f74860c4a7b1af7e677f62f4d59a1571c2758ff427ba8627607dcd8c8536
