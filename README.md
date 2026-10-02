# Representation Inequivalence in Open-Shell ADAPT-VQE

This repository contains the code, fixed Hamiltonians, checkpoints, validation records, and figure source data for

**Representation Inequivalence in Open-Shell ADAPT-VQE: Projected Optimization versus Full-Space Implementation**

by Lucia Malíčková and Petr Klenovský.

The current paper is being prepared as a Research Paper for *Electronic Structure*. The repository is the data/reproducibility archive referenced by the manuscript. The earlier JCP submission files are retained under `paper/jcp_submission/` for provenance only and are not the current manuscript.

## Current scientific result

A projected/reduced-space ADAPT-VQE optimization and the corresponding unprojected parent-generator implementation are not automatically the same variational construction. For a parent anti-Hermitian generator `A` and target-space projector `P`, generator-by-generator equivalence for arbitrary amplitudes requires target-space invariance,

```
Q A P = 0,    Q = I - P.
```

The invariant-subspace algebra is standard; the contribution of the work is the explicit representation-consistency framework, finite-sequence/state-specific diagnostics, and controlled numerical demonstration of the resulting mismatch in open-shell ADAPT-VQE workflows.

The exact projector-based diagnostics are small-system validation tools, not a universal scalable certification protocol. For the spin-resolved workflows studied here, preservation of `M_S` together with `[A,S^2]=0` provides a projector-free sufficient generator-level test of target-spin invariance.

## Key validated results

- Canonical Cu CAS(15e,9o) projected-doublet ADAPT: **1.5719 mEh** with 121 operators.
- Same parent sequence and final optimized amplitudes in full space: **464.730 mEh**, projected-state fidelity **0.07343**, doublet weight **0.45549**.
- Best tested target-spin-constrained full-space reoptimization of the same fixed sequence: **8.4815 mEh** at `p_D=0.999`; no global optimum is claimed.
- Independent NO benchmark: projected **1.4286 mEh**, same-angle parent **24.537 mEh**, fidelity **0.97205**.
- OH near-equivalent counterexample: projected **1.0551 mEh**, same-angle parent **1.9328 mEh**, fidelity **0.999229**.
- Four-site Hubbard stress tests confirm that the representation issue is not specific to the Cu molecular Hamiltonian.
- Representation-faithful positive controls: **1.4927 mEh** (18 qubits) and **1.1813 mEh** (24 qubits), with reduced/full-space agreement to numerical precision.
- Fixed-seed selective 18q Qiskit validation: synthesis error **0.0950208 mEh**, circuit-to-fermionic fidelity **0.99958246**, **19,660 CX** gates, depth **23,016**; all four predeclared circuit-equivalence criteria pass.

The Cu prefix plots are **final-angle prefix reconstructions**: prefix `k` uses the first `k` factors and first `k` entries of the final globally optimized 121-parameter vector. They are diagnostics of the final ansatz, not the per-iteration ADAPT convergence history, and therefore need not be monotonic.

## Scientific scope

The work does **not** claim:
- a new spin-adapted ADAPT algorithm;
- a new invariant-subspace theorem;
- a chemically converged blue-copper prediction;
- quantum advantage;
- hardware-level resource bounds;
- a globally optimal circuit or fixed-sequence parameter solution.

The 24-qubit Hamiltonian is retained as a fixed algorithmic FCIDUMP benchmark because the dedicated 12-orbital selection driver was not preserved; no geometry-to-Hamiltonian reproducibility or active-space-convergence claim is made for that case.

## Repository layout

```text
src/                         implementation and validation scripts
data/18q/                    canonical CAS(15e,9o) Hamiltonian and checkpoints
data/24q/                    fixed CAS(21e,12o) algorithmic benchmark
benchmarks/                  Cu protocol and independent NO/OH benchmarks
results/                     validated outputs and audit records
paper/                       figure source data and deterministic figure generator
paper/jcp_submission/        archived previous JCP submission materials
docs/                        validation and data-audit notes
legacy/                      superseded workflow retained for provenance
```

See `paper/CURRENT_MANUSCRIPT.md` for the current manuscript-level summary and `REPRODUCIBILITY.md` for rerun commands.

Canonical 18q FCIDUMP SHA256:

```text
4fe2f74860c4a7b1af7e677f62f4d59a1571c2758ff427ba8627607dcd8c8536
```
