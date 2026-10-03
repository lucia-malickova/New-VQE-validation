# Electronic Structure submission state

**Article type:** Research Paper  
**Target journal:** *Electronic Structure* (IOP Publishing)  
**Submission-ready date:** 3 October 2026

## Title

**Representation Inequivalence in Open-Shell ADAPT-VQE: Projected Optimization versus Full-Space Implementation**

## Abstract

Adaptive variational quantum eigensolver workflows may optimize projected generators and then implement the corresponding unprojected parent generators with the same amplitudes, although the two constructions need not be equivalent. We formulate target-subspace invariance, QAP = 0, as a representation-consistency criterion and derive finite-sequence diagnostics. In a CAS(15e,9o) Cu(II) benchmark, projected-doublet ADAPT-VQE reaches 1.5719 mEh, whereas the same amplitudes in the parent sequence give 464.730 mEh and fidelity 0.07343. Reoptimizing the same fixed parent sequence under a target-doublet-weight constraint p_D >= 0.999 lowers the best tested error to 8.482 mEh, identifying parameter transfer as the dominant source of the severe failure without proving global sequence insufficiency. Independent NO and OH benchmarks show that the representation mismatch can range from substantial to nearly negligible, while Hubbard-chain tests demonstrate that the mechanism is not specific to molecular Hamiltonians. Spin-preserving full-space controls recover representation equivalence to numerical precision. Explicit Suzuki/Qiskit tests further show that circuit synthesis forms a separate validation layer. These results provide a practical framework for distinguishing reduced-space variational accuracy from the accuracy of a proposed full-space implementation.

## Reproducibility mapping

- Canonical Cu 18q Hamiltonian: data/18q/active.FCIDUMP
- Cu 18q Hamiltonian/geometry protocol: benchmarks/cu18/
- Cu projected negative-control checkpoint/audit: results/revision_v6/cu/
- Fixed-parent energy reoptimization: results/revision_v9/cu/
- Fixed-parent fidelity reoptimization: results/revision_v10/cu/
- NO and OH independent benchmark files: benchmarks/no/, benchmarks/oh/
- Selective-Qiskit validation: results/revision_v12/selective_suzuki_qiskit/
- Figure source data: paper/source_data_*.csv
- Current figure regeneration: python paper/make_all_figures.py

The prior JCP submission directory is archival provenance only.
