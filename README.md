# Representation-Equivalence Audits for Open-Shell ADAPT-VQE

This repository contains the code, fixed Hamiltonians, validation data, and final pre-submission sources for

**From Projected Subspaces to Full-Space Implementations: Representation-Equivalence Audits for Open-Shell ADAPT-VQE**

by Lucia Malíčková and Petr Klenovský.

## Current submission state

Prepared for initial submission to **The Journal of Chemical Physics (JCP)** as a regular Article.

Canonical sources:
- `paper/manuscript.tex` — AIP/JCP REVTeX main manuscript;
- `paper/supporting_information.tex` — separate JCP supplementary material;
- `paper/references.bib`;
- `paper/jcp_submission/` — mirrored self-contained submission-support source set.

The current pre-submission state is **v20**. Main/JCP mirrors and SI/JCP mirrors are identical. CAS, ROHF, AVAS, CASSCF, and SLSQP are now defined explicitly at first use in the main manuscript; the standalone Supplement independently defines ROHF, AVAS, CAS, and SLSQP. Every figure/table is explicitly called out and labeled. The final technical wording makes clear that the projector-free total-spin sufficiency statement is used for the $M_S$-preserving generators considered in this work.

## Key validated results

- Canonical Cu projected-doublet ADAPT: 1.5719 mEh with 121 operators.
- Same parent sequence and amplitudes in full space: 464.730 mEh, projected-state fidelity 0.07343.
- Best tested full-space reoptimization of the same fixed sequence under $p_D\ge0.999$: 8.4815 mEh.
- Representation-faithful positive-control ansätze: 1.4927 mEh (18q) and 1.1813 mEh (24q).
- Fixed-seed selective 18q Qiskit compilation: 0.0950208 mEh synthesis error, fidelity 0.99958246, 19,660 CX, depth 23,016; all four predeclared circuit-equivalence gates PASS.

## Scientific scope

The work does **not** claim a new spin-adapted ADAPT algorithm. It audits the narrower identification step between a reduced/projected variational construction and a proposed parent full-space implementation/resource estimate.

The Cu model is a chemically motivated, exactly diagonalizable algorithmic benchmark. No chemically converged blue-copper prediction, quantum-advantage claim, physical-QPU final-energy result, globally optimal circuit, global fixed-sequence representability theorem, or basis-independent minimal ansatz is claimed.

## Repository layout

```text
src/                       implementation and validation scripts
data/18q/                  canonical CAS(15e,9o) Hamiltonian and 18q checkpoints
data/24q/                  fixed CAS(21e,12o) algorithmic benchmark
benchmarks/                geometry/protocol and independent NO/OH benchmarks
results/                   validated outputs and revision-specific audits
paper/                     canonical JCP sources, figure source data, figure generator
paper/jcp_submission/      mirrored JCP submission-support materials
docs/                      final submission/validation/data-audit notes
legacy/                    superseded workflow retained only for provenance
```

See `REPRODUCIBILITY.md` for exact rerun commands.

Canonical 18q FCIDUMP SHA256:

```text
4fe2f74860c4a7b1af7e677f62f4d59a1571c2758ff427ba8627607dcd8c8536
```

The 24q Hamiltonian is retained as a fixed FCIDUMP benchmark because the dedicated 12-orbital selection driver was not preserved.


Revision v20 condenses the main-text scope section so that it states the positive scope of the claims rather than repeating detailed exclusions already documented elsewhere.
