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

Every figure and table in the main manuscript and supplementary material has an explicit textual callout and LaTeX label.
Conflict of Interest and CRediT statements are inserted. Both authors must approve the final manuscript and individual CRediT roles before actual submission.

A graphical abstract is not required for initial JCP submission; AIP recommends a highlight image at revision.

## Referee-hardening revision v12

This revision sharpens the JCP submission and closes the previously open 18q circuit-validation item. It (i) distinguishes the representation-equivalence audit from spin projection, symmetry-preserving/spin-adapted pool design, contextual-subspace VQE, and exact spin-adapted factorization; (ii) adds a four-layer validation schematic to the main manuscript; (iii) states explicitly that the generalized spin-preserving pool is a positive control rather than a claimed new scalable algorithm; (iv) adds targeted symmetry-projection and symmetry-preserving-circuit references; and (v) records a fixed-seed Qiskit re-transpilation of the selectively refined 18q Suzuki-2 schedule.

The confirmed 18q selective schedule uses $r_9=4$, $r_{21}=r_{25}=r_{34}=2$, and $r=1$ elsewhere. With Qiskit 2.2.3, optimization level 1, basis `{rz,sx,x,cx}`, and `seed_transpiler=9272026`, it passes all four predeclared circuit-equivalence gates: synthesis error 0.0950208 mEh, circuit-to-fermionic fidelity 0.99958246, $S^2=0.75072418$, and target-$M_S$ weight 0.9999999999994. The abstract compiled resources are 19,660 CX gates and depth 23,016.

During this rerun, two repository provenance/interface defects were exposed and repaired: the committed 18q FCIDUMP was not the canonical validated file, and the standalone helper incorrectly assumed that the intentionally compact best-50 checkpoint already contained serialized generators. Revision v12 restores the canonical FCIDUMP, adds the archived full serialized best-50 checkpoint, lets the selective validator accept either checkpoint form, and restores Python-3.9 compatibility.

## Scientific scope

The work does **not** claim a new spin-adapted ADAPT algorithm.
It audits the narrower failure mode that arises when a projected/reduced-space variational optimization is identified with an unprojected parent full-space implementation without establishing representation equivalence.

The canonical Cu same-angle mismatch is severe (1.5719 mEh projected versus 464.730 mEh full-space).
Full-space energy and direct-fidelity reoptimization recover much of the target state, so the result is interpreted primarily as a representation/parameter-transfer failure rather than proof that the parent sequence is intrinsically useless.

Independent NO, OH, and Hubbard tests, representation-faithful positive controls, circuit validation, and final-energy measurement audits define the scope.

## Repository layout

```text
src/                       implementation and validation scripts
data/18q/                  canonical CAS(15e,9o) Hamiltonian and compact ansatz
data/24q/                  fixed CAS(21e,12o) algorithmic benchmark
results/                   validated numerical outputs and revision audits
paper/                     canonical JCP manuscript sources and figure data
paper/jcp_submission/      mirrored JCP submission-support materials
docs/                      reproducibility and revision notes
legacy/                    superseded workflow retained only for provenance
```

See `REPRODUCIBILITY.md` for commands and environment.

Canonical 18q FCIDUMP SHA256:
```text
4fe2f74860c4a7b1af7e677f62f4d59a1571c2758ff427ba8627607dcd8c8536
```

The 24q Hamiltonian is a fixed FCIDUMP benchmark because the dedicated 12-orbital selection driver was not preserved.

No physical-QPU final-energy result, quantum-advantage claim, chemically converged blue-copper prediction, global fixed-sequence representability theorem, or basis-independent minimal-ansatz claim is made.
