# Submission revision v12 — referee hardening

This revision is deliberately conservative: it strengthens positioning and validation language without inventing new numerical results.

## Changes already incorporated

1. **Novelty positioning.** The Introduction now separates the present representation-equivalence audit from:
   - symmetry projection / symmetry-adapted VQE;
   - symmetry-preserving state-preparation circuits and symmetry-aware or symmetry-complete pools;
   - contextual-subspace / projected reduced-space VQE;
   - exact or closed-form spin-adapted unitary factorization.

2. **Operational failure path.** The Introduction now states explicitly how the ambiguity can enter when projected/compressed operators are optimized in a reduced space and parent labels/amplitudes are later reused to construct a full-space fermionic or qubit implementation.

3. **Positive-control scope.** The generalized spin-preserving pool is explicitly identified as a representation-faithful positive control, not as a claimed new or preferred scalable spin-adapted ADAPT algorithm. The same clarification appears in the SI.

4. **Literature coverage.** Added targeted references for symmetry-preserving state preparation and symmetry-projected VQE:
   - Gard et al., npj Quantum Information 6, 10 (2020);
   - Seki, Shirakawa, and Yunoki, Phys. Rev. A 101, 052340 (2020);
   - Seki and Yunoki, Phys. Rev. A 105, 032419 (2022).

5. **Conceptual figure.** Added a four-layer validation-ladder schematic:
   reduced-space optimization -> parent full-space sequence -> synthesized qubit circuit -> measurement protocol.
   Its caption states the exact first-arrow condition (QA_kP=0).

6. **Cover letter.** Strengthened the novelty distinction without claiming that prior projected or spin-adapted methods are incorrect.

## Intentionally pending numerical item

The SI currently reports that a nonuniform 18q Suzuki-2 schedule passes the strict state-level criteria in a determinant-equivalent emulator, but it does **not** report a re-transpiled Qiskit CX/depth result for that schedule.

A separate rerun package fixes:
- schedule: r_9=4, r_21=r_25=r_34=2, all other r=1;
- Qiskit 2.2.3 / Aer 0.17.2 / Nature 0.7.2;
- basis {rz,sx,x,cx};
- optimization level 1;
- seed_transpiler=9272026.

The manuscript must not be changed to claim a passing compiled 18q circuit until the returned validation JSON contains `"PASS": true`.

## Static audit

At creation of this revision:
- all manuscript/SI citation keys resolve in `paper/references.bib`;
- all manuscript `\ref` targets resolve;
- no duplicate manuscript labels were detected.
