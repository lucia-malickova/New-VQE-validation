# Submission v17 — final pre-submission technical state

This note records the final source-level technical tightening before JCP submission.

- The main manuscript explicitly restricts the projector-free $[A,\hat S^2]=0$ sufficiency statement to the $M_S$-preserving generators considered in this work.
- The Supplement states the corresponding combined condition: the production pool uses $M_S$-preserving excitation patterns together with $[B,\hat S^2]\simeq0$, while exact-projector leakage remains a classical validation oracle.
- The manuscript distinguishes amplitude-independent generator-level equivalence from accidental agreement at an isolated finite parameter vector.
- The practical significance paragraph explicitly avoids claiming that prior projected-space studies necessarily make the representation-identification error.
- The 0.1 mEh synthesis threshold is justified as 6.25% of the 1.6 mEh model-space benchmark and was fixed before selective schedule search.
- Full Lucia Malíčková affiliation and Euro-Q-Exa acknowledgement are restored.
- Fixed-seed selective 18q Qiskit validation remains PASS; numerical scientific results are otherwise unchanged.
- Final source audit is enforced by `paper/audit_submission.py`.
- `.github/workflows/package_jcp.yml` regenerates figures, runs the audit, and creates the self-contained submission bundle.
