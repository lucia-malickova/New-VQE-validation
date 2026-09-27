# Submission revision v6

This revision addresses the remaining detailed-referee concerns after v5.

- Adds the exact canonical Cu first-shell geometry and the archived ROHF/AVAS/CASSCF generation protocol.
- States explicitly that the dense exact projector is a small-system classical validation oracle, not a scalable quantum primitive.
- Corrects lifted reduced-space versus full-space notation.
- Distinguishes method comparison from a one-variable-at-a-time ablation.
- Quantifies generalized-pool null-space-basis sensitivity with three independent O(2)-rotated pool bases.
- Reassesses the Suzuki commutator heuristic against an amplitude-only baseline and narrows the claim: noncommutativity is a structural separator; amplitude dominates the ranking among noncommuting factors.
- Clarifies the phase-minimized state-distance diagnostic used alongside the D_seq vector-norm certificate.
- Records that seed_transpiler was not pinned in the stored Qiskit resource campaign.
- Documents the missing dedicated 24q active-space-selection driver instead of reconstructing it from inference.
