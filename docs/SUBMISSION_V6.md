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


## Additional v7 clarification

- The main text now distinguishes the fully archived geometry-to-FCIDUMP protocol for the canonical 18q Cu Hamiltonian from the fixed 24q FCIDUMP, whose dedicated 12-orbital selection driver was not preserved.
- The relation-to-prior-work section now states explicitly that projected-pool methods are valid reduced-space methods; the audited failure concerns assigning an unprojected parent-generator implementation/resource count without proving representation equivalence.
- The SI writes the equivalence with the explicit lift (U\exp[\theta(U^\dagger A U)]U^\dagger), removing any dimensional ambiguity between reduced- and full-space operators.
