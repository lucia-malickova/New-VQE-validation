# Current manuscript status

## Submission state

**Journal target:** *Electronic Structure* (IOP Publishing)  
**Article type:** Research Paper  
**Submission-ready state:** 3 October 2026

## Title

**Representation Inequivalence in Open-Shell ADAPT-VQE: Projected Optimization versus Full-Space Implementation**

Authors: Lucia Malíčková and Petr Klenovský.

The repository is the reproducibility/data archive cited by the current manuscript. Earlier JCP-formatted materials are preserved under jcp_submission/ for provenance and are not the current article.

## Current scientific framing

The manuscript studies a representation-consistency question: when an ADAPT-VQE calculation selects and optimizes projected/reduced-space operators, under what condition may the resulting parent labels and amplitudes be identified with an unprojected full-space implementation?

For a determinant excitation T, the corresponding anti-Hermitian parent generator is A=T-T†. For target-space projector P and Q=I-P, arbitrary-amplitude generator-level equivalence requires QAP=0. The invariant-subspace algebra itself is standard. The methodological contribution is to use it as an explicit cross-representation consistency check, formulate finite-sequence/state-specific diagnostics, and quantify the mismatch on controlled open-shell benchmarks.

The target-doublet weight used in constrained full-space reoptimization is p_D=||U†Psi||^2.

## Headline numerical results

| Test | Projected/reduced result | Parent/full-space result |
|---|---:|---:|
| Cu CAS(15e,9o) | 1.5719 mEh | 464.730 mEh, F=0.07343, pD=0.45549 |
| NO CAS(7e,6o) | 1.4286 mEh | 24.537 mEh, F=0.97205 |
| OH CAS(7e,5o) | 1.0551 mEh | 1.9328 mEh, F=0.999229 |

For Cu, target-spin-constrained reoptimization of the same fixed 121-generator parent sequence reaches 8.4815 mEh in the best tested basin at pD=0.999; no global-optimum or impossibility result is claimed.

The OH result is an explicit near-equivalent counterexample: nonzero invariance defect does not by itself imply a large practical error. Four-site Hubbard tests demonstrate that the mechanism is not molecule specific.

Representation-faithful spin-preserving positive controls recover reduced/full-space agreement to numerical precision. The selective 18q Qiskit validation passes all four predeclared circuit-equivalence gates at 0.0950208 mEh synthesis error and 0.99958246 circuit-to-fermionic fidelity.

## Figure and prefix provenance

The benchmark-system overview uses the archived Cu geometry as an x-z projection, with y perpendicular to the panel; NO/OH distances are those stored in the independent benchmark metadata. The Hubbard panel is schematic and labels the uniform nearest-neighbor hopping t and uniform on-site interaction U.

The Cu prefix plots are diagnostic reconstructions of the final ansatz. Prefix k uses the first k factors and first k amplitudes of the final globally optimized 121-parameter vector. They are not the per-iteration ADAPT convergence history and need not be monotonic.

## Scope limitations

- The Cu active-space Hamiltonian is an algorithmic benchmark, not a chemically converged blue-copper prediction.
- The 24q Hamiltonian has incomplete geometry-to-active-space generation provenance.
- Exact projector-based diagnostics are benchmark-scale classical oracles, not universal scalable quantum primitives.
- The generalized spin-preserving pool is used as a positive control, not claimed as a new preferred ADAPT algorithm.
- Circuit counts are abstract-basis compilation results without a device coupling map.
