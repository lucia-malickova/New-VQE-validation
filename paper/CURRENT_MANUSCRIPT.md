# Current manuscript status

## Title

**Representation Inequivalence in Open-Shell ADAPT-VQE: Projected Optimization versus Full-Space Implementation**

Authors: Lucia Malíčková and Petr Klenovský.

The current journal-formatted submission is being prepared as a Research Paper for *Electronic Structure*. This repository is the reproducibility/data archive cited by the manuscript. Earlier JCP-formatted sources are preserved under `jcp_submission/` as provenance and should not be treated as the current article text.

## Current scientific framing

The manuscript studies a representation-consistency question: when an ADAPT-VQE calculation selects and optimizes projected/reduced-space operators, under what condition may the resulting parent labels and amplitudes be identified with an unprojected full-space implementation?

For an anti-Hermitian parent generator `A`, target-space projector `P`, and `Q=I-P`, the arbitrary-amplitude generator-level equivalence condition is `QAP=0`. The invariant-subspace algebra itself is standard. The methodological contribution is to use it as an explicit cross-representation consistency check, formulate finite-sequence/state-specific diagnostics, and quantify the mismatch on controlled open-shell benchmarks.

## Headline numerical results

| Test | Projected/reduced result | Parent/full-space result |
|---|---:|---:|
| Cu CAS(15e,9o) | 1.5719 mEh | 464.730 mEh, F=0.07343, pD=0.45549 |
| NO CAS(7e,6o) | 1.4286 mEh | 24.537 mEh, F=0.97205 |
| OH CAS(7e,5o) | 1.0551 mEh | 1.9328 mEh, F=0.999229 |

For Cu, target-spin-constrained reoptimization of the same fixed 121-generator parent sequence reaches 8.4815 mEh in the best tested basin at `pD=0.999`; no global-optimum or impossibility result is claimed.

The OH result is an explicit near-equivalent counterexample: nonzero invariance defect does not by itself imply a large practical error. Four-site Hubbard tests demonstrate that the mechanism is not Cu-specific.

Representation-faithful spin-preserving positive controls recover reduced/full-space agreement to numerical precision. The selective 18q Qiskit validation passes all four predeclared circuit-equivalence gates at 0.0950208 mEh synthesis error and 0.99958246 circuit-to-fermionic fidelity.

## Prefix-plot provenance

The Cu prefix plots are diagnostic reconstructions of the final ansatz. Prefix `k` uses:
1. the first `k` factors of the final 121-factor sequence, and
2. the first `k` amplitudes of the final globally optimized 121-parameter vector.

They are **not** the per-iteration ADAPT convergence history and need not be monotonic.

## Scope limitations

- The Cu active-space Hamiltonian is an algorithmic benchmark, not a chemically converged blue-copper prediction.
- The 24q Hamiltonian has incomplete geometry-to-active-space generation provenance.
- Exact projector-based diagnostics are benchmark-scale classical oracles, not universal scalable quantum primitives.
- The generalized spin-preserving pool is used as a positive control, not claimed as a new preferred ADAPT algorithm.
- Circuit counts are abstract-basis compilation results without a device coupling map.
