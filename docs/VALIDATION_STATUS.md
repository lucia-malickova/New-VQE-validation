# Validation status

## Passed

- 18q exact physical best-50 reconstruction
- 24q exact physical best-61 reconstruction
- generalized physical spin-pool audit
- projected-vs-bare equivalence falsification test
- 18q and 24q five-operator Qiskit smoke tests
- 24q full Suzuki-2 r=1 strict circuit-equivalence validation
- multi-start energetic robustness checks
- one-particle RDM / natural-occupation audits
- greedy-QWC measurement-resource audits

## Qualified result

18q full Suzuki-2 r=4:
- fidelity 0.999484
- synthesis-energy error 0.135 mEh
- total target-doublet error 1.628 mEh
- 61,482 CX
- depth 71,636

This is a close circuit-to-fermionic reproduction but fails the deliberately strict 0.1 mEh synthesis-energy gate and lies slightly above the 1.6 mEh model-space benchmark.

## Excluded

The preliminary Cu-S geometry scan failed its own validation gate because the doublet/quartet CASSCF calculations did not converge across the five points and active-orbital continuity failed the preset threshold.

## Not claimed

No physical-QPU final-energy result is reported.
