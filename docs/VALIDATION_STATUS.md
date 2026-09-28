# Final validation status before JCP submission

## Passed

- canonical 18q FCIDUMP provenance and hash audit;
- deterministic 18q best-50 reconstruction and exact physical validation;
- deterministic 24q best-61 reconstruction and exact physical validation;
- projected-vs-parent representation-equivalence negative control;
- fixed-parent full-space energy and fidelity reoptimization controls;
- independent NO, OH, and Hubbard representation tests;
- full generalized spin-preserving pool audit;
- full/reduced positive-control prefix equivalence for 18q and 24q;
- normalization-controlled null-space basis-sensitivity audit;
- 18q/24q Qiskit prefix smoke tests;
- 24q full Suzuki-2 r=1 strict circuit-equivalence validation;
- fixed-seed selective 18q Qiskit validation: PASS;
- multi-start, 1-RDM/natural-occupation, and final-energy measurement-resource audits;
- main/SI labels, references, bibliography, mirror synchronization, and local figure-file audit.

## Selective 18q Qiskit PASS

Schedule: $r_9=4$, $r_{21}=r_{25}=r_{34}=2$, all others $r=1$.

- Qiskit 2.2.3;
- optimization level 1;
- abstract basis `{rz,sx,x,cx}`;
- `seed_transpiler=9272026`;
- synthesis error 0.0950208 mEh;
- circuit-to-fermionic fidelity 0.99958246;
- $\langle S^2\rangle=0.75072418$;
- target-$M_S$ weight 0.9999999999994;
- 19,660 CX, depth 23,016.

The resource numbers are abstract-basis compiler outputs without a device coupling map and are not hardware-resource bounds.

## Qualified comparison

Uniform 18q Suzuki-2 r=4 gives fidelity 0.999484, synthesis error 0.135 mEh, 61,482 CX, and depth 71,636. It is a close circuit-to-fermionic reproduction but fails the predeclared 0.1 mEh synthesis-energy gate.

## Excluded

The exploratory Cu-S geometry scan failed its pre-specified CASSCF-convergence and active-orbital-continuity gates and is excluded from scientific claims.

## Not claimed

No physical-QPU final-energy result, chemically converged Cu prediction, quantum advantage, global optimizer proof, globally optimal circuit, or basis-independent minimal ansatz is claimed.
