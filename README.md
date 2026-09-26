# Spin-Faithful ADAPT-VQE for an Open-Shell Copper Center

Code, checkpoints, and validation data for the manuscript **From Projected Subspaces to Physical Circuits: Spin-Faithful ADAPT-VQE for an Open-Shell Copper Center** by Lucia Malíčková and Petr Klenovský.

## What this repository tests

For an open-shell VQE calculation, four layers must be distinguished:

1. the target spin sector of the electronic Hamiltonian;
2. the state optimized by the variational algorithm;
3. the physical fermionic unitary sequence;
4. the synthesized qubit circuit and its measurement cost.

The repository contains an explicit counterexample showing that an accurate state obtained after exact projection into a doublet subspace need not be represented by the corresponding bare fermionic sequence.

## Final validated benchmarks

| Model | Qubits | Target-doublet energy (Eh) | Compact physical ansatz |
|---|---:|---:|---:|
| CAS(15e,9o) | 18 | -2518.989067369946 | 50 operator applications |
| CAS(21e,12o) | 24 | -2518.989131438323 | 61 operator applications |

Both active-space Hamiltonians contain a lower quartet. The doublet is the **target state** because the benchmark is intended to represent oxidized Cu(II); the quartet is therefore called non-target, not mathematically unphysical. The calculated ordering is not claimed to be a converged spin-state prediction for a real blue-copper protein.

Main numerical findings:

- projected-doublet ADAPT: 1.572 mEh with 123 operators;
- the same angles with bare generators: fidelity 0.07416 to the projected state, <S^2>=2.381596, error 464.711 mEh;
- restricted physical S+D: 153.629 mEh;
- restricted physical S+D+T: 141.392 mEh;
- restricted physical S+D+T with repeats: 13.166 mEh;
- generalized physical spin-preserving ADAPT, 18q: 1.492674 mEh, fidelity 0.982614, <S^2>=0.75;
- generalized physical spin-preserving ADAPT, 24q: 1.181293 mEh, fidelity 0.990681, <S^2>=0.75.

Circuit synthesis is validated separately. The 24q Suzuki-2 r=1 circuit passes the strict circuit-equivalence gate with 0.0906 mEh synthesis error and 0.999694 fidelity, but still requires 22,574 CX gates at abstract-gate-set depth 26,409. The 18q Suzuki-2 r=4 circuit reaches 0.999484 fidelity and 0.135 mEh synthesis error at 61,482 CX and depth 71,636; its total error is 1.628 mEh.

The commonly used 1.6 mEh value is treated only as a **model-space algorithmic benchmark**. It is not a claim of 1 kcal/mol predictive accuracy for a real protein, laccase, or the STO-3G model.

## Layout

```text
src/                       physical-spin implementation and validators
src/search/                adaptive search, reoptimization, and pruning scripts
data/18q/                  CAS(15e,9o) Hamiltonian and best-50 checkpoint
data/24q/                  CAS(21e,12o) Hamiltonian and best-61 checkpoint
results/final/             validated exact and Qiskit/Aer outputs
results/ablation/          restricted-pool/repeat ablations
results/excluded/          failed validation extensions kept for transparency
paper/                     manuscript, bibliography, source data, figure builder
docs/                      audit and validation notes
legacy/initial_repository_snapshot/
                           complete snapshot of the superseded initial repository
```

The `legacy/` tree is retained for provenance and contains the previous projected-ADAPT workflow and exploratory PES material. It should not be used to reproduce the new manuscript.

## Reproduction

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

The colleague-tested stack is pinned in `requirements.txt`. No physical-QPU final-energy result is claimed; stored circuit results are Qiskit/Aer statevector validation and abstract-basis transpilation results.

The final 18q FCIDUMP SHA256 is:

```text
4fe2f74860c4a7b1af7e677f62f4d59a1571c2758ff427ba8627607dcd8c8536
```

The preliminary Cu-S geometry scan failed its own validation gates and is explicitly excluded from geometry-dependent scientific claims.
