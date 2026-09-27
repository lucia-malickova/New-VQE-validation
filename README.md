# Representation-Equivalence Audits for Open-Shell ADAPT-VQE

This repository contains the code and validation data for the draft manuscript

**From Projected Subspaces to Physical Circuits: Representation-Equivalence Audits for Open-Shell ADAPT-VQE**

by Lucia Malíčková and Petr Klenovský.

## Scientific scope

The repository is deliberately organized around an end-to-end validation question:

> Does the state optimized in the emulator correspond to the physical fermionic
> sequence, the synthesized qubit circuit, and the measurement claim that is
> ultimately reported?

The answer is not automatic for open-shell VQE.

The benchmark uses an idealized type-1 (blue-copper) Cu(II) center and two
classically exact active-space Hamiltonians:

| Model | Qubits | Target doublet energy (Eh) | Compact physical ansatz |
|---|---:|---:|---:|
| CAS(15e,9o) | 18 | -2518.989067369946 | 50 operator applications |
| CAS(21e,12o) | 24 | -2518.989131438323 | 61 operator applications |

Both active-space Hamiltonians also contain a lower quartet. The doublet is the
**target state** because the model is intended to represent oxidized Cu(II);
the quartet is therefore called *non-target*, not mathematically unphysical.
The calculated spin ordering is **not** claimed to be a converged prediction
for the real protein.

## Main validated findings

1. **Projected-subspace accuracy is not the same as physical-circuit accuracy.**
   A freshly regenerated canonical exact-doublet projected ADAPT sequence reaches
   1.5719 mEh in 121 operators, but the same amplitudes applied to the bare
   fermionic sequence give only 0.07343 fidelity with the projected state,
   <S^2>=2.383537, and 464.730 mEh error. The submission-v5 manuscript also
   reports an independent NO benchmark and model stress tests; the associated
   revision-specific peer-review archive accompanies the manuscript.

2. **Physical spin preservation alone is not sufficient.**
   Restricted occupied-to-virtual spin-adapted pools remain far from the target:
   153.629 mEh (S+D) and 141.392 mEh (S+D+T). Allowing repeats improves the
   restricted S+D+T sequence to 13.166 mEh, but still does not reach the
   1.6 mEh model-space benchmark.

3. **Generalized physical spin-preserving ADAPT succeeds in the emulator.**
   The final fermionic ansätze reach:
   - 18 qubits: 1.492674 mEh, fidelity 0.982614, <S^2>=0.75;
   - 24 qubits: 1.181293 mEh, fidelity 0.990681, <S^2>=0.75.

4. **Circuit synthesis adds another approximation.**
   - 18q, Suzuki-2 r=4: synthesis error 0.135 mEh, fidelity 0.999484,
     61,482 CX, depth 71,636; total doublet error 1.628 mEh.
     This narrowly misses the 1.6 mEh model-space benchmark and also fails the
     stricter 0.1 mEh synthesis-error gate.
   - 24q, Suzuki-2 r=1: synthesis error 0.0906 mEh, fidelity 0.999694,
     22,574 CX, depth 26,409; total doublet error 1.272 mEh.

5. **Measurement cost is severe even before device noise.**
   For the chosen greedy QWC grouping and variance-optimal allocation,
   a 1.6 mEh statistical uncertainty requires approximately
   1.43e8 shots (18q) and 5.46e8 shots (24q).

These numbers are validation results for the specified active-space
Hamiltonians. The commonly used 1.6 mEh value is treated here only as a
**model-space algorithmic benchmark**, not as a claim of 1 kcal/mol accuracy
for the real blue-copper protein or laccase.

## Repository layout

```text
src/                       final physical-spin implementation and validators
data/18q/                  CAS(15e,9o) Hamiltonian and best-50 checkpoint
data/24q/                  CAS(21e,12o) Hamiltonian and best-61 checkpoint
results/final/              validated exact and Qiskit/Aer outputs
results/ablation/           restricted-pool/repeat ablation data
results/excluded/           failed validation extensions kept for transparency
paper/                      manuscript, bibliography, source data, figure builder
docs/                       audit and reproducibility notes
legacy/initial_repository_snapshot/
                           complete snapshot of the superseded initial repository
```

The `legacy/` tree is retained only for provenance. It contains the previous
projected-ADAPT workflow and exploratory PES material and should not be used to
reproduce the new manuscript.

## Quick start

See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for exact commands.

The final colleague-tested environment was:

```text
numpy==2.0.2
scipy==1.13.1
qiskit==2.2.3
qiskit-aer==0.17.2
qiskit-nature==0.7.2
```

No physical-QPU energy result is claimed in this repository. Qiskit/Aer
statevector circuit validation is separated explicitly from device-specific
hardware execution.

## Current manuscript sources

`paper/manuscript.tex`, `paper/supporting_information.tex`, `paper/references.bib`,
and the CSV figure source data correspond to the current submission-v8 refinement.
Generated PDF figures are reproduced locally with `python paper/make_figures.py`.

## Manuscript figures

The manuscript and SI figures can be regenerated from the CSV source data in `paper/`
with:

```bash
python paper/make_figures.py
```

## Data integrity

The final 18q FCIDUMP SHA256 used by the validation campaign is:

```text
4fe2f74860c4a7b1af7e677f62f4d59a1571c2758ff427ba8627607dcd8c8536
```

The repository intentionally records failed or excluded branches instead of
silently converting them into positive results. In particular, the preliminary
Cu-S geometry scan failed its own CASSCF-convergence/orbital-continuity gates
and is not used for a geometry-dependent scientific claim.


## Submission revision v6

The current manuscript revision adds three deliberately critical checks:

- the exact geometry and archived ROHF/AVAS/CASSCF generation protocol for the
  canonical 18q Cu Hamiltonian;
- a null-space-basis sensitivity test showing that the ADAPT trajectory is not
  invariant under orthogonal rotations within multidimensional spin-preserving
  pool subspaces;
- a Suzuki diagnostic baseline showing that constituent noncommutativity is a
  useful structural separator, while optimized amplitude magnitude already
  carries most of the error ranking among noncommuting factors.

The dense exact target-space projector is used only as a classical validation
oracle in these small benchmarks; it is not proposed as a scalable quantum
primitive.  The expanded 24q FCIDUMP is retained as an exactly validated fixed
algorithmic Hamiltonian, but its dedicated 12-orbital active-space selection
driver was not preserved, so geometry-to-FCIDUMP reproducibility is not claimed
for that case.

Revision-specific scripts are under `src/revision_v6/`, benchmark metadata
under `benchmarks/`, and audit summaries/source data under
`results/revision_v6/` and `paper/`.


## Submission revision v8

The v8 referee audit supersedes the raw O(2)-mixing interpretation used in the
earlier v6 basis-sensitivity probe. The deterministic sparse null-space basis
is individually normalized but not mutually orthogonal, so raw pair mixing can
rescale generators and hence their ADAPT gradients.

The current manuscript uses a normalization-controlled audit instead: 100
random bases are evaluated at all 51 stored Cu prefix states under both
coefficient-vector and determinant-space Frobenius normalization. The
highest-gradient null-space group changes at only 4/51 prefixes under either
convention, with maximum changed-basis fractions 0.14 and 0.12. The old v6
rotation table is retained only for provenance and is not used for a scientific
claim in the v8 manuscript.
