# Submission revision v12 — referee hardening and circuit closure

Revision v12 is a referee-hardening and reproducibility revision for the planned
JCP submission. It does not change the central representation-equivalence
result. It sharpens novelty positioning, clarifies the role of the
spin-preserving pool as a positive control, adds a four-layer validation
schematic, and closes the previously incomplete 18q circuit-synthesis layer.

## Manuscript changes

1. **Novelty positioning.** The Introduction separates the present
   representation-equivalence audit from symmetry projection / symmetry-adapted
   VQE, symmetry-preserving state-preparation circuits and symmetry-aware pools,
   contextual-subspace VQE, and exact or closed-form spin-adapted factorization.
2. **Operational failure path.** The Introduction states explicitly how a
   reduced-space optimization can be incorrectly identified with a parent
   full-space implementation when labels and amplitudes are transferred without
   an invariance audit.
3. **Positive-control scope.** The generalized spin-preserving pool is explicitly
   a representation-faithful positive control, not a claimed new scalable
   spin-adapted ADAPT algorithm.
4. **Literature coverage.** Targeted references were added for
   symmetry-preserving state preparation and symmetry-projected VQE.
5. **Conceptual figure.** The main manuscript now contains a four-layer
   validation ladder:
   reduced-space optimization -> parent full-space sequence -> synthesized
   qubit circuit -> measurement protocol.
6. **Cover letter.** The novelty distinction is made explicitly without claiming
   that prior projected or spin-adapted methods are incorrect.

## Confirmed 18q selective Suzuki-2 circuit

The determinant-equivalent emulator identified the nonuniform schedule

- application 9: r=4
- applications 21, 25, 34: r=2
- all other applications: r=1

An independent rerun re-transpiled exactly this schedule with Qiskit 2.2.3,
Qiskit Aer 0.17.2, Qiskit Nature 0.7.2, optimization level 1, abstract basis
{rz,sx,x,cx}, and seed_transpiler=9272026.

Verified result:

- |delta E_synth| = 0.09502078273726511 mEh
- full circuit-to-fermionic fidelity = 0.9995824567353195
- <S^2> = 0.7507241840360085
- target-Ms weight = 0.9999999999993576
- doublet weight within target Ms = 0.9997586053213303
- CX = 19,660
- depth = 23,016
- strict circuit-equivalence gate: **PASS**

The corresponding uniform-r=4 reference has 61,482 CX, depth 71,636, and
0.1354618902951188 mEh synthesis error, so it fails the strict energy gate.
The selective compiled realization uses about 68.0% fewer CX gates and 67.9%
less depth than that uniform-r=4 realization while passing all four state-level
criteria.

The received validation ZIP SHA256 is
`cf78e4b11859a758bfeef97d9dc9b91681404e49a475143e76d5dc49324bf09d`.
The returned QPY SHA256 is
`ab917b9c9b91413a366707b7d33043655e50591166c8d47b90da222fd8a61c76`.

## Repository defects found and repaired

The independent rerun exposed three genuine reproducibility problems.

1. The committed `data/18q/active.FCIDUMP` did not match the canonical
   Hamiltonian used by the validated 18q campaign. Revision v12 restores the
   canonical file with SHA256
   `4fe2f74860c4a7b1af7e677f62f4d59a1571c2758ff427ba8627607dcd8c8536`.
2. `best50_compact.json` is intentionally a compact checkpoint containing
   pool indices and amplitudes rather than serialized generator definitions.
   The failed helper incorrectly assumed a `selected_generators` field.
   Revision v12 keeps the compact form, adds the archived full serialized
   checkpoint `best50_generalized_spin_adapted.json`, and makes the new
   selective validator accept either form.
3. The helper used `int.bit_count()`, which is incompatible with the archived
   Python 3.9.6 environment. The rerun path now uses a Python-3.9-compatible
   population count and includes explicit finite-value checks.

A direct audit confirms that the compact and full checkpoints have identical
50 selected pool indices and identical optimized amplitudes, and that the
restored FCIDUMP and full checkpoint exactly match the archived validated
inputs.

The colleague's macOS run emitted BLAS-related RuntimeWarnings during some
matrix multiplications. Explicit finite-value checks confirm that the final
energies, spin diagnostics, overlaps, and weights used for acceptance are all
finite; the warnings are therefore recorded as an environment diagnostic, not
silently discarded.

## Scope

The compiled resource counts remain abstract-basis results. No device coupling
map or calibration data are used, so no hardware-resource or physical-QPU claim
is made. The selective schedule is the minimum-proxy passing schedule within
the explicitly searched grid; no global resource-optimality claim is made.
