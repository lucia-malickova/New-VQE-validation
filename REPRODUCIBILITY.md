# Reproducibility guide

Run commands from the repository root.

## Environment

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
export PYTHONPATH="$PWD/src"
```

On Windows PowerShell use `$env:PYTHONPATH="$PWD\\src"`.
The colleague's stored environment record reports Python 3.9.6.

## Materialize the deterministic physical checkpoints

The repository stores compact checkpoints as pool indices plus optimized angles.
Materialize the exact serialized physical generators first:

```bash
python src/materialize_checkpoint.py data/18q/best50_compact.json --fcidump data/18q/active.FCIDUMP --output best50_full.json
python src/materialize_checkpoint.py data/24q/best61_compact.json --fcidump data/24q/active.FCIDUMP --output best61_full.json
```

The pool construction is deterministic. It was checked against the original
serialized campaign checkpoints: all selected generator specifications and
coefficients are reproduced exactly for both 18q and 24q.

## Exact physical ansatz validation

```bash
python src/validate_physical_ansatz.py best50_full.json --fcidump data/18q/active.FCIDUMP --output reproduced_18q_exact.json
python src/validate_physical_ansatz.py best61_full.json --fcidump data/24q/active.FCIDUMP --output reproduced_24q_exact.json
```

Expected errors are approximately 1.492674 mEh and 1.181293 mEh, respectively,
with <S^2>=0.75 and unit doublet weight to numerical precision.

## Qiskit circuit validation

18q smoke test:

```bash
python src/circuit_validate_physical.py best50_full.json --fcidump data/18q/active.FCIDUMP --prefix 5 --method suzuki2 --reps 1 --strict --output reproduced_18q_prefix5.json
```

18q full synthesis scan:

```bash
for reps in 1 2 4; do
  python src/circuit_validate_physical.py best50_full.json --fcidump data/18q/active.FCIDUMP --method suzuki2 --reps "$reps" --output "reproduced_18q_full_r${reps}.json"
done
```

The r=4 result has fidelity 0.999484 and 0.135 mEh synthesis error, but it does
not pass the deliberately strict 0.1 mEh synthesis-error criterion.

The submission-v12 selective 18q circuit can be reproduced directly from the
compact checkpoint; the validator materializes its generators deterministically
when `selected_generators` are absent:

```bash
python src/revision_v12/run_selective_18q_qiskit.py \
  data/18q/best50_compact.json \
  --fcidump data/18q/active.FCIDUMP \
  --seed-transpiler 9272026 \
  --optimization-level 1 \
  --basis-gates rz,sx,x,cx \
  --qpy-output reproduced_18q_selective.qpy \
  --strict \
  --output reproduced_18q_selective.json
```

The archived reference run used Qiskit 2.2.3, Qiskit Aer 0.17.2, Qiskit Nature
0.7.2, NumPy 2.0.2, SciPy 1.13.1, and Python 3.9.6. It gives 19,660 CX gates,
depth 23,016, 0.0950208 mEh synthesis error, fidelity 0.99958246,
$\langle S^2\rangle=0.75072418$, and target-$M_S$ weight
0.9999999999994. The result passes all four predeclared circuit-equivalence
gates. The full serialized checkpoint is also retained at
`data/18q/best50_generalized_spin_adapted.json`.

The canonical 18q FCIDUMP must have SHA256
`4fe2f74860c4a7b1af7e677f62f4d59a1571c2758ff427ba8627607dcd8c8536`.
The revision-v12 validator checks this hash by default. Its Python-popcount path
uses `bin(i).count("1")`, so the recorded Python 3.9.6 environment is supported.

24q full strict validation:

```bash
python src/circuit_validate_physical.py best61_full.json --fcidump data/24q/active.FCIDUMP --method suzuki2 --reps 1 --strict --output reproduced_24q_full_r1.json
```

The stored reference run passes the strict circuit-equivalence gate.

Resource-only scans can use:

```bash
python src/circuit_resource_physical.py best61_full.json --fcidump data/24q/active.FCIDUMP --method suzuki2 --reps 2 --output reproduced_24q_resources_r2.json
```

These are abstract-basis resources unless transpilation is replaced by an actual
device backend and coupling map.

## Figures

```bash
python paper/make_figures.py
```

## Excluded extension

The preliminary Cu-S scan failed its own CASSCF-convergence and orbital-continuity
gates. Its log is retained under `results/excluded/`, but it is not used for a
geometry-dependent scientific claim.


## Additional referee checks

Canonical Cu geometry and the archived 18q generation protocol are stored in:

```text
benchmarks/cu18/T1_model_geometry.xyz
benchmarks/cu18/T1_18q_electronic_structure_protocol.md
```

The generalized-pool basis-sensitivity test can be rerun with:

```bash
PYTHONPATH="$PWD/src" python src/revision_v6/test_pool_basis_sensitivity_v6.py   --fcidump data/18q/active.FCIDUMP --seeds 1 2 3   --output reproduced_pool_basis_sensitivity_v6.json
```

The Suzuki diagnostic baseline can be regenerated with:

```bash
python src/revision_v6/compare_suzuki_diagnostics_v6.py   paper/source_data_local_suzuki_predictor_v5.csv   --output reproduced_suzuki_diagnostic_v6.csv
```

The independent NO and OH FCIDUMPs and metadata used in the manuscript are under
`benchmarks/no/` and `benchmarks/oh/`.

The exact projector used in the representation audit is a classical small-system
oracle, not a proposed scalable circuit primitive.


## Normalization-controlled basis audit

The v6 raw O(2)-mixing probe is retained for provenance but is superseded for
scientific interpretation because the canonical two-dimensional null-space
basis vectors are individually normalized but not mutually orthogonal.

Use the normalization-controlled audit instead:

```bash
PYTHONPATH="$PWD/src" python src/revision_v8/audit_normalized_pool_basis_sensitivity_v8.py \
  --fcidump data/18q/active.FCIDUMP \
  --checkpoint data/18q/best50_compact.json \
  --n-random 100 \
  --output-prefix reproduced_pool_basis_normalized_v8
```

The committed summary data are under paper/source_data_pool_basis_*_v8.csv
and results/revision_v8/cu/.
