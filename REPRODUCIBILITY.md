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
