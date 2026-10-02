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

On Windows PowerShell use `$env:PYTHONPATH="$PWD\src"`.

The archived circuit-validation campaign reports Python 3.9.6, Qiskit 2.2.3, Qiskit Aer 0.17.2, Qiskit Nature 0.7.2, NumPy 2.0.2, and SciPy 1.13.1.

## Canonical Hamiltonians

18q:
```text
data/18q/active.FCIDUMP
SHA256 4fe2f74860c4a7b1af7e677f62f4d59a1571c2758ff427ba8627607dcd8c8536
```

24q:
```text
data/24q/active.FCIDUMP
```

The 24q file is a fixed algorithmic benchmark. Its dedicated geometry-to-12-orbital selection driver was not preserved.

## Materialize deterministic physical checkpoints

```bash
python src/materialize_checkpoint.py data/18q/best50_compact.json --fcidump data/18q/active.FCIDUMP --output best50_full.json
python src/materialize_checkpoint.py data/24q/best61_compact.json --fcidump data/24q/active.FCIDUMP --output best61_full.json
```

## Exact positive-control validation

```bash
python src/validate_physical_ansatz.py best50_full.json --fcidump data/18q/active.FCIDUMP --output reproduced_18q_exact.json
python src/validate_physical_ansatz.py best61_full.json --fcidump data/24q/active.FCIDUMP --output reproduced_24q_exact.json
```

Expected energy errors are approximately 1.492674 mEh and 1.181293 mEh, respectively, with `<S^2>=0.75` and unit doublet weight to numerical precision.

## Projected-to-parent representation audit

The canonical 121-factor Cu checkpoint is under:

```text
results/revision_v6/cu/canonical_projected_121_checkpoint_compact.json
results/revision_v6/cu/canonical_projected_121_audit_summary.json
```

Regenerate the state-specific certificate with:

```bash
PYTHONPATH="$PWD/src" python src/revision_v6/audit_state_specific_representation.py \
  results/revision_v6/cu/canonical_projected_121_checkpoint_compact.json \
  data/18q/active.FCIDUMP \
  --output reproduced_state_specific_representation_certificate.json
```

The committed source data for the final-angle prefix reconstruction are:

```text
paper/source_data_projected_negative_control_prefix_v5.csv
paper/source_data_state_specific_certificate_v5.csv
```

Important: manuscript prefix plots are **not ADAPT optimization-history curves**. Prefix `k` uses the first `k` factors and first `k` entries of the final globally optimized 121-parameter vector.

## Qiskit circuit validation

18q uniform Suzuki-2 scan:

```bash
for reps in 1 2 4; do
  python src/circuit_validate_physical.py best50_full.json \
    --fcidump data/18q/active.FCIDUMP \
    --method suzuki2 --reps "$reps" \
    --output "reproduced_18q_full_r$reps.json"
done
```

None of the tested uniform 18q schedules `r=1,2,4` passes all four strict validation criteria. The `r=4` result is the closest tested uniform comparator: fidelity 0.999484 and 0.135 mEh synthesis error, narrowly above the 0.1 mEh synthesis-error threshold.

Fixed-seed selective 18q validation:

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

Archived reference result:
- 19,660 CX gates;
- depth 23,016;
- synthesis error 0.0950208 mEh;
- circuit-to-fermionic fidelity 0.99958246;
- `<S^2>=0.75072418`;
- target-`M_S` weight 0.9999999999994.

## Figures

The figure source data are committed under `paper/source_data_*.csv`. Run:

```bash
python paper/make_figures.py
```

The current generator writes:
```text
paper/fig1_validation_hierarchy.pdf
paper/fig2_single_generator_scaling.pdf
paper/fig3_projected_prefix_divergence.pdf
paper/figS1_prefix_state_diagnostics.pdf
paper/figS2_hubbard_representation_audit.pdf
paper/figS3_local_suzuki_amplitude.pdf
paper/figS4_measurement_shots.pdf
paper/figS5_natural_occupation_errors.pdf
paper/figS6_multistart_robustness.pdf
```

## Independent benchmarks and additional audits

- NO and OH FCIDUMPs/metadata: `benchmarks/no/`, `benchmarks/oh/`.
- Pool-basis normalization/sensitivity data: `paper/source_data_pool_basis_*_v8.csv`.
- Fixed-parent energy/fidelity reoptimization data: `paper/source_data_fixed_parent_*.`
- Selective-Qiskit result: `paper/source_data_selective_suzuki_qiskit_v12.csv`.

The exact-projector machinery is a small-system classical validation oracle, not a proposed scalable circuit primitive.
