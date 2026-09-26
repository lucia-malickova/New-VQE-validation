# Audit of the two additional colleague ZIP files

## results_v2_for_professor.zip
This is numerically redundant with the final campaign for the final 18-qubit
workflow. Best-50 validation, generalized-pool audit, prefix-5 r=1, and the full
50-operator r=1/r=2/r=4 Qiskit results match the corresponding files in
results_final_campaign.zip.

## results_for_professor.zip
This is an earlier/intermediate campaign. It adds useful ablation data:
- restricted physical S+D: 35 ops, 153.628889 mEh;
- restricted physical S+D+T: 147 ops, 141.392349 mEh;
- restricted physical S+D+T with repeats: 250 applications, 13.166089 mEh,
  <S^2>=0.75, exact-doublet fidelity 0.878314, residual 0.041204 Eh;
- restricted 147-op Qiskit circuit at Suzuki-2 r=1: 218,502 CX, depth 261,514,
  synthesis error 0.01331 mEh, while the underlying fermionic state remains
  141.392 mEh above the target doublet.

The repeat-only 250-application checkpoint was independently revalidated against
the final 18q FCIDUMP before being included in the manuscript.

Important provenance issue:
The old package's adapt_physical_spinadapted_SD_validation.json used a different
FCIDUMP hash and reports a 4-operator/31.45 mEh result. It is inconsistent with
the 35-operator checkpoint and must not be used. Revalidation against the final
18q FCIDUMP gives 35 operators and 153.628889 mEh. The final campaign already
contains this corrected validation.

## Geometry scan
The final campaign geometry-generation validation is FAIL. Doublet and quartet
CASSCF did not converge across the five scan points and active-orbital continuity
failed the preset threshold. The scan is therefore excluded from the scientific
results, except as an explicitly stated failed validation extension.