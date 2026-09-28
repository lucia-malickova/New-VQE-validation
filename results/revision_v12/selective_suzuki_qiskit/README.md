# Submission-v12 selective 18q Qiskit validation

A colleague reran the nonuniform Suzuki-2 schedule for the canonical 18-qubit
Cu positive-control ansatz with Qiskit 2.2.3.

The run used:
- Python 3.9.6
- Qiskit 2.2.3
- Qiskit Aer 0.17.2
- Qiskit Nature 0.7.2
- NumPy 2.0.2
- basis gates {rz, sx, x, cx}
- optimization_level = 1
- seed_transpiler = 9272026
- r_9 = 4, r_21 = r_25 = r_34 = 2, all other applications r = 1

All four predeclared state-level gates pass:
- |delta E_synth| = 0.09502078273726511 mEh <= 0.1 mEh
- F(circuit, fermionic) = 0.9995824567353195 >= 0.999
- |<S^2>-0.75| = 0.0007241840360085 <= 0.001
- target-Ms weight = 0.9999999999993576 >= 0.999

Abstract compiled resources:
- CX = 19,660
- depth = 23,016

The received result ZIP had SHA256
cf78e4b11859a758bfeef97d9dc9b91681404e49a475143e76d5dc49324bf09d.
The archived transpiled QPY file had SHA256
ab917b9c9b91413a366707b7d33043655e50591166c8d47b90da222fd8a61c76.

## Reproducibility issue found during the run

The original helper package pinned repository commit 02c04e0, but that commit
contained a noncanonical data/18q/active.FCIDUMP and the helper script assumed
a fully serialized checkpoint although data/18q/best50_compact.json is
intentionally compact.

Submission-v12 repairs this by:
1. restoring the canonical 18q FCIDUMP,
2. adding the archived full serialized best-50 checkpoint,
3. adding a validator that accepts either the compact or full checkpoint,
4. restoring Python-3.9 compatibility, and
5. adding explicit finite-value checks.

The original colleague run was therefore a manual local reconstruction using
the canonical input files. The scientific PASS is supported by the recorded
input hashes and output diagnostics, while the repaired repository path should
be used for future reruns.
