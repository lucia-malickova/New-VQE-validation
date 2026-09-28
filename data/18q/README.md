# 18-qubit canonical Cu benchmark inputs

This directory contains two forms of the same validated 50-application
representation-faithful ansatz.

- `active.FCIDUMP` is the canonical CAS(15e,9o) Hamiltonian used by the
  validated 18q campaign. Its SHA256 is
  `4fe2f74860c4a7b1af7e677f62f4d59a1571c2758ff427ba8627607dcd8c8536`.
- `best50_compact.json` is intentionally compact: it stores the 50
  `selected_pool_indices` and optimized amplitudes, but not serialized
  generator definitions. Use `src/materialize_checkpoint.py` or the
  submission-v12 selective-Qiskit validator to reconstruct those generators
  deterministically from the committed pool construction.
- `best50_generalized_spin_adapted.json` is the full serialized checkpoint
  used in the archived validation campaign, including `selected_generators`.
  Its SHA256 is
  `594214e506c7e2411e1b3e35d9419ad52bf95cc8d60a0ca4e08292f0ba04e4db`.

The compact and full checkpoints contain the same selected pool indices and
the same 50 optimized amplitudes. The full file is retained as an independent
provenance/reference artifact; the compact form remains useful for testing the
deterministic materialization path.
