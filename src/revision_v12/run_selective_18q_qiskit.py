#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import sys

import numpy as np
from scipy.sparse.linalg import expm_multiply

from spin_fci import (
    read_fcidump,
    pure_spin_hamiltonian,
    hf_determinant,
    spin_square_expectation,
)
from physical_pool import (
    build_pool,
    entry_from_dict,
    fermionic_terms_for_entry,
)

CANONICAL_FCIDUMP_SHA256 = "4fe2f74860c4a7b1af7e677f62f4d59a1571c2758ff427ba8627607dcd8c8536"
DEFAULT_SCHEDULE = {9: 4, 21: 2, 25: 2, 34: 2}
DEFAULT_SEED_TRANSPILER = 9272026


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def require_finite(label, *values):
    arr = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise RuntimeError("%s contains non-finite values: %r" % (label, values))


def materialize_generators(checkpoint, basis, U, splus, norb, nalpha, nbeta):
    if "selected_generators" in checkpoint:
        defs = checkpoint["selected_generators"]
        mode = "serialized_selected_generators"
    elif "selected_pool_indices" in checkpoint:
        pool = build_pool(
            basis, U, splus, norb, nalpha, nbeta, max_rank=2
        )
        indices = [int(x) for x in checkpoint["selected_pool_indices"]]
        if not indices:
            raise RuntimeError("selected_pool_indices is empty")
        if min(indices) < 0 or max(indices) >= len(pool):
            raise RuntimeError(
                "selected_pool_indices outside deterministic pool of size %d" % len(pool)
            )
        defs = [pool[i].as_dict() for i in indices]
        mode = "deterministically_materialized_from_selected_pool_indices"
    else:
        raise RuntimeError(
            "checkpoint must contain selected_generators or selected_pool_indices"
        )
    return defs, mode


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Validate the canonical 18q Cu nonuniform Suzuki-2 schedule "
            "with a fixed Qiskit transpiler seed."
        )
    )
    ap.add_argument("checkpoint")
    ap.add_argument("--fcidump", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--qpy-output", default=None)
    ap.add_argument(
        "--seed-transpiler", type=int, default=DEFAULT_SEED_TRANSPILER
    )
    ap.add_argument(
        "--optimization-level", type=int, default=1, choices=(0, 1, 2, 3)
    )
    ap.add_argument("--basis-gates", default="rz,sx,x,cx")
    ap.add_argument(
        "--expected-fcidump-sha256", default=CANONICAL_FCIDUMP_SHA256
    )
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    import qiskit
    import qiskit_aer
    import qiskit_nature
    from qiskit import QuantumCircuit, transpile
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.synthesis import SuzukiTrotter
    from qiskit_nature.second_q.mappers import JordanWignerMapper
    from qiskit_nature.second_q.operators import FermionicOp
    from qiskit_aer import AerSimulator

    actual_fcidump_sha = sha256(args.fcidump)
    if (
        args.expected_fcidump_sha256
        and actual_fcidump_sha.lower()
        != args.expected_fcidump_sha256.lower()
    ):
        raise RuntimeError(
            "FCIDUMP SHA256 mismatch: got %s expected %s"
            % (actual_fcidump_sha, args.expected_fcidump_sha256)
        )

    checkpoint = json.load(open(args.checkpoint, encoding="utf-8"))
    metadata = checkpoint.get("metadata", {})
    metadata_hash = metadata.get("fcidump_sha256")
    if metadata_hash and metadata_hash.lower() != actual_fcidump_sha.lower():
        raise RuntimeError(
            "checkpoint metadata FCIDUMP SHA256 %s does not match input %s"
            % (metadata_hash, actual_fcidump_sha)
        )

    d = read_fcidump(args.fcidump)
    nalpha = (d.nelec + d.ms2) // 2
    nbeta = d.nelec - nalpha
    ms = 0.5 * (nalpha - nbeta)

    basis, H, U, Hd, splus = pure_spin_hamiltonian(
        d, nalpha, nbeta, target_s=0.5
    )
    e0 = float(np.linalg.eigvalsh(Hd)[0])
    require_finite("exact target reference", e0)

    defs, checkpoint_mode = materialize_generators(
        checkpoint, basis, U, splus, d.norb, nalpha, nbeta
    )
    theta = np.asarray(checkpoint["theta"], dtype=float)
    if len(defs) != 50 or len(theta) != 50:
        raise RuntimeError(
            "expected canonical 50-application checkpoint; got %d generators and %d amplitudes"
            % (len(defs), len(theta))
        )
    require_finite("theta", *theta.tolist())

    entries = [entry_from_dict(x, basis) for x in defs]

    # Exact full-space fermionic target state for the stored positive-control ansatz.
    position = {int(x): i for i, x in enumerate(basis)}[
        hf_determinant(d.norb, nalpha, nbeta)
    ]
    psi = np.zeros(len(basis), dtype=complex)
    psi[position] = 1.0
    for t, entry in zip(theta, entries):
        psi = expm_multiply(float(t) * entry.matrix, psi)
    psi /= np.linalg.norm(psi)

    exact_energy = float(np.vdot(psi, H @ psi).real)
    exact_s2 = spin_square_expectation(psi, splus, ms=ms)
    require_finite("exact fermionic ansatz", exact_energy, exact_s2)

    nso = 2 * d.norb
    if nso != 18:
        raise RuntimeError("this validator is specific to the 18-qubit benchmark")

    hf_det = hf_determinant(d.norb, nalpha, nbeta)
    qc = QuantumCircuit(nso)
    for q in range(nso):
        if (hf_det >> q) & 1:
            qc.x(q)

    mapper = JordanWignerMapper()
    schedule = []
    total_pauli_terms = 0
    for one_based_index, (definition, t) in enumerate(zip(defs, theta), start=1):
        reps = int(DEFAULT_SCHEDULE.get(one_based_index, 1))
        schedule.append(reps)
        fop = FermionicOp(
            fermionic_terms_for_entry(definition),
            num_spin_orbitals=nso,
        )
        qop = mapper.map(fop)
        Hh = ((-1j) * qop).simplify(atol=1e-12)
        total_pauli_terms += len(Hh)
        qc.append(
            PauliEvolutionGate(
                Hh,
                time=-float(t),
                synthesis=SuzukiTrotter(order=2, reps=reps),
            ),
            range(nso),
        )

    basis_gates = [
        x.strip() for x in args.basis_gates.split(",") if x.strip()
    ]
    tqc = transpile(
        qc,
        basis_gates=basis_gates,
        optimization_level=args.optimization_level,
        seed_transpiler=args.seed_transpiler,
    )
    depth = int(tqc.depth())
    counts = {str(k): int(v) for k, v in tqc.count_ops().items()}
    two_qubit = int(
        sum(
            value
            for gate, value in counts.items()
            if gate in ("cx", "cz", "ecr", "iswap", "rzz", "rxx", "ryy")
        )
    )

    qpy_record = {"status": "not_requested"}
    if args.qpy_output:
        try:
            import qiskit.qpy as qpy

            with open(args.qpy_output, "wb") as fh:
                qpy.dump(tqc, fh)
            qpy_record = {
                "status": "saved",
                "path": args.qpy_output,
                "sha256": sha256(args.qpy_output),
            }
        except Exception as exc:
            qpy_record = {"status": "failed", "error": repr(exc)}

    simulator = AerSimulator(method="statevector")
    run = tqc.copy()
    run.save_statevector()
    statevector = np.asarray(
        simulator.run(run).result().get_statevector(run), dtype=complex
    )

    indices = np.asarray([int(x) for x in basis], dtype=np.int64)
    qms = statevector[indices]
    pms = float(np.vdot(qms, qms).real)

    if pms <= 0.0 or not np.isfinite(pms):
        raise RuntimeError("target-Ms weight is non-positive or non-finite")

    conditioned = qms / np.sqrt(pms)
    circuit_energy = float(np.vdot(conditioned, H @ conditioned).real)
    circuit_s2 = spin_square_expectation(conditioned, splus, ms=ms)
    doublet_coeff = U.conj().T @ conditioned
    p_doublet = float(np.vdot(doublet_coeff, doublet_coeff).real)
    fidelity_conditioned = float(abs(np.vdot(psi, conditioned)) ** 2)
    fidelity_full = float(abs(np.vdot(psi, qms)) ** 2)

    probabilities = np.abs(statevector) ** 2
    p_number = float(
        sum(
            probabilities[i]
            for i in range(len(probabilities))
            if bin(int(i)).count("1") == d.nelec
        )
    )

    synthesis_error_mEh = float((circuit_energy - exact_energy) * 1e3)
    require_finite(
        "compiled-circuit diagnostics",
        p_number,
        pms,
        circuit_energy,
        synthesis_error_mEh,
        circuit_s2,
        p_doublet,
        fidelity_full,
        fidelity_conditioned,
    )

    criteria = {
        "fidelity_full_min": 0.999,
        "abs_energy_error_mEh_max": 0.1,
        "abs_S2_minus_0p75_max": 1e-3,
        "target_Ms_weight_min": 0.999,
    }
    passed = bool(
        fidelity_full >= criteria["fidelity_full_min"]
        and abs(synthesis_error_mEh)
        <= criteria["abs_energy_error_mEh_max"]
        and abs(circuit_s2 - 0.75)
        <= criteria["abs_S2_minus_0p75_max"]
        and pms >= criteria["target_Ms_weight_min"]
    )

    output = {
        "purpose": (
            "18q Cu selective Suzuki-2 Qiskit re-transpilation "
            "and state-level validation"
        ),
        "checkpoint": args.checkpoint,
        "checkpoint_sha256": sha256(args.checkpoint),
        "checkpoint_mode": checkpoint_mode,
        "fcidump": args.fcidump,
        "fcidump_sha256": actual_fcidump_sha,
        "versions": {
            "python": sys.version,
            "qiskit": qiskit.__version__,
            "qiskit_aer": qiskit_aer.__version__,
            "qiskit_nature": qiskit_nature.__version__,
            "numpy": np.__version__,
        },
        "n_qubits": nso,
        "n_operators": len(defs),
        "selective_suzuki2_reps_1based": schedule,
        "nonuniform_entries_1based": DEFAULT_SCHEDULE,
        "seed_transpiler": int(args.seed_transpiler),
        "optimization_level": int(args.optimization_level),
        "basis_gates": basis_gates,
        "sum_pauli_terms_before_synthesis": int(total_pauli_terms),
        "transpiled_depth": depth,
        "transpiled_gate_counts": counts,
        "two_qubit_gate_count_proxy": two_qubit,
        "qpy": qpy_record,
        "exact_full_space_fermionic_ansatz": {
            "energy_Eh": exact_energy,
            "dev_vs_exact_doublet_mEh": float((exact_energy - e0) * 1e3),
            "S2": float(exact_s2),
        },
        "actual_transpiled_circuit": {
            "total_N_weight": p_number,
            "target_Ms_weight": pms,
            "energy_Eh_conditioned_on_target_Ms": circuit_energy,
            "energy_error_vs_exact_fermionic_ansatz_mEh": synthesis_error_mEh,
            "S2_conditioned_on_target_Ms": circuit_s2,
            "doublet_weight_within_target_Ms": p_doublet,
            "fidelity_full_with_exact_fermionic_ansatz": fidelity_full,
            "fidelity_conditioned_on_target_Ms": fidelity_conditioned,
        },
        "validation_criteria": criteria,
        "PASS": passed,
        "resource_note": (
            "Abstract gate-set resources only; no device coupling map "
            "and no physical-QPU claim."
        ),
    }

    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(output, fh, indent=2)

    print("SELECTIVE18Q_QISKIT_VERSION=%s" % qiskit.__version__)
    print("SELECTIVE18Q_SEED_TRANSPILER=%d" % args.seed_transpiler)
    print("SELECTIVE18Q_CX=%d" % counts.get("cx", 0))
    print("SELECTIVE18Q_DEPTH=%d" % depth)
    print(
        "SELECTIVE18Q_DELTA_E_MILLIHARTREE=%r"
        % synthesis_error_mEh
    )
    print("SELECTIVE18Q_FIDELITY_FULL=%r" % fidelity_full)
    print("SELECTIVE18Q_MS_WEIGHT=%r" % pms)
    print("SELECTIVE18Q_S2=%r" % circuit_s2)
    print("SELECTIVE18Q_DOUBLET_WEIGHT=%r" % p_doublet)
    print("SELECTIVE18Q_PASS=%s" % ("PASS" if passed else "FAIL"))

    if args.strict and not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
