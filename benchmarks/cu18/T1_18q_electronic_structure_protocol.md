# Canonical 18q T1 Hamiltonian generation protocol

This records the archived geometry-to-FCIDUMP driver used for the canonical CAS(15e,9o) benchmark.

- Geometry: T1_model_geometry.xyz
- Total charge: +1
- Spin input: 2S = 1
- AO basis: STO-3G
- Reference: PySCF ROHF
- ROHF level shift: 0.2
- ROHF maximum cycles: 300; Newton fallback if needed
- AVAS labels: ["Cu 3d", "9 S 3p"]
- AVAS threshold: 0.2
- In the zero-based atom list, atom 9 is the short Cu-S donor
- AVAS result used by the canonical file: CAS(15e,9o)
- CASSCF maximum macro cycles: 50
- CASSCF maximum micro cycles: 20
- CASSCF energy tolerance: 1e-6
- CASSCF gradient tolerance: 1e-4
- FCI-solver tolerance: 1e-7
- Spin constraint: fix_spin_(ss=0.75)
- The archived generator explicitly permits using the best orbitals if the macro cap is reached. A separate convergence flag was not serialized with the FCIDUMP, so the file is treated as a fixed algorithmic Hamiltonian rather than as a chemically converged CASSCF prediction.

The expanded CAS(21e,12o) FCIDUMP is retained as a fixed 24q algorithmic benchmark from the same model campaign. The separate 12-orbital selection driver was not preserved, so geometry-to-FCIDUMP reproducibility is not claimed for the 24q case.
