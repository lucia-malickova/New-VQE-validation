# Final JCP submission-support materials — v17

Canonical sources are `paper/manuscript.tex` and `paper/supporting_information.tex`; exact mirrored copies are retained in this directory as `manuscript_JCP.tex` and `supplementary_material_JCP.tex`.

REVTeX mode: AIP/JCP preprint.

Formal audit:
- main manuscript: 4/4 figure/table objects explicitly cited;
- supplementary material: 14/14 figure/table objects explicitly cited;
- bibliography keys resolved;
- no duplicate labels or unresolved `\\ref`/`\\eqref` targets;
- figure PDFs are generated deterministically by `paper/make_figures.py` and included in the submission bundle;
- main/JCP and SI/JCP mirrors are identical.

Scientific/reproducibility status:
- canonical 18q FCIDUMP provenance repaired and documented;
- fixed-seed selective 18q Qiskit validation PASS;
- final $M_S$-preserving/$[A,S^2]$ technical wording synchronized between main and SI.

Conflict of Interest and CRediT statements are inserted. Both authors must approve the final manuscript and individual CRediT roles before upload.
