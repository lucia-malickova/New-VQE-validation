#!/usr/bin/env python3
"""Static pre-submission audit for the JCP manuscript source tree."""
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
main = (HERE / "manuscript.tex").read_text(encoding="utf-8")
si = (HERE / "supporting_information.tex").read_text(encoding="utf-8")
main_jcp = (HERE / "jcp_submission" / "manuscript_JCP.tex").read_text(encoding="utf-8")
si_jcp = (HERE / "jcp_submission" / "supplementary_material_JCP.tex").read_text(encoding="utf-8")
bib = (HERE / "references.bib").read_text(encoding="utf-8")
bib_jcp = (HERE / "jcp_submission" / "references.bib").read_text(encoding="utf-8")

errors = []
if main != main_jcp:
    errors.append("main manuscript mirror differs")
if si != si_jcp:
    errors.append("supplementary mirror differs")
if bib != bib_jcp:
    errors.append("bibliography mirror differs")

combined = main + "\n" + si
labels = re.findall(r"\\label\{([^}]+)\}", combined)
refs = re.findall(r"\\(?:ref|eqref)\{([^}]+)\}", combined)
for label in sorted(set(labels)):
    if labels.count(label) > 1:
        errors.append("duplicate label: " + label)
for ref in sorted(set(refs)):
    if ref not in labels:
        errors.append("unresolved reference: " + ref)

cites = []
for group in re.findall(r"\\cite\{([^}]+)\}", combined):
    cites.extend(x.strip() for x in group.split(","))
bibkeys = set(re.findall(r"@\w+\{([^,]+),", bib))
for cite in sorted(set(cites)):
    if cite not in bibkeys:
        errors.append("missing bibliography key: " + cite)

for directory, texfiles in [
    (HERE, [HERE / "manuscript.tex", HERE / "supporting_information.tex"]),
    (HERE / "jcp_submission", [HERE / "jcp_submission" / "manuscript_JCP.tex", HERE / "jcp_submission" / "supplementary_material_JCP.tex"]),
]:
    used = []
    for tex in texfiles:
        used.extend(re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex.read_text(encoding="utf-8")))
    for fig in used:
        if not (directory / fig).exists():
            errors.append(f"missing graphic in {directory.name}: {fig}")
    for pdf in directory.glob("fig*.pdf"):
        if pdf.name not in used:
            errors.append(f"unreferenced figure PDF in {directory.name}: {pdf.name}")

required_phrases = {
    "full Lucia affiliation": "Modelos Inteligencia Artificial S.L., Cl. Tajinaste 54, 386 20 San Miguel de Abona, Santa Cruz de Tenerife, Spain",
    "Lucia Euro-Q-Exa acknowledgement": "Euro-Q-Exa",
    "final M_S-preserving wording": "For the $M_S$-preserving generators considered here",
    "main CAS definition": "complete active space, CAS(15e,9o), comprising 15 active electrons in 9 active spatial orbitals",
    "main ROHF definition": "restricted open-shell Hartree--Fock (ROHF)",
    "main AVAS definition": "atomic-valence active-space (AVAS)",
    "main CASSCF definition": "complete-active-space self-consistent-field (CASSCF)",
    "main SLSQP definition": "Sequential Least Squares Programming (SLSQP)",
}
for label, phrase in required_phrases.items():
    if phrase not in main:
        errors.append("missing required content: " + label)
if "production physical pool uses $M_S$-preserving excitation patterns" not in si:
    errors.append("missing required content: SI M_S-preserving wording")
for label, phrase in {
    "SI ROHF definition": "restricted open-shell Hartree--Fock (ROHF)",
    "SI CAS definition": "primary complete active space, CAS(15e,9o), contains 15 active electrons in 9 active spatial orbitals",
    "SI AVAS definition": "atomic-valence active-space (AVAS)",
    "SI SLSQP definition": "Sequential Least Squares Programming (SLSQP)",
}.items():
    if phrase not in si:
        errors.append("missing required content: " + label)

if errors:
    print("SUBMISSION_AUDIT=FAIL")
    for err in errors:
        print(" -", err)
    sys.exit(1)

print("SUBMISSION_AUDIT=PASS")
print("MAIN_OBJECTS=", len(re.findall(r"\\begin\{(?:figure\*?|table)\}", main)), sep="")
print("SI_OBJECTS=", len(re.findall(r"\\begin\{(?:figure\*?|table)\}", si)), sep="")
print("CITATION_KEYS=", len(set(cites)), sep="")
