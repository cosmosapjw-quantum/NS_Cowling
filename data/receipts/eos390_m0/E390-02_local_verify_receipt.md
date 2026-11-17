# E390-02 local native verify receipt (Haiku run, 2026-11-11)

Binary: gfortran 13.3.0 build of fortran_eos at a97582f base (build.sh), flags -std=f2018 -O2 -fcheck=all.
Command per model: compose_eos_fortran verify <raw eos_data/MODEL> <out.json>. Raw inputs read only.

| Model | exit | Result |
|---|---|---|
| APR_APR | 0 | FAIL C10 (known acausal tail, matches v3.2 baseline) |
| BBB_BHF2 | 0 | PASS |
| GM_GM1 | 0 | FAIL C3,C4 (not in v3.2 baseline notes; needs investigation) |
| GPPVA_* (8 models) | 3 | fail-closed: composition fractions outside physical range |
| SPG_M1..M5 (5 models) | 3 | fail-closed: composition fractions outside physical range |
| RG_KDE0v1, RG_SLY4 | 3 | fail-closed: thermo table is not a single selected T/Yq slice |

Unit tests: polynomial extrema PASS (exit 0); CSS segment analytic PASS (exit 0).

Interpretation: 12 of 18 raw providers are NOT readable by the current cold 1D reader.
Two causes: (a) multidimensional T/Yq layout (E390-10 scope); (b) composition-fraction
range rejection, consistent with the unresolved Iphase/Nadd semantics recorded in C6B.
Per DAG stop rule (unknown source-field semantics), E390-02 is NOT closed. No semantics were guessed.
