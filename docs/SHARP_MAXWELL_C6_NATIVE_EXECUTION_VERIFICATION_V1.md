# Sharp Maxwell C6: Native Execution Verification (v1)

Status: **C6_PARTIAL**. Evidence: `artifacts/sharp_maxwell/C6_NATIVE_EXECUTION_20261110T093714Z/`.

## Scope
Execution check of the C5 read-only native probe on SPG_M1 and GPPVA_FSU2. Not a Maxwell certification, pair admission, root search, TOV or Cowling run.

## Identity
Branch research/sharp-maxwell-cohort-admission-v1, head 43eea8bd9bee7ebb222e3c66babf5078b398d7ae (base a09ec30d, C5 IO commit 5ecbca8e). All five C5 blobs matched the handoff. gfortran 13.3.0, Python 3.12.3 (differs from the C5 sandbox values; real values recorded).

## Results
| Gate | Result |
|---|---|
| Native build (official fortran_eos/build.sh) | PASS |
| Committed pytest, exact test file | PASS, 6/6 |
| SPG_M1 native probe | runner PASS; read/write/re-read PASS; verifier executed; scientific FAIL (C3, C4) |
| GPPVA_FSU2 native probe | same |
| Yh = Y_I mapping | PASS on raw rows |
| Raw source SHA-256 before/after (18 files) | identical |
| Legacy outputs (234 files) before/after | identical |
| 65 pairs | unchanged, all BLOCKED_BEFORE_ROOT_SEARCH |

## Open issue (BLOCKED_REVIEW_REQUIRED)
Native C3 (charge closure Yp + (Zbar/Abar) Yh = Ye + Ymu) and C4 (Yn + Yp + Yh = 1) fail on every Iphase=7 row (SPG_M1 148 of 365, max C4 residual 0.978; GPPVA_FSU2 52 of 269, 0.989) and pass on every Iphase=0 row. In those rows Y_I < 1 and Yn = Yp = 0 in the source. So the C5 statement that Yh = Y_I is the same baryon fraction used by C3/C4 holds only for Iphase=0 rows. Either the C3/C4 formulas assume the wrong quantity for Iphase=7 groups, or provider semantics of Iphase=7 rows are unresolved (already listed as open in C5). Per the C6 rules no physical definition, verifier or tolerance was changed. C7 entry is not recommended until this is reviewed.

Also: manifest source_catalog.composition_fields still reads Y_h=Abar*Y_I (pre-correction); left unmodified, compouse_group_fraction_mapping is authoritative.

## Not run
sharp_transition_search, tov_cowling, C7, C8, C9, C10, 65-pair sweep.

## Reproduce
```
bash fortran_eos/build.sh /tmp/c6build
export COMPOSE_EOS_FORTRAN_EXE=/tmp/c6build/compose_eos_fortran
python -m pytest -q python_frontend/tests/test_sharp_maxwell_source_probe.py
python scripts/sharp_maxwell_source_probe.py --manifest data/sharp_maxwell_65_cohort_manifest.v1.json --source-id SPG_M1 --source-id GPPVA_FSU2 --output out.json
```
