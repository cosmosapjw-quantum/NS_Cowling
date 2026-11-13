# Sharp Maxwell C6C: provider clarification readiness (v1)

Date: 2026-11-10. Status: **C6C_READY_FOR_PROVIDER_QUERY** (preparation only). All messages: **DRAFT_NOT_SENT**.

## 1. State preserved

```text
C6: C6_PARTIAL
C6B: SPG_M1 PROVIDER_SEMANTICS_UNRESOLVED__BLOCKED; GPPVA_FSU2 PROVIDER_SEMANTICS_UNRESOLVED__BLOCKED
C6C: provider inquiry preparation only; messages DRAFT_NOT_SENT
C7: NOT_RUN / NOT_RECOMMENDED
65 candidate pairs: BLOCKED_BEFORE_ROOT_SEARCH x65
legacy maxwell: LEGACY_NON_MAXWELL_READ_ONLY
sharp_transition_search: NOT_RUN; tov_cowling: NOT_RUN
```

C6C does not resolve C6B. Nothing here authorizes a production change.

## 2. Identity and integrity (C6C-0/1)

- PR #4 open, draft, unmerged, no auto-merge; head at start `58eab171acd18fb85c7c86b1903487d42859d2d4`, base branch `research/maxwell-native-fortran-v1`; working tree clean.
- C6B evidence `SHA256SUMS` verified (11/11 OK); C5 probe + C6B audit tests: 16 passed; C6B audit rerun byte-identical (both row-audit CSV/JSON files and closure summary).
- Source authority commit `a09ec30d74b363652d2d1aefe990557c468679dd`; raw `eos.compo` blobs SPG_M1 `ab79da51...`, GPPVA_FSU2 `fa1ad2b8...` re-verified.

## 3. Public-source delta search (C6C-2)

No new primary source defines the fourth nuclear-group field, the group index, Aav/Zav scope or SPG gas recording. Details in `artifacts/sharp_maxwell/C6C_PROVIDER_CLARIFICATION_20261110T103036Z/public_source_delta_search.md`. Notable items:

- Data sheets and the SPG paper re-downloaded; hashes identical to C6B.
- NEW, RELATED only: GPPVA(FSU2) unified product (CompOSE /eos/215) repeats rows 1-269 identically and has 157 core rows (Iphase=6) with Aav=Zav=Y_I=0 and Yn, Yp nonzero.
- Zenodo record 11077759 (from a search summary) returned HTTP 403; not bypassed; UNVERIFIED lead.

## 4. Contact routes (C6C-3)

See `artifacts/sharp_maxwell/C6C_PROVIDER_CLARIFICATION_20261110T103036Z/provider_contacts.md`. GPPVA_FSU2: Helena Pais, `hpais@uc.pt`, VERIFIED (data sheet + CFisUC profile). SPG_M1: Luigi Scurto, data-sheet address `lscurto@student.uc.pt` (LISTED, possibly stale: CFisUC lists him as former PhD student) and CFisUC-listed `luigiscurto.pr@gmail.com` (LISTED); Helena Pais proposed as Cc (verified route, role in the export unconfirmed). CompOSE core-team address is shown only through an obfuscated mail link and was not transcribed. The person who ran each export is unresolved.

## 5. Dossiers and attachments (C6C-4..6)

Separate drafts for each source, with a shared background note, matrix and attachments, are in `artifacts/sharp_maxwell/C6C_PROVIDER_CLARIFICATION_20261110T103036Z/`. What only providers can resolve: the fourth-field meaning, Aav/Zav scope, SPG gas recording, whether Yn=Yp=0 is deliberate, and the baryon normalization. Builder: `scripts/build_sharp_maxwell_c6c_provider_package.py`; test: `python_frontend/tests/test_sharp_maxwell_c6c_provider_package.py`.

## 6. Constraints honoured

No change to `compose_eos_io.f90`, `compose_eos_verify.f90`, C3/C4, tolerances, raw data, the 65 pair records or legacy outputs; no `Yh=1`, no `Aav` scaling, no phase-dependent mapping; no root search, sweep, new sharp EOS, TOV/Cowling or C7. No message, e-mail, issue or form was sent or filed. Dropbox Paper not modified; C6B evidence not modified.
