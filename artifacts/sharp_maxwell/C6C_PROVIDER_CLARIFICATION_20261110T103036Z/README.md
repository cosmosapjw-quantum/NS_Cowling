# C6C provider clarification readiness (2026-11-10)

Status: **C6C_READY_FOR_PROVIDER_QUERY**. All inquiry messages are **DRAFT_NOT_SENT**. C6B remains `PROVIDER_SEMANTICS_UNRESOLVED__BLOCKED` for SPG_M1 and GPPVA_FSU2.

| File | Purpose |
|---|---|
| `evidence_summary_onepage.md` | One-page issue summary, equation sheet, hashes, aggregates, links (attachment) |
| `sample_rows_SPG_M1.csv`, `sample_rows_GPPVA_FSU2.csv` | Raw lines and audited fields for rows 1, 217, 218, a mid-inner row, and the last row (attachments) |
| `aggregate_statistics.json` | Per-Iphase aggregates with pinned blob SHA-1 / SHA-256 (attachment) |
| `public_source_delta_search.md` | New-source search record |
| `provider_contacts.md` | Contact route verification |
| `question_matrix.md`, `shared_background_note.md` | Matrix and shared note |
| `draft_SPG_M1_inquiry.md`, `draft_GPPVA_FSU2_inquiry.md` | Separate inquiry drafts, DRAFT_NOT_SENT |
| `repository_identity.json`, `SHA256SUMS` | Identity and integrity |

Reproduce: `python scripts/build_sharp_maxwell_c6c_provider_package.py --out-dir <dir>` (requires the authority commit in the local git object store and the C6B evidence directory; output is deterministic and byte-identical to `sample_rows_*.csv` and `aggregate_statistics.json` here). Full row-level audit stays in the C6B directory (not duplicated).

No production code, verifier, tolerance, raw data, pair record or legacy output was changed. No inquiry has been sent.
