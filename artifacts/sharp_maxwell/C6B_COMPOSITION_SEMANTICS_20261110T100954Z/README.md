# C6B composition semantics evidence

Status: PROVIDER_SEMANTICS_UNRESOLVED__BLOCKED for both SPG_M1 and GPPVA_FSU2. No production code, verifier, tolerance, raw data, pair record or legacy output was changed.

Reproduce: `python scripts/audit_sharp_maxwell_composition_semantics.py --out-dir <DIR>` (deterministic; fails if a raw blob differs from the manifest pin). Tests: `python -m pytest -q python_frontend/tests/test_sharp_maxwell_composition_semantics_audit.py`.
Files: provider_provenance.md, hypothesis_adjudication.md, row_audit_*.csv/json (row-wise candidates and residuals), closure_summary.json (schema checks + per-phase aggregates), source_identity.json, repository_identity.json, execution_report.json. Verify: `sha256sum -c SHA256SUMS`.
Downloaded datasheet/paper PDFs are referenced by sha256 in provider_provenance.md, not copied here.
