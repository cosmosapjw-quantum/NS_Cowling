# Validation matrix

| Gate | Status | Evidence |
|---|---|---|
| Software | PASS | frontend/pipeline 1242 passed, 9 skipped; generator 142 passed; theory 6 passed; four native controls exit 0 |
| Scientific | PASS_SCOPED | One CSS cell passed first-law, stability, causality, jump and admission checks; APR boundary law remains HOLD |
| Numerical | PASS_SCOPED | 252-row readback exact; analytic CSS max relative error 4.35e-15; first-law integral max absolute residual 2.80e-13 MeV |
| Operational | PASS_SCOPED | One allowed run consumed; requested/eligible/attempted/admitted = 390/1/1/1 |
| Reproducibility | PASS | Manifest, claim, executable, execution, census, consistency and admission hashes are in `CLOSEOUT.json` |
| Independent review | PASS_SCOPED | Post-execution raw-artifact review; stable-star/Cowling and APR consumer claims remain unapproved |
