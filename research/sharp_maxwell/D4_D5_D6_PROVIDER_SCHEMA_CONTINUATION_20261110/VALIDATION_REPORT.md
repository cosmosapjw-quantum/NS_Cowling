# Validation report

Validation performed 2026-11-10.

| Check | Result | Qualification |
|---|---|---|
| D4 synthetic tests | PASS, 6 tests | Synthetic fixtures only. |
| D5 synthetic tests | PASS, 5 tests | Synthetic fixtures only. |
| Combined tests | PASS, 11 tests | Python pytest. |
| Python syntax/AST compilation | PASS | Both audit scripts compile. |
| Raw D4 source receipt | 4,975/4,975 parsed | Dropbox-extracted text; zero malformed rows and zero row-index mismatches. Raw files are not copied into the bundle. |
| D5 | DELTA REPLAY ONLY | 15 discrete-state boundary cells independently evaluated/excluded; all remaining classifications inherited from prior 5,550-cell baseline. |
| Native Fortran/project integration | NOT RUN | No usable local checkout/native executable. |
| D6 density drift/Maxwell roots | NOT RUN | Stop condition retained. |
| TOV/Cowling | NOT RUN | Downstream DAG not reached. |

Reproduce synthetic tests with PYTHONDONTWRITEBYTECODE=1 python -m pytest -q -p no:cacheprovider tests. The scripts and full JSON/CSV receipts are in the Dropbox backup bundle.