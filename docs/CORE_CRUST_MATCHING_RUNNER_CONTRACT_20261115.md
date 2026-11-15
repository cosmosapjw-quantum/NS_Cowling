# Core--crust matching runner contract and readback audit

**Repository:** cosmosapjw-quantum/NS_Cowling  
**Audit baseline:** 69e9752023e3f6c93552e99f0c5e455560b5a568  
**Scope:** `scripts/match_core_crust.py`, native Fortran matching/write/read path, and survivor manifest.

## Physical contracts

The runner contains two distinct constructions; they must not be described as interchangeable.

### Sharp Maxwell

The low- and high-density phases coexist at common pressure and baryon chemical potential. Baryon density and energy density may have a finite jump, subject to the cold first-law jump identity. The construction is local-neutrality Maxwell, not a Gibbs mixed phase.

### Single-speed CSS bridge

One constant (c_s^2) is used over the CSS bridge. At the crust--CSS interface, pressure and (mu_B) are continuous while (n_B) and (epsilon) have a positive jump satisfying
[
Deltaepsilon=mu_tDelta n_B.
]
At the CSS--core endpoint, (n_B,P,epsilon,mu_B) are continuous. This construction is **not** a density-continuous join over the entire crust-plus-bridge domain.

The runner is thermo-only: composition is dropped. The construction contract does not imply C1 continuity of (c_s^2), nor does it automatically establish perturbation-interface validity.

## Survivor admission contract

A generated EOS is a survivor only if all of the following are true:

1. Native matching returns status `ok` and a plausible point count.
2. Native CompOSE writing returns status `ok`.
3. Native CompOSE readback returns status `ok`.
4. Native-reported readback row count equals the matched row count.
5. The pre-write and post-read internal CSV files both exist and parse without malformed rows or duplicate columns.
6. Their column sets agree and include `nb`, `Q1`--`Q7`, and `Iphase`.
7. Every numerical column is compared row-by-row. This includes optional q_extra columns, composition payload values, row status, presence masks, and C3/C4 applicability flags.
8. Integer CSV metadata columns match exactly; finite floating-point values match within the declared `rtol=5e-13`, `atol=1e-14` criterion.
9. The complete scalar `eos_data.nml` metadata set is compared as well: rest-mass conventions, grid/index metadata, temperature/charge-grid values, composition-presence flags, and beta-equilibrium flag. Integers and logicals match exactly; floating-point metadata obeys the same tolerance.

A status-only or row-count-only readback is insufficient. The report contains per-column maximum absolute and relative residuals. Any discrepancy rejects that pair and prevents it from entering the survivor list.

## Provenance

For surviving pairs the manifest records:

- source CompOSE file SHA256 hashes for crust and core;
- output CompOSE file SHA256 hashes;
- native executable SHA256;
- native match and write reports;
- readback validation result and residuals;
- UTC build time and schema version.

The hashes establish byte identity only. They do not establish the correctness of the source physics.

## Validation status

The tests added in `pipeline/tests/test_core_crust_matching_pipeline.py` are intended to demonstrate that:
- a single changed (Q_i) fails readback;
- a point-count mismatch fails readback;
- the runner cannot classify a mismatch as a survivor.

The CI run on test-only commit `69e9752023e3f6c93552e99f0c5e455560b5a568` reproduced these three expected failures (269 passed, 9 skipped, 3 failed). The implementation commit must rerun the exact PR head and record its exact SHA and CI result before claiming verification.
