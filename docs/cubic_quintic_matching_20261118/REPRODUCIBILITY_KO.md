# Cubic/Quintic Hermite research: reproducibility and scope

Base source commit: 2b32672da2b9b602942a362ff7208510708be459
Branch: research/cubic-quintic-matching-modes-20261118
Requested density sources: 10 fixed RMF/core native model inputs (25,612 states) and 5 TF N128 droplet models (30 genuinely new stationary density holdouts).
Quintic derivative: PCHIP derivative of n(mu) from the sparse training set only; no withheld native states enter fit.
Actual native matching: 4 source pairs × cubic_F, quintic_F, cubic_P, quintic_P, CSS = 20 serialized and verified EOSs.
Maxwell: independent original-source root re-search on 4 pairs, 2 PK1 physical crossings and 2 SLy4 no-crossing cases.
Star metrics: 4 H sharp physical Maxwell EOS variants and 2 SPG→PK1 match pairs; mass 1.4 Msun, ℓ=2 equilibrium Cowling, slow/rapid interfaces.
P1: 32 bounded multi-root node-classification tasks and 5 targeted strict-accuracy refinements; unresolved cases are not called missing physical modes.
Stellar native quintic_F for two PK1 models is NOT whole-EOS admitted because normal smooth interval reconstruction fails, despite native matching readback success.
Numerical comparisons are either relative errors versus same-model full direct native EOS, or method shifts at identical physical mass; not observational errors.
Scope exclusions: global pasta geometry, frozen composition, elastic crust, full GR damping, complete mode spectrum and physical support beyond source.
Verification: exact full suite exit 0; 1508 passed, 9 skipped, 4 subtests passed; 228 source/test Python files same hashes.
Independent read-only audit: 403 checks, 0 fail.
Raw failures: dense v1 JSON serialisation, partial cancelled original full spectrum scan, blocked nodes, quintic_F whole-star rejection retained without fake PASS.

Principal files: data/cubic_quintic_matching_20261118/SUMMARY.json, INDEPENDENT_AUDIT.json, FULL_SUITE_RECEIPT.json, native_matching/CAMPAIGN.json, maxwell/CAMPAIGN.json, stellar_bounded/PARALLEL_RECEIPT.json, native_stellar/PARALLEL_RECEIPT.json, docs/cubic_quintic_matching_20261118/RESULTS_KO.md.
