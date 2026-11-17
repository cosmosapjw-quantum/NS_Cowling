# Phase continuation experiment — 2026-11-17

Base: 4a44333377fc6252ad6bcc0c1024614277199f31; remote CI #630 passed.
This is an unfinished local experiment, not an admission or publication.
Original checkout was not modified.

Using the retained SPG_M1 and full APR generated inputs from
`data/research_matching_20261117`, the baseline native `match-core` CLI was
called directly with CSS, 40 bridge points, and explicit composition dropping.
The core had phase labels 6, 8, 7. The command returned success and its 288-row
internal output had only Iphase 0. `baseline.log` preserves the raw response.
This bypassed the public frontend's existing SUPPORT guard solely to diagnose
native behavior; it is not an approved APR EOS or stellar input.

`experimental_guard.patch` records an attempted fail-closed check in the
native matcher and two direct-CLI tests. With the guard the same CLI returned
exit 3 (`guarded.log`). This guard is overbroad: the existing GPPVA_FSU2H
crust also has multiple labels. The focused test run produced 17 failures,
18 passes (`focused_tests.log`); failures include prior matching cases and
masked validation errors. The shell wrapper printed the report after pytest,
so its shell exit status is not the pytest result. Do not publish this patch.
The worktree retains this experimental diff for continuation.

Remaining: preserve source phase identities and their provenance through
Fortran assembly and CompOSE serialization, define allowed segment and
boundary handling from source support without inferring transition physics from
labels, and implement matching stellar intake. Re-run actual generation,
readback and independent physical checks before any commit/push and CI claim.
Source fidelity, independent physics, numerical support and CompOSE comparison
remain separate. The 18 native routes are not proof of full original EOS
reproduction.
