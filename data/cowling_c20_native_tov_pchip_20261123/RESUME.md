# C20 scope and resume checkpoint — 2026-11-23 Asia/Seoul

Baseline: 122c502d1ef5e55fd890f47d86444df725987032 (C19).
Branch: research/cowling-c20-native-tov-pchip-20261123.
Objective: move native sharp TOV PCHIP slopes/dense output and batched
metric construction to Fortran, preserving the source-locked physical problem.
Python adds array marshaling, phase routing, test oracles and reporting only.
Existing Python source-contract EOS state evaluation remains a future migration.

Protected: original source archive SHA256
299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04;
Maxwell slow/rapid laws; cold Gamma1=Gamma_eq; DVODE tolerances; root/join
1e-5; node identity; phase overlap0.95; C17/C19 canonical evidence.
DVODE uses isolated OS processes, never shared concurrent threads.

Predeclared validation: native checked/trap build and O2, independent SciPy1.17.0
oracle for generic two-field interpolation, no-SciPy construction from original
M2/M3 TOV start, existing C13-C17 source-path regressions. Full original
40-density/4440-shot replay: exact branch statuses/node identities, root/join gates,
source/phase provenance and scan status/frequency/restart identity. Compare
root frequencies at abs<=0.005Hz (existing C13 scoped migration criterion),
mass/radius/interface geometry rel<=2e-7 (existing C15 native TOV criterion).
These are migration-equivalence limits, not an absolute physical accuracy claim.
Report actual deltas, raw scan changes and adaptive work counters.
1/2-process outputs of C20 must be byte-identical. Cold build prewarmed.

Progress: meaningful RED for original _NativePhase SciPy dependency observed,
and original M2/M3 no-SciPy TOV construction fails as intended. Frontend and native production implementation complete. Final integration:73 passed,8 skipped; C19/C20 serial/parallel science byte-identical. Publication state is recorded by the Git branch and backup HANDOFF.
Wolfram exact Hermite-vs-power identity returned zero, endpoint values/slopes
returned {y0,y1,d0,d1}, affine reproduction residual zero.

Native delegated tasks (PhysMath1.3 selected C1-C8 support, one integration review):
/root/fortran_migration_survey: read-only survey then owned new frontend test file.
/root/native_tov_pchip_impl: owns native generic PCHIP, TOV geometry module/tests.
/root/physical_metric_derivation: read-only same-background kinetic product derivation.
Parent owns frontend integration, scientific replay, review, backup and publication.

Final native geometry SHA256637d6f8514f55a01ccb90352395a352fe0052f95a4d5fc302f2bf677d17416d2.
Remaining: native original EOS state/grid migration and mode conditioning study.
