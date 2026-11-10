# Sharp-Maxwell C5 read-only native source probe

## Purpose

This C5 deliverable implements a runnable vertical slice from immutable CompOSE source files through the native Fortran reader, the native Compose writer/read-back path, and the existing C1-C11 verifier. It is separate from the sharp-Maxwell constructor: it performs no physical source repair, normalization, projection, interpolation or matching.

The default cohort slice is SPG_M1 (Nadd=0) and GPPVA_FSU2 (Nadd=1). Both have source file blob identities pinned in the 65-pair cohort manifest. The Nadd=1 input exercises the provider-defined extra thermodynamic column without assigning it a physical meaning.

## CompOSE nuclear-group abundance semantics: provisional mapping

The previous version of this section made an incorrect semantic inference from Eq. (4.22): it described the stored \(Y_I\) as the baryon/mass fraction itself. That claim is superseded.

Under the standard CompOSE definitions:

- Eq. (3.10): \(Y_i=n_i/n_b\), the particle number fraction.
- Eq. (3.12): \(X_i=B_iY_i\), the particle mass-number (baryon) fraction.
- Eqs. (4.20)-(4.22): for a nuclear group, \(A_{\rm av}\) and \(Z_{\rm av}\) are abundance-weighted averages and \(Y_I=\sum_jY_j\) is the combined group number fraction.

Thus the standard baryon/mass fraction carried by the nuclear group is \(X_h=A_{\rm av}Y_I\). The native verifier's internal \`Yh\` is used as a baryon/mass fraction: C3 evaluates \(Y_p+(Z_{\rm av}/A_{\rm av})Y_h-Y_e-Y_\mu\), and C4 evaluates \(Y_n+Y_p+Y_h-1\).

The C5 reader/writer currently maps \`Yh = YI\` and \`YI = Yh\`. That mapping is retained as the committed raw-field-preservation path exercised by the read/write/read test. **It is not a validated physical conversion and is not endorsed by the CompOSE equations.** Round-trip equality proves serialization consistency only.

The C6 result exposed the unresolved provider-data issue:

- Both official data sheets identify \`Iphase=0\` as outer crust and \`Iphase=7\` as inner crust. The GPPVA(FSU2) sheet additionally identifies \`Iphase=6\` as homogeneous matter. The phase label \`7\` is therefore known; the unresolved issue is the physical meaning and consistency of the composition values in those records.
- An outer-crust example from each source has \((A_{\rm av},Z_{\rm av},Y_I,Y_e)=(56,26,1,0.46429)\). The direct mapping \`Yh=YI\` closes C3/C4 on these rows; the standard product \(A_{\rm av}Y_I=56\) cannot be a normalized mass fraction.
- At the first inner-crust row (row 218), SPG_M1 reports \((82.107,24.387,0.06087,0.06087)\), while GPPVA_FSU2 reports \((740,36,0.04865,0.04865)\), with the listed free-neutron, free-proton and muon fractions zero. The current native verifier fails C3/C4 on all inner-crust rows in both sources.

These observations do not establish a provider-specific replacement formula. In particular, neither applying \(A_{\rm av}\) blindly nor assuming that the field named \(Y_I\) is already a mass fraction resolves both sources and phases. The provenance of the provider's group-abundance values, possible source-table inconsistency or alignment issues, and the applicability of the verifier's composition closure must be reviewed against the original provider workflow.

No raw source file, parser mapping, verifier equation or tolerance is changed by this note. Until the evidence establishes a unique interpretation, the composition mapping remains provisional and physical admission remains blocked.

## Execution

From the repository root, with Python dependencies and gfortran available:

    python scripts/sharp_maxwell_source_probe.py \
      --manifest data/sharp_maxwell_65_cohort_manifest.v1.json \
      --source-id SPG_M1 \
      --source-id GPPVA_FSU2 \
      --output dist/sharp_maxwell_c5_source_probe.json

The runner calls compose_eos.api.read_compose, write_compose and check_consistency. It checks source Git blob SHA-1 and SHA-256 before reading, independently decodes the raw text, then performs native read -> write to a temporary directory -> read. The raw source directories are never used as output destinations.

## Preservation checks

- eos.nb: I_l, point count, and every density knot.
- eos.thermo: every Q1-Q7 value, row-wise Nadd, and row-wise q_extra.
- eos.compo: every Iphase, particle fraction, average A/Z, and the combined nuclear-group fraction Yh=YI.
- Manifest: point count, neutron/proton masses, Nadd, phase histogram/boundary, composition presence/row count and q_extra category histogram.
- Native round-trip: all EOS arrays including Iphase and q_extra, Nadd metadata and composition flags.

Decimal thermodynamic and composition arrays are compared with absolute tolerance 5e-13. q_extra is compared exactly after decimal-to-binary parsing. Raw source identity is exact by Git blob and SHA-256.

## Receipt semantics

The JSON receipt schema is sharp_maxwell_c5_native_source_probe_receipt.v1. Separate statuses are recorded for:
- source_identity_status: pinned original-file identity;
- raw_source_parse_status: independent read-only row decoder;
- reader_status: native Fortran reader execution;
- structural_status: row-wise and manifest preservation checks;
- compose_roundtrip_status: native writer then native reader;
- verifier_execution_status: native C1-C11 execution;
- scientific_status: native physical PASS/FAIL, not runner execution.

A scientific FAIL remains a scientific result and does not abort processing of the second source. The default exit code is zero only when source identity, native reading, structural checks, Compose round-trip and verifier execution all complete. The strict-physics option returns 3 if the scientific verifier fails; input, parser, runtime, structural or round-trip failures return 2.

The native verifier hard gate is C1-C6/C10-C11; C7-C9 are finite-difference diagnostics. It is a knot-level check, not whole-Hermite-interval stability or causality.

## Scope boundary

This probe does not establish full physical source admission. It does not determine provider Data Sheet semantics of Iphase or Nadd, perform whole-interval Hermite extrema, choose a density-reconstruction tolerance, align rest-mass references, construct the multi-segment phase-chain adapter, test crust-core pair compatibility, search Maxwell roots, or run TOV/Cowling.

All 65 pair records remain BLOCKED_BEFORE_ROOT_SEARCH; historical maxwell output remains LEGACY_NON_MAXWELL_READ_ONLY. No new sharp_maxwell result is produced.

## C5 implementation and C6 execution status

C5 committed the native reader/writer mapping and the read-only source probe, with six pure-helper regression tests. The exact committed test file and native integration path were not run at C5 time; that historical status is preserved in the manifest.

C6 later executed the official Fortran build and the exact committed test file on a Linux/gfortran host. The committed tests passed 6/6. Native read/write/read, verifier execution and round-trip completed for SPG_M1 and GPPVA_FSU2. The runner completed, but the scientific composition checks failed: C3 and C4 fail on every \(Iphase=7\) row (148/148 for SPG_M1 and 52/52 for GPPVA_FSU2), while every \(Iphase=0\) row passes. The C6 report and complete evidence are in [SHARP_MAXWELL_C6_NATIVE_EXECUTION_VERIFICATION_V1.md](SHARP_MAXWELL_C6_NATIVE_EXECUTION_VERIFICATION_V1.md) and \`artifacts/sharp_maxwell/C6_NATIVE_EXECUTION_20261110T093714Z/\`.

The C6 run was executed on commit \`43eea8bd9bee7ebb222e3c66babf5078b398d7ae\`; the branch subsequently advanced to include the C6 documentation and manifest updates. Do not describe the native run as having been executed on the later branch head.

Status: \`C6_PARTIAL\`; composition semantics \`BLOCKED_REVIEW_REQUIRED\`. This note supersedes the earlier assertion that CompOSE Eq. (4.22) directly establishes \`Yh=YI\` as a baryon/mass fraction. No C7, root search, pair sweep, TOV or Cowling run is authorized by this result.

