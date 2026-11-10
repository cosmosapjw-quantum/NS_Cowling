# Sharp-Maxwell C5 read-only native source probe

## Purpose

This C5 deliverable implements a runnable vertical slice from immutable CompOSE source files through the native Fortran reader, the native Compose writer/read-back path, and the existing C1-C11 verifier. It is separate from the sharp-Maxwell constructor: it performs no physical source repair, normalization, projection, interpolation or matching.

The default cohort slice is SPG_M1 (Nadd=0) and GPPVA_FSU2 (Nadd=1). Both have source file blob identities pinned in the 65-pair cohort manifest. The Nadd=1 input exercises the provider-defined extra thermodynamic column without assigning it a physical meaning.

## CompOSE nuclear-group fraction correction

The input review exposed a concrete reader/writer defect in the current Fortran IO layer. CompOSE Reference Manual v3.01, Eq. (4.22), defines the final quantity in a nuclear-group quadruple as

    Y_I = (sum_j A_j Y_j) / A_av = sum_j Y_j.

This is already the combined group fraction. It is not A_av times the stored value. The source rows in these crust tables include A_av=56 and Y_I=1; multiplying them gives 56 and then fails the reader's own Y_h <= 1 bound. The internal Yh representation is the same combined baryon fraction used by C3/C4 composition closure.

The Fortran IO correction therefore makes the mapping symmetric:
- read: Yh = YI;
- write: YI = Yh.

Both old expressions (Yh=Abar*YI on read and YI=Yh/Abar on write) were inconsistent with the CompOSE definition. This correction does not change raw source files or infer the provider-specific identity of Iphase / Nadd.

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

## C5 implementation and execution status

The runner syntax check and six isolated helper tests passed. They cover Git blob identity, raw density and thermodynamic parsing, CompOSE group-fraction semantics, Iphase knot indexing, and refusal to coerce non-integer extra-column values. These helper tests do not execute the native EOS reader.

The actual native read/write/read/verifier run remains NOT_RUN in this checkpoint: the connected development device is offline and the current local runtime cannot fetch the repository source tree. No GitHub Actions result is claimed. C6 must execute this runner on a Linux/gfortran host and preserve its real receipt, classifying scientific FAIL separately from parser/runtime failure.
