# Native sharp-Maxwell research DAG and next gates

Authority: `research/maxwell-native-fortran-v1` at verified head `d3649906292b0ec8f1caf98961f12436fa00871b`; PR #3 is open as draft and is not merged.

## Current dependency graph

```mermaid
flowchart TD
  A[Two one-sided EOS tables] --> B[Maxwell coexistence audit at Pt]
  B --> C[Core TOV to exact Pt]
  C --> D[Envelope TOV from shared rt mt nu_t]
  D --> E[Separate phase backgrounds and lapse normalization]
  E --> F[Core Cowling integration]
  F --> G[Slow or rapid interface transfer]
  G --> H[Envelope Cowling integration]
  H --> I[Single-frequency surface residual]
  I --> J[Frequency bracket scan and mode continuation]
  D --> K[Central-pressure sequence and M-R branches]
  K --> L[Radial stability and twin-star classification]
  J --> M[Benchmarked mode observables]
  L --> N[Hybrid-EOS parameter analysis]
  M --> N
```

## Gate ledger

| Node | Status | What is actually established |
|---|---|---|
| A-B: two-phase coexistence input | PASS_SYNTHETIC_TABLE | Two normalized T14 phase tables were loaded and admitted; both one-sided queries at Pt satisfy pressure/chemical-potential coexistence, positive jump ordering and cold first-law closure. Tables are generated analytic fixtures, not a real CompOSE hybrid pair. |
| C-E: phase-segmented TOV | PASS_SYNTHETIC_TABLE | Callback and normalized-table fixtures verify exact Pt endpoint/restart, shared geometric state and lapse normalization. Production CLI gives M=9.6357e-4 nominal solar masses, R=642.793 m, r_t=332.172 m for the artificial fixture only. |
| F-I: native Cowling interface shot | PASS_SYNTHETIC_TABLE | `ell=2`, one 1800 Hz trial, both slow and rapid maps, plus slow-conversion production CLI shot. Table-backed test pressure-perturbation continuity is below 1e-9; CLI sides agree at about 3.2e-16 relative. No spectrum search. |
| CLI / table provider | PASS_SYNTHETIC_TABLE | `maxwell_table_cowling` executed through the real `case_driver`; returned `status=ok`, `result_valid=true`, `cowling_shot_valid=true`. This remains a generated analytic two-table smoke fixture. |
| J: spectrum scan / mode continuation | NOT_STARTED | Requires bracket/root scan, node-count or invariant-based mode identity, and convergence tests. |
| K-L: M-R branches / radial stability | NOT_STARTED | Seidov ratio currently provides only a local onset diagnostic. |
| M-N: data benchmarks / inference | BLOCKED_BY_UPSTREAM | No realistic hybrid EOS benchmark or verified mode data cohort is attached to this implementation yet. |
| Real CompOSE high-/low-phase pair | NOT_STARTED | Requires provenance, units, crust/core matching, a source-backed Pt coexistence construction and full EOS-domain admission. |

## Next implementation order

1. Add negative provider/CLI tests for out-of-range Pt, mismatched chemical potential, reversed jump ordering and first-law closure failure. The valid two-table integration and production CLI path already pass.
2. Cross-check the synthetic table-backed native TOV result against the existing Python/SciPy phase-aware reference using identical pressure knots, constants, cutoff and central pressure. Compare interface radius/mass, total M/R and one-sided metric quantities.
3. Add tolerance-refinement tests for the phase-segmented TOV profile and Cowling state at the interface. Record a convergence envelope, not only one residual value.
4. Implement a frequency-bracket scanner around `shoot_maxwell_cowling`, keeping slow and rapid conversion spectra separate. Require bracket existence, stable root identity and refinement agreement before naming an `f`, `p1` or interface mode.
5. Add central-pressure continuation for M-R branches. Only then evaluate turning points and radial stability; keep the Seidov onset diagnostic as a separate local criterion.
6. Run a realistic hybrid EOS pair only after source provenance, units, crust/core matching, transition pressure and chemical-potential closure have been audited. Do not reuse the 234 smooth-matching spectrum rows as true Maxwell results.

## Stop conditions

- Stop before TOV integration if the two input phases do not satisfy coexistence at the declared `P_t`.
- Stop before interpreting a Cowling root if the background is not independently verified or the interface map fails its jump conditions.
- Do not call a single-frequency residual a mode or frequency spectrum.
- Do not infer stable branches from the Seidov ratio alone.
- Do not combine this branch's discontinuous phase results with the earlier smooth `maxwell` bridge cohort.
- Do not apply the 33-path v5 overlay without checking its recorded base SHA and all old-blob identities against the current branch.

## Current stop point

The native analytic and synthetic normalized-table implementations have passed exact-head tests and the production CLI smoke. The next unit of work is independent Python/SciPy cross-validation and negative table-admission tests. A source-backed realistic EOS cohort, full mode spectra and observational inference remain downstream.