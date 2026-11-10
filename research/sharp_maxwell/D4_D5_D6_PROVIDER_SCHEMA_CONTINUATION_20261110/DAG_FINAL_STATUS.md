# C7X DAG final status (2026-11-10)

| Node | Status | Finding |
|---|---|---|
| D0 | PASS_SCOPED_16_OF_18 | 16 selected sources; 46 raw-file Git-blob identities and row counts matched; two quarantined candidates remain excluded. |
| D1 | BLOCKED_NATIVE_EXECUTION | No usable local project checkout/native Fortran executable. |
| D2 | PASS_SCOPED_WITH_SOURCE_CAVEATS | Thermo/density alignment and identified phase boundaries checked; provider semantics separately gated. |
| D3 | PARTIAL_FAIL | Strict Q3/Euler fails 7/16; Q5 fails for RG_SLY4 and RG_KDE0v1. |
| D4 | BLOCKED_PROVIDER_OR_EXPORT_SCHEMA_MISMATCH | 13 eos.compo files, 4,975 rows parse; all 1,151 GPPVA/SPG inner-crust rows have Y_quad=Y_e and fail canonical and mass-diagnostic closure. No provider-confirmed alternate meaning. |
| D5 | DELTA_CORRECTED_EXPLORATORY_NOT_ADMISSION | 5,535 cells; 1,153 causality violations; 1,110 dn/dmu_E <= 0; zero nonpositive interpolated density; 4,382 selected interval passes. Boundary-delta correction only. |
| D6 | STOP_BEFORE_ROOT_SEARCH | 55/55 pair support intersections; 55/55 rest-mass rebasing required; 0/55 both-side D5 pass; D4 blocks all crust pairs. |
| D7-D10 | NOT_RUN | No adapter, pilot, sweep, TOV, or Cowling results. |

## Evidence classes

Directly observed: row/file structure, phase counts, species index positions, exact Y_quad=Y_e in 1,151 inner-crust rows, q_extra histograms, closure residuals from printed decimals.

Literature-supported: CompOSE particle fractions are Y_i=n_i/n_b; the group quadruple's fourth field is the combined group number fraction Y_I; the mass fraction is distinct. Model data sheets define model-specific phase/index and extra-field meanings.

Diagnostic only: GPPVA/SPG outer-crust Y_quad=1 and many RG heavy-nucleus rows exhibit mass-fraction-like closure. This is not permission to reinterpret Y_I.

Unresolved: export stage responsible for mismatch; exact group-field semantics in inner crust; RG_SLy4 row 143; RG_KDE0v1 q_extra; complete per-model Nadd maps where individual sheets were not inspected.

Not run: native Fortran, full raw-tree D5 replay, D6 density drift/common-pressure roots, D7-D10, provider inquiry.

**Stop rule:** do not advance to Maxwell root search before D4 clarification/corrected raw tables, D3 resolution, full D5 replay/admission, and native execution.