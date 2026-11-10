# R16 cellwise-minimal Hermite repair: matched-pressure TOV pilot

**Status:** `PASS_R16_MINIMAL_REPAIR_TOV_PILOT`  
**Workflow:** [GitHub Actions run 42752885366](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42752885366)  
**Code commit:** `5fc8abe69c9ed8eeab51fcf6ee05068b4b176771`  
**Frozen spectrum-summary SHA-256:** `ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965`  
**R15 reference-result SHA-256:** `0b157864698af4756576c85e856bf4acc708ad2e759be337f1ca97be1089e962`

## 1. Question and comparison contract

This bounded pilot compares the cellwise-minimal SI-native Hermite pressure repair against the frozen R15 q-linear projected solver input. It uses three crossed EOS pairs:
- GPPVA_FSU2 × BBB_BHF2 (source-family/structural stress pair);
- SPG_M2 × BBB_BHF2 (crust-family contrast at fixed core);
- GPPVA_FSU2 × GM_GM1 (core-family contrast at fixed crust).

Every pair has the six historical smooth bridge tokens and eight native source-row markers `(1, 9, 17, 23, 52, 73, 116, 217)`: 3 × 6 × 8 = 144 integrations.

The central pressure (P_c) is taken from the same frozen 234-table spectrum summary for each EOS name. The absolute cutoff pressure (P_{\rm cut}) is taken from the corresponding R15 case row and is passed unchanged to the repaired-EOS TOV solve. The comparison therefore matches both (P_c) and (P_{\rm cut}) exactly; the source-row index is only a density marker. Since repaired and projected pressures at the same source row need not agree, the cutoff in the repaired table is generally an interpolated isobar rather than that exact native knot.

The same frozen v5 native Fortran `cell_solve task=star` backend is used in both runs. All 144 case records report the same backend build key as R15: `a65563fabee0e662430050b6`.

## 2. Derived EOS construction and source preservation

Each of the 18 normalized-source tables was copied into a separate derived diagnostic CSV. The SI-native repair clips only source pressure increments outside the exact cubic-Hermite admissible interval:
[
L_i < \Delta P_i < U_i,
]
with a (10^{-8}) relative interior-width margin; increments already accepted by the native constructor are preserved exactly. Pressure is accumulated from the original first pressure knot and energy density is reconstructed as (epsilon=\mu_b n_b-P). No row was added, removed or reordered.

Across the 18 derived tables:
- 8,034 native cells total;
- 1,586 pressure increments adjusted;
- 6,448 pressure increments preserved exactly;
- zero synthetic rows;
- all 18 derived-table SHA-256 values match the receipts;
- all 18 normalized-source and all 18 projected-reference files remained unchanged.

The procedure is cellwise-minimal in the sense of nearest admissible pressure increment for each cell. It is **not** a globally optimal variational reconstruction of the full EOS.

## 3. TOV acceptance and independent artifact verification

| Acceptance item | Result |
|---|---:|
| Derived EOS tables | 18 / 18 |
| TOV cases | 144 / 144, all `ok` |
| Background profiles | 144 / 144 |
| Profile SHA-256 and data-row count | 144 / 144 |
| Unique case IDs | 144 / 144 |
| Exact absolute (P_{\rm cut}) match to R15 | 144 / 144 |
| Original normalized source mutations | 0 |
| Frozen projected-reference mutations | 0 |
| Pilot receipt failures | 0 |
| Radial stability / Cowling modes | Not computed |

The downloaded 56.5-MB artifact passed ZIP integrity checking. Its SHA-256 is `d8acfd1dad803e6580c97c2fcd1084d6c002ddd043abd6883fdad143d51717bc`, equal to the GitHub artifact digest. The manifest receipt, cases CSV, tables CSV and pair-cutoff summary CSV hashes all independently match. Every profile file's SHA and data-row count was recomputed independently.

Artifact ID: `13734297622`.

## 4. Matched-pressure structural response

Percentages below compare repaired EOS values to the corresponding frozen R15 projected-EOS value at the same (P_c) and same absolute (P_{\rm cut}).

| Pair | Mean (Delta M/M_{\rm R15}) | Mean (Delta r_{\rm cut}/r_{\rm R15}) |
|---|---:|---:|
| GPPVA_FSU2 × BBB_BHF2 | -0.263900% | -0.335693% |
| SPG_M2 × BBB_BHF2 | -0.256384% | -0.335217% |
| GPPVA_FSU2 × GM_GM1 | -0.006094% | -0.000610% |
| All 144 cases | -0.175459% | -0.223840% |

Across all cases, the mass response ranges from (-0.264961\%) to (-0.005612\%), and the internal-isobar radius response from (-0.337241\%) to (-0.000220\%). Across the eight cutoffs, the cohort-mean mass response is nearly constant (about (-0.175461\%) to (-0.175450\%)); the cohort-mean radius response ranges from (-0.224467\%) to (-0.221554\%).

The result is strongly core-conditioned in this small crossed pilot: the two BBB_BHF2 pairs show (sim0.26\%) mass and (sim0.335\%) internal-radius shifts, while the GM_GM1 pair shows (sim0.0061\%) and (sim0.00061\%), respectively. Method-to-method differences are considerably smaller than this pilot's core contrast.

The largest mismatch between the repaired table's pressure at a selected source row and the common R15 cutoff pressure is (0.616596\%). This is expected from the matched-absolute-isobar design: source-row labels locate the low-density sampling sequence, but (P_{\rm cut}) itself remains the frozen R15 pressure, not the repaired table's pressure at the same row index.

## 5. What the TOV contrast measures

This result must **not** be interpreted as the effect of correcting only the rejected low-density cells. The comparison is between two complete EOS pressure-closure policies over the whole table:

- the current q-linear projector replaces every source pressure increment by the exact integral of piecewise-linear (q=n_b/\mu_b), including cells whose source pressure increment already passes the native Hermite test;
- the minimal repair preserves source pressure increments in already-admissible cells and modifies only increments rejected by the native Hermite constructor.

The difference is especially visible near the central-pressure region in representative CSS tables (nearest source row to (P_c)):
- GPPVA_FSU2 × BBB_BHF2: q-linear (P) is about (0.6190\%) below the minimal-repair (P), while q-linear (epsilon) is about (0.09255\%) higher;
- SPG_M2 × BBB_BHF2: q-linear (P) is about (0.4941\%) below minimal-repair (P), while q-linear (epsilon) is about (0.09584\%) higher;
- GPPVA_FSU2 × GM_GM1: the corresponding differences are only about (+0.000785\%) in (P) and (-0.0000845\%) in (epsilon).

These are representative endpoint comparisons, not a complete central-density sensitivity analysis. They explain why the TOV response differs by core family, and they show that the pilot probes a global thermodynamic reconstruction choice, not solely an outer-crust patch.

## 6. Interpretation boundary

The reported (M) and (r) are internal finite-pressure isobar values (M(r_{\rm cut})) and (r_{\rm cut}). They are not complete-star (P=0) mass/radius differences, physical crust mass/thickness, or surface-radius corrections. The pilot did not compute radial stability, (\omega_0^2), Cowling (f/p_1) modes, or a physical crust-core transition. MD-31 (sourced EOS/surface to (P=0)) and MD-32 (justified perturbation boundary condition) remain unresolved.

The result is a bounded three-pair structural pilot only. It does not authorize promotion of the minimal-repair closure as a physically preferred EOS, nor an automatic 234-table / 1,872-case replacement of the frozen projected cohort.

## 7. Evidence and handoff

- Full workflow: [42752885366](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42752885366)
- Artifact ID: `13734297622`
- Artifact SHA-256: `d8acfd1dad803e6580c97c2fcd1084d6c002ddd043abd6883fdad143d51717bc`
- Local verified artifact: `ns-cowling-r16-minimal-repair-tov-pilot-42752885366.zip`
- Derived EOS directory within artifact: `eos_repaired/`
- Case CSV: `r16_minimal_repair_tov_pilot_cases.csv`
- Pair-cutoff summary: `r16_minimal_repair_tov_pilot_summary.csv`
- Profile directory: `backgrounds/`

### Next gate

Continue with R18: localize raw pressure-secant decreases and their source/matching provenance. Keep raw CompOSE, normalized matched inputs, current q-linear projected solver inputs, and the new minimal-repair derived EOS as distinct identities. Do not scale this TOV contrast to the full 234-table cohort until R18 establishes which pressure closure is supported by the underlying source physics and transformation history.
