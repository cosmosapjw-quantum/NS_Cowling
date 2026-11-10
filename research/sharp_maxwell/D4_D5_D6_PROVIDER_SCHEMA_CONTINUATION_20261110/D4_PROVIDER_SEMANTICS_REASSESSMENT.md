# D4 provider semantics reassessment: CompOSE crust composition fields

**Date:** 2026-11-10  **Status:** D4_PARTIAL_BLOCKED_PROVIDER_SCHEMA_MISMATCH  **Scope:** 13 EOS / 4,975 eos.compo rows.

## Canonical schema

CompOSE v3.01 defines particle fractions as Y_i=n_i/n_b and imposes sum_i B_i Y_i=1. A nuclear-group quadruple (I,A_av,Z_av,Y_I) stores Y_I=sum_{j in M_I}Y_j=(sum_j A_jY_j)/A_av. Its baryon and charge contributions are A_av Y_I and Z_av Y_I. The separate mass fraction is X_I=A_av Y_I. Model-specific freedom covers the group-index mapping and Iphase interpretation, not the standard meaning of the fourth field.

## Raw-data finding

All 13 files parse with zero malformed rows and zero inb-index sequence mismatches. In the 11 GPPVA/SPG crust sources (3,538 rows), there are 2,387 Iphase=0 outer-crust rows and 1,151 Iphase=7 inner-crust rows. Outer crust has Y_n=Y_p=0, Y_quad=1, and Y_e approximately Z_av/A_av: a mass-fraction-like diagnostic closes to printed precision, but canonical number-fraction closure fails by about A_av-1. In every one of the 1,151 inner-crust rows, Y_quad equals Y_e exactly at printed precision and Y_n=Y_p=0. Neither canonical number-fraction closure nor the simple mass-fraction diagnostic closes.

RG(SLy4) and RG(KDE0v1) exhibit a different, mass-fraction-like closure signature in most/all inhomogeneous rows. RG(SLy4) has one terminal heavy-nucleus outlier at row 143. These observations do not establish provider semantics and do not authorize a production remapping.

## Source documentation

The GPPVA FSU2H data sheet documents its extra eos.thermo quantity as a pasta-geometry index and lists phase 0/7/6 as outer/inner/homogeneous. The SPG M2 data sheet likewise lists four particle pairs, one heavy-nucleus quadruple, and phase 0/7. RG(SLy4) documents phase 4 as heavy nucleus present and phase 3 as homogeneous matter. None of these inspected sheets redefines the nuclear-group fourth field as mass fraction, electron fraction, volume fraction, or cell occupancy.

## Interpretation and decision

The evidence supports a schema/export mismatch or missing provider documentation, but it does not identify which software stage introduced the issue. Do not assign the inner-crust Y_quad an inferred physical meaning. Do not remap it to a mass fraction or electron fraction. Keep D4 blocked and quarantine composition from production adapters pending provider/CompOSE clarification.

## Clarification request

Ask each provider to specify (1) exact quantity and denominator in the fourth quadruple slot by phase; (2) whether free n/p pair values are global baryon fractions, local cell fractions, or otherwise normalized; (3) baryon and charge closure equations used; (4) whether inner-crust Y_quad=Y_e is intentional or a faulty export; (5) meaning of RG_SLy4 terminal row 143; (6) per-model GPPVA Nadd geometry and RG_KDE0v1 q_extra; (7) whether corrected eos.compo files/data sheets will be issued. No inquiry has been sent.