# D5 topology-aware correction (2026-11-10)

This is a correction of the earlier Iphase-only exploratory Hermite replay for the 16-source C7X subset. It does not grant source admission and does not run Maxwell roots, TOV, or Cowling. The local raw-source project tree and native Fortran parser/build remain unavailable.

## Segmentation correction

The GPPVA FSU2H data sheet identifies the final Nadd=1 eos.thermo field as a discrete pasta-geometry index: 0 outer crust, 1 droplets, 2 rods, 3 slabs, 4 tubes, 5 bubbles, 6 homogeneous core. Seven GPPVA crust tables use 0-3. The earlier D5 replay split only at Iphase and therefore incorrectly included 14 within-inner-crust geometry-transition cells. RG_KDE0v1 has a discrete q_extra transition 0 to 1 between knots 1220 and 1221. Its physical meaning remains unknown; it is treated only as an opaque segmentation boundary. RG_SLY4 has q_extra=0 throughout.

Additional boundaries: DD2 275/276 and 280/281; DDME2 276/277 and 284/285; FSU2H 257/258 and 283/284; FSU2R 253/254 and 276/277; NL3wrL55 279/280 and 290/291; TM1e 259/260 and 290/291; TW 276/277 and 284/285; RG_KDE0v1 1220/1221. Including Iphase boundaries gives 45 segments and 15 excluded boundary cells.

## Exact boundary-delta counts

The 15 boundary cells were individually evaluated with the same cubic-Hermite polynomial and kernel-scaled causality statistic, then excluded from interpolation. Seven would have failed; eight would have passed. All other interval classifications are inherited unchanged from the earlier 5,550-cell D5 replay. This is an exact boundary-delta correction, not a fresh one-shot replay of every retained source cell. The updated script is syntax- and synthetic-test checked but has not been run against a local raw-source tree.

- 5,535 topology-preserving cells.
- 1,153 causal-polynomial violations.
- 1,110 cells with dn/dmu_E <= 0.
- 0 nonpositive interpolated-density cells.
- 4,382 cells pass these three selected interval conditions.

The previous source-wide minima are not carried forward as topology-aware global extrema; a full source replay is required. This is a representation-level interpolation test, not a physical instability claim about the microscopic EOS.

## DAG disposition

D5 remains exploratory and not an official admission gate. D3 Q3/Q5 failures, D4 schema mismatch, and lack of native Fortran execution remain. D6-D10 remain blocked/not run. No density drift, common-pressure equality, Maxwell root, TOV sequence or Cowling spectrum was computed.