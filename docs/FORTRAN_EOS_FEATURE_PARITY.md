# Python EOS package to Fortran feature-parity matrix

Status values: PORT_VERIFIED means the native Fortran implementation was compiled and tested; WRAPPER_VERIFIED means the Python API invokes the native routine; PARTIAL means a bounded subset is implemented and verified while distinct legacy behaviour remains unavailable; NOT_PORTED means deliberately unavailable and no Python numerical fallback is invoked.

| Historic capability | Fortran status | Python frontend | Evidence / limitation |
|---|---|---|---|
| CompOSE eos.nb/t/yq/thermo reader | PORT_VERIFIED | WRAPPER_VERIFIED | Reads 1-D cold slices, preserves index origin and Nadd |
| Optional eos.compo reader | PORT_VERIFIED | WRAPPER_VERIFIED | Supports particle IDs 10,11,0,1 and one nuclear group |
| EOSData Q1-Q7 thermodynamic properties | PORT_VERIFIED | WRAPPER_VERIFIED | Exact source-array round trip in differential corpus |
| CompOSE writer | PORT_VERIFIED | WRAPPER_VERIFIED | 18/18 table numerical columns reproduced exactly in round trip |
| C1-C6, C10-C11 consistency verdicts | PORT_VERIFIED | WRAPPER_VERIFIED | Differential hard verdicts compared against Python reference |
| C7-C9 derivative diagnostics | PORT_VERIFIED | WRAPPER_VERIFIED | Finite-difference residuals; final parity recheck recorded separately |
| Density-gap detection | PARTIAL | Python utility | Largest adjacent density ratio is selected in the Python EOSData array helper |
| Master-function density-gap fill | PORT_VERIFIED | WRAPPER_VERIFIED | Fortran PCHIP of log(Q1), integral correction enforcing F'=Q1, causal Hermite fallback; synthetic and GPPVA differential tests |
| Cold master-function reconstruction v1 | PARTIAL (native kernel PORT_VERIFIED) | WRAPPER_VERIFIED | F(log n_B) local weighted-cubic/PCHIP or unsmoothed PCHIP; cold Q6=Q7 and CompOSE Q3/Euler gates; Il=0 Q3 convention; first-law closure; whole-interval PCHIP core and quartic outer-bridge extrema; matched outer-bridge join derivative; composition/outer-crust subset; no global-spline or microscopic-model parity |
| Historic global smoothing/candidate-selection policy and full master-function repair | NOT_PORTED | Explicitly unavailable | Legacy global UnivariateSpline smoothing-factor policy, global candidate selection and complete anomaly-repair behaviour remain unported |
| Crust-core smooth matching: maxwell/css/cubic_F/cubic_P/quintic_F/quintic_P | PARTIAL (six thermo-only methods PORT_VERIFIED) | WRAPPER_VERIFIED | All six smooth strict-thermo methods compiled and passed native-versus-legacy differential tests on GPPVA_FSU2H + GM_GM1. Composition-aware matching and literature_beta_eq remain unported. The legacy `maxwell` token is a smooth near-pressure plateau, not a density jump |
| Literature-specific RMF/Skyrme/APR solvers | PARTIAL | Not exposed as complete parity | Existing Fortran nuclear/liquid-drop/beta-equilibrium examples are not equivalent to the literature solver suite |
| Genuine one-charge Maxwell transition with local neutrality | PORT_VERIFIED (thermo-only v1) | WRAPPER_VERIFIED | Exact-head push/PR CI passed; caller confirms local neutrality; cold beta-equilibrium; one isolated transition; no Seidov/TOV interface verdict; not complete CompOSE composition output |
| Gibbs mixed phase / global electric neutrality | NOT_PORTED | Explicitly unavailable | Requires phase charge densities, independent charge chemical potential and volume fraction absent from 1-D barotropes |
| Causality profile/completion | PARTIAL | Native CSS completion wrapper; profile remains Python-only | First resolved single-phase c_s^2=1 crossing; native constant-c_s^2 tail; no phase-jump interpolation; explicit composition omission; independent profile remains unported |
| TOV star | PORT_VERIFIED in solver workspace | WRAPPER_ADDED | Native real_table_star callback |
| TOV sequence | Fortran solve per point | WRAPPER_ADDED | Python orchestrates pressure scan; each TOV integration is Fortran |
| Relativistic Cowling f/p1 | PORT_VERIFIED for solver | WRAPPER_ADDED | Uses native cell_solve; cohort results tracked separately |
| Plotting | Python presentation layer | RETAINED | Not part of thermodynamic core |

This matrix is intentionally conservative. Native master reconstruction and the one-charge Maxwell kernel are verified only within their stated v1 contracts. The historical global-spline selection, microscopic solvers, Gibbs surfaces and phase-interface stellar stability remain distinct unported or unverified capabilities; a successful CompOSE round-trip is not a complete EOS-generator parity claim.


2026-11-08: this feature remains PARTIAL; historic global-spline selection, microscopic solvers and complete BPS nuclide generation remain unimplemented.
