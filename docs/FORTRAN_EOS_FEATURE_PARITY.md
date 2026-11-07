# Python EOS package to Fortran feature-parity matrix

Status values: PORT_VERIFIED means the native Fortran implementation was compiled and differentially tested; WRAPPER_VERIFIED means the Python API invokes the native routine; PARTIAL means only a subset is implemented; NOT_PORTED means deliberately unavailable and no Python numerical fallback is invoked.

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
| Full master-function smoothing/candidate selection | NOT_PORTED | Explicit unavailable error | Historic global smoothing, anomaly repair and candidate-selection policy are still not ported |
| Crust-core matching: maxwell compatibility bridge | NOT_PORTED | Explicit unavailable error | Do not confuse historic smooth near-plateau method with a physical Maxwell jump |
| Crust-core matching: cubic_P/quintic_P/css/cubic_F/quintic_F | NOT_PORTED | Explicit unavailable error | Full six-method bridge selection and slope/area gates remain |
| Literature-specific RMF/Skyrme/APR solvers | PARTIAL | Not exposed as complete parity | Existing Fortran nuclear/liquid-drop/beta-equilibrium examples are not equivalent to the literature solver suite |
| Genuine Maxwell/Gibbs phase transitions | PARTIAL | Not exposed as complete parity | Existing Fortran interface toys do not implement every Python phase-transition/TOV treatment |
| Causality profile/completion | PARTIAL | Through table admission only | No full causal-completion generator parity |
| TOV star | PORT_VERIFIED in solver workspace | WRAPPER_ADDED | Native real_table_star callback |
| TOV sequence | Fortran solve per point | WRAPPER_ADDED | Python orchestrates pressure scan; each TOV integration is Fortran |
| Relativistic Cowling f/p1 | PORT_VERIFIED for solver | WRAPPER_ADDED | Uses native cell_solve; cohort results tracked separately |
| Plotting | Python presentation layer | RETAINED | Not part of thermodynamic core |

This matrix is intentionally conservative. It prevents a successful CompOSE table round-trip from being misrepresented as a complete port of the full Python EOS generation and matching research pipeline.
