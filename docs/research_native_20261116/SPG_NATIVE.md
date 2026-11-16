# Native SPG calculation and endpoint verification

All physical evaluation, constrained cell minimization, neutron-drip search and upper endpoint search execute in Fortran. Python only transports source parameters and records validation evidence. The surface coefficients remain explicit candidates, not recovered original author calibrations.

## Computed endpoint contract

The source prescription compares clustered and homogeneous energies at the same baryon density. The native transition routine scans stationary clustered states over 0.002–0.12 fm^-3, brackets a sign change of epsilon_cluster−epsilon_uniform and bisects it. It does not use CompOSE transition rows. The separately exposed common-tangent routine is an alternative thermodynamic diagnostic; it is not substituted for the source condition.

| Model | Native drip (fm^-3) | Fixed-density crossing (fm^-3) | Rows | Maximum projected KKT (MeV) |
|---|---:|---:|---:|---:|
| SPG_M1 | 0.0001940206412 | 0.0746578940153 | 129 | 1.995e-08 |
| SPG_M2 | 0.000205001488786 | 0.0879137359858 | 129 | 1.837e-08 |
| SPG_M3 | 0.00019386210629 | 0.0791172665358 | 129 | 1.991e-08 |
| SPG_M4 | 0.00020210860012 | 0.0973067723513 | 129 | 1.966e-08 |
| SPG_M5 | 0.000202958941605 | 0.0718003231287 | 129 | 1.846e-08 |

## Validation and scope

- The native campaign generated 985 stationary rows: five AME2408 candidate tables with 129 rows and twenty alternative-surface tables with 17 rows. All 25 native drip-to-crossing runs passed baryon and charge conservation, positive pressure, and projected KKT < 1e-5 MeV.
- Dedicated tests compare all five bulk models and analytic cell gradients with the retained research implementation; independent finite differences of reoptimized energy verify mu_B=d epsilon/d n_B and P=n_B mu_B−epsilon at three interior densities per model.
- The finite-dimensional search checks local curvature and gas-proton/muon insertion margins. These checks do not prove a global minimum or recover an unavailable original surface calibration.
- Two M5 geometric-surface variants exposed a distinction between general minimum selection and source endpoint branch comparison. Collapsed attempts on the artificial ni≈n_B,u→1 boundary or vanishing cluster fraction t→0 are analytic uniform limits. The endpoint branch routine retains these attempts but compares a separately stationary finite cluster against independently solved uniform matter. Other lower-energy unconverged interior trials remain a failure. The general cell-selection policy is unchanged.
- The failing pre-repair campaign and an explicit M5 attempt dump are preserved. Repetitive compiler runtime stderr is losslessly gzipped; STDERR_ARCHIVE.json records original byte count/hash and verifies decompressed byte identity.
- Source provenance and executable/source hashes are in data/research_native_20261116/spg/campaign/RECEIPT.json. Native readback through the public CLI is a separate parent integration check.
