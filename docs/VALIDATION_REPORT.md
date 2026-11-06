# v3 direct-EOS integration validation

## Evidence status

**Python tests:** 8 passed in the full `pipeline/tests` suite; the focused direct-provider tests are 4/4 passing.

**Fortran build:** all 62 support/solution modules compiled with gfortran in Fortran 2018 mode, runtime checks enabled, and both `case_driver` and `cell_solve` linked. The repository's declared FPM harness was not run because FPM 0.13.0 is absent on the target machine. Do not interpret the manual compiler build as an FPM test pass.

**Fortran production EOS cohort:** 234/234 original source tables admitted; 109,776 cells total; 0 rejected cells; 0 raw first-law gate failures. Maximum source knot closure residual (6.1990\times10^{-16}). Maximum relative knot difference between reconstructed and source (n_B), (9.8825\times10^{-4}). Machine-readable reports are in `data/eos_direct_provider_v3/fortran_direct_eos_cohort.csv` and `.json`.

**Negative fixture:** `oscillation_source/data/normalized/compose_134_si.csv) is rejected because its raw first-law knot residual is (1.1327081\times10^{-4}), above the (2\times10^{-12}) admission limit. The continuous provider does not conceal inconsistent input by silently reconstructing its missing thermodynamics.

**Representative mode calculation:** `GPPVA_FSU2H__BBB_BHF2__quintic_P), direct provider, target mass (2.644606257493485\times10^{30}\,\mathrm{kg}). The actual mass is (2.644606257493893\times10^{30}\,\mathrm{kg}) (relative residual (1.54\times10^{-13})); radius (11.2249075\,\mathrm{km}); (f=2.517378706\,\mathrm{kHz}); (p_1=7.205657666\,\mathrm{kHz}). Matching residuals are (4.16\times10^{-9}) and (5.61\times10^{-12}). Node identity remains (f:(0,0)), (p_1:(1,1)) in the inner segment and ((0,0)) in the outer segment.

## Scope restrictions

The existing v2 234-row mode cohort remains the q-linear-adapter cohort. Only the representative EOS above has been recalculated with v3; the other 233 v3 spectra are not claimed complete. No full-GR metric perturbation, gravitational-wave damping, frozen-composition response, superfluid modes, sharp-interface conversion, or zero-pressure extrapolation is included. The surface remains at a finite positive pressure cutoff.
