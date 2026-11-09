# Native Fortran sharp-Maxwell TOV/Cowling path (v1)

## Scope

This branch adds a native Fortran path for a single sharp first-order transition. It is independent of the earlier 234-row smooth-matching cohort and does not reinterpret the historical maxwell bridge as a density jump.

The implementation accepts two independently admitted normalized EOS tables: a high-density provider for the stellar core, queried only at/above the transition, and a low-density provider for the envelope, queried only at/below the transition.

The production-facing case-driver names are maxwell_table_star and maxwell_table_cowling. The second also evaluates one Cowling shooting residual at a requested frequency and phase-conversion prescription.

## Thermodynamic interface contract

At the transition pressure P_t, the two providers must satisfy

\[ P_-=P_+=P_t, \qquad \mu_{B,-}=\mu_{B,+}=\mu_t, \]

with positive density and energy jumps

\[ \Delta n_B=n_+-n_->0, \qquad \Delta\epsilon=\epsilon_+-\epsilon_->0, \]

and cold first-law closure

\[ \Delta\epsilon=\mu_t\Delta n_B. \]

The native TOV entry point checks pressure equality, chemical-potential equality, positive phase ordering and the first-law identity before returning a valid star. It does not build a coexistence construction from arbitrary EOS tables: the input high/low tables must already represent the two coexisting phases at the supplied P_t.

## Numerical construction

The high-density core and low-density envelope are integrated on separate logarithmic pressure coordinates. The core is stopped at the exact endpoint P_t; the envelope is restarted at P_t using the same r_t, m_t and nu_t. Each phase is retained as its own monotone background table, so no interpolation or derivative is evaluated across the discontinuity.

The TOV equations use SI units and the finite-pressure-cutoff surface convention already used by the native smooth solver. The lapse is normalized only after the envelope integration reaches the declared pressure cutoff.

The Cowling shooting path integrates (W,V) independently through the core and envelope. At the interface it applies one of two maps:

- slow/material conversion: [xi^r]=0, [Delta P]=0;
- rapid conversion: [Delta P]=0, [xi^r-Delta P/P'_0]=0.

The code keeps the inner and outer one-sided backgrounds separate and uses the native Fortran interface-transfer routines. The returned frequency residual is a single-frequency shooting diagnostic, not a complete f/p1 spectrum or an eigenmode-family identification.

## Independent validation plan

tests_fortran/maxwell_checks.f90 uses an analytic cold quadratic barotrope with a constructed, nonzero Maxwell density jump. It checks:

1. exact arrival at P_t from the core and exact restart at P_t in the envelope;
2. continuity of radius and enclosed mass at the interface;
3. chemical-potential continuity and the cold first-law jump identity;
4. the Seidov rapid-onset diagnostic for a deliberately small stable-onset jump;
5. native slow- and rapid-conversion Cowling shots and conservation of Delta P;
6. the small-jump limit against the pre-existing single-EOS Fortran TOV implementation.

The analytic test is deliberately separate from any student EOS cohort. It supplies a reproducible numerical validation fixture, not a phenomenological claim about realistic hybrid stars.

## Build and test

From oscillation_source/solutions:

    fpm test --compiler gfortran --flag '-std=f2018 -O0 -g -fcheck=all -ffree-line-length-none -fimplicit-none'
    fpm build --compiler gfortran --flag '-std=f2018 -O0 -g -fcheck=all -ffree-line-length-none -fimplicit-none'

The GitHub Actions workflow .github/workflows/maxwell-native-fortran.yml runs these commands on the exact pushed branch head and uploads the complete test log.

## Table-driven use

Build the case_driver executable, then pass maxwell_table_star or maxwell_table_cowling with a NAMELIST file. The parameter names are:

    &request
      request_token = 'maxwell-native-validation'
    /
    &parameters
      core_table_path = 'path/to/high_density_phase.csv'
      envelope_table_path = 'path/to/low_density_phase.csv'
      pc_pa = 1.4e35
      pt_pa = 1.0e35
      p_cut_pa = 1.0e30
      r0_m = 1.0
      r_limit_m = 1.0e6
      frequency_hz = 1800.0
      conversion = 'slow'
    /

The numerical values above are placeholders; the two tables and transition pressure must be mutually consistent and use the normalized table wire format accepted by ns_table_eos. frequency_hz=0 skips Cowling even for maxwell_table_cowling.

## Explicit limits

- One sharp Maxwell interface only; no multiple transitions or Gibbs mixed phase.
- The core/envelope tables are supplied separately; no automatic phase coexistence root search is performed by this module.
- Cowling support is ell=2, barotropic and single-frequency shooting. Root scanning, mode continuation and robust identification of interface g-modes are not yet implemented.
- The Seidov criterion is an onset diagnostic only; it is not a substitute for a complete radial-stability sequence.
- No full-GR metric perturbations, rotation, magnetic fields, elasticity, superfluidity, composition stratification or gravitational-wave damping.
- The surface remains a finite positive-pressure cutoff.
- This branch does not contain the previously prepared 33-path v5 source overlay. The v5 sensitivity report and this physical phase-transition solver remain separate research tracks.

A passing CI run establishes compilation and the checks in this fixture only. It does not establish agreement with realistic CompOSE hybrid EOS, a published mode benchmark, or a production M(R)/spectrum survey.

## Exact-head verification (2026-11-09)

Verified source head: `98d2f9ce6033015109286e5713639c6626a6f7fb` on `research/maxwell-native-fortran-v1`, tracked by draft [PR #3](https://github.com/cosmosapjw-quantum/NS_Cowling/pull/3). Both exact-head GitHub Actions events passed:

- PR run [42581660400](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42581660400): success; artifact ID `13656748986`, SHA-256 `083f3eabf9215ac69cde06677f06c0fc58746bb35b24e5559c3a6184b5e7c5f6`.
- Push run [42581655510](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42581655510): success; artifact ID `13656434783`, SHA-256 `b6274813105f8e41f9eec269ffbe747ee3c5b594aea3a7bf36c3427a0cfb665e`.
- The native fixture reports 12 existing cell checks and 18 Maxwell checks passed. Its small-jump TOV comparison gives `dM/M=9.1343e-8` and `dR/R=2.8981e-8`.
- The exact-head workflow compiled with GNU Fortran 15.2.0, ran FPM 0.13.0 with bounds/runtime checks, built the production case driver, and uploaded source-identity evidence.

These results establish the bounded analytic fixture and build path only. They do not establish a real normalized two-table EOS execution, a published hybrid-star benchmark, a complete mode spectrum, a stable-branch survey or observational inference. See `docs/ADVERSARIAL_AUDIT_MAXWELL_NATIVE_FORTRAN_20261109.md` and `docs/MAXWELL_NATIVE_FORTRAN_DAG_20261109.md` for the audit and remaining gates.
## Synthetic normalized-table and production CLI verification (2026-11-09)

The later integration gate closes the previous normalized-table runtime gap on a generated analytic fixture. Final source head: `d3649906292b0ec8f1caf98961f12436fa00871b`.

- `tests_fortran/maxwell_table_checks.f90` writes two independent 19-column T14 normalized CSV tables, with `g=epsilon-P` linear in `ln(P)` so the current monotone PCHIP representation is exact and the first-law density quadrature is self-consistent.
- Both tables load through `load_normalized_table`, pass `init_table_eos`, and return one-sided states at the same `P_t`. Chemical-potential coexistence, positive density/energy jumps and cold first-law closure are checked.
- Table-backed native phase-segmented TOV and both slow/rapid Cowling shots pass. The table fixture reports 74 checks passed.
- The GitHub Actions workflow then invokes the actual production `case_driver` with `maxwell_table_cowling` and the generated NAMELIST. The case returns `status=ok`, `result_valid=true`, and `cowling_shot_valid=true`.
- For this synthetic fixture only: `M=9.6356877e-4` nominal solar masses, `R=642.793 m`, and `r_t=332.172 m`; the Seidov ratio is `0.21276596`, below its `0.50159574` critical value.
- At the 1800 Hz Cowling trial, the surface residual is `0.9935236`, so this is not a mode root. The interface pressure perturbations agree across the interface to approximately `3.2e-16` relative.

Exact-head GitHub Actions:

- Push run [42582837414](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42582837414): success; artifact ID `13657672304`, SHA-256 `ec8333d0afc1fb65b57277c669333d403fbbbd6e18a19fb701aa369397bc2515`.
- Pull-request run [42582845328](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42582845328): success; artifact ID `13658265995`, SHA-256 `5d5a1683e4ba708f246c275e225b462bfc6968f201cf02d4956c9ceecf2cc278`.
- Both events ran 12 existing cell checks, 18 analytic Maxwell checks and 74 normalized-table checks, built the production driver, ran the table-driven CLI smoke, and recorded exact source identities.

Important boundary: these CSV tables are generated analytic fixtures, not source-backed CompOSE or other realistic hybrid EOS data. This verifies the normalized-table/provider/CLI integration contract, not realistic hybrid-star phenomenology. The old 33-path v5 overlay remains unapplied.

## Production generic dispatch extension (2026-11-09)

The low-level kernels remain separately callable, but production callers now have one typed generic API in `src/production/ns_production_solvers.f90`:

```fortran
call solve_tov_star(eos, pc, options, smooth_profile)
call solve_tov_star(eos_high, eos_low, pc, pt, options, maxwell_star)

call solve_cowling_spectrum(smooth_background, omega_grid, mode_options, smooth_spectrum)
call solve_cowling_spectrum(smooth_profile, omega_grid, mode_options, smooth_spectrum)
call solve_cowling_spectrum(maxwell_star, omega_grid, mode_options, 'slow', maxwell_spectrum)
call solve_cowling_spectrum(maxwell_star, omega_grid, mode_options, 'rapid', maxwell_spectrum)
```

The smooth TOV overload also preserves the optional `eos_events(:)` smooth-stop argument. The Maxwell overload still requires two callbacks whose already-admitted tables satisfy coexistence at the supplied `P_t`; it does not construct coexistence or extrapolate either phase table.

The Maxwell spectrum implementation lives in `src/modes/ns_maxwell_mode_spectrum.f90`. It scans the supplied positive, strictly increasing angular-frequency grid, retains failed trials as barriers, finds candidates using the existing bracket builder, refines sign/near brackets with bisection, and reruns the shot at the candidate root. The conversion prescription is fixed for one scan (`slow` or `rapid`). Root-family identity is always reported as `unresolved`; the implementation does not yet assign f/p1/interface labels. A finite scan is only a bounded search over the supplied grid and is not a completeness theorem.

Production case adapters now route existing smooth TOV, normalized-table TOV, and smooth Cowling spectrum calls through the generic API. The Maxwell case adapter routes TOV through the same generic API and adds `maxwell_table_spectrum`. Current exact-head CI uses a generated analytic two-table fixture and checks the scan/no-root path; positive-root bisection, realistic CompOSE data, root-family identity and physical spectrum convergence remain open.


## Exact-head generic-dispatch and spectrum CLI update (2026-11-09)

The updated production smoke now runs both `maxwell_table_cowling` and `maxwell_table_spectrum` through `case_driver`. The finite frequency grid is 1799, 1800 and 1801 Hz; all 3 shots return status `ok`, with residuals (0.989787), (0.993524) and (0.996295). There are no sign brackets, sample candidates or roots. Thus the production scanner's successful no-root path is verified, but no physical eigenfrequency has been found and the bisection-positive-root branch is not exercised by this fixture.

The same exact-head CI also checks 12 cell, 22 analytic Maxwell and 74 normalized-table assertions, and both push and pull-request events pass:
- Push [42585358415](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42585358415), artifact 13658722733, SHA-256 `26313f4315af74718f5fa9d1687658ef603c405807b2f7d44d27057fa38242a1`.
- Pull request [42585366571](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42585366571), artifact 13659056245, SHA-256 `12afce4e910c53b9a0dda180826bea2a98a5986562df0866bb7f9a26da219412`.

CI-validated workflow head: `22b0faa3c69150f3ae4fffcd2a6cae5d8231bb6c`. Code-bearing implementation head: `f8678cb13a18f63bad9f8db65c4cd3e66e491ad3`. The normalized tables remain generated analytic fixtures; realistic CompOSE data, positive-root convergence and mode identity remain unverified.


## Exact-head positive-root validation for both conversion prescriptions (2026-11-09)

The normalized-table fixture now writes separate slow- and rapid-conversion NAMELIST inputs. Both use a frequency grid corresponding to 100-20,000 Hz with 41 uniformly spaced samples, and both run through the actual production `case_driver maxwell_table_spectrum` route. The bisection routine refines each sign bracket and the scanner reruns the Cowling shot at each resulting root.

Slow conversion:
- Root 1: (1764.8989380477\,\mathrm{Hz}), final surface residual (-6.9575\times10^{-7}).
- Root 2: (4937.6449322607\,\mathrm{Hz}), final surface residual (9.3334\times10^{-8}).

Rapid conversion:
- Root 1: (4937.0162338857\,\mathrm{Hz}), final surface residual (8.5331\times10^{-9}).
- Root 2: (19493.6326018255\,\mathrm{Hz}), final surface residual (3.8951\times10^{-7}).

All four roots report `root_valid=true`, status `ok`, and finite final residuals below the CLI fixture's (10^{-5}) acceptance threshold. Their family labels are intentionally `unresolved`. These are numerical roots for a generated analytic two-phase regression EOS, not realistic neutron-star frequencies and not f/p1/interface-mode identifications.

The strict workflow requires `spectrum_roots_found=2` and `spectrum_root_valid_00001/00002=true` independently for slow and rapid conversions. Exact code/test head: `9e92655a8277922f312f4c31df8e331e8c12df1b`. Exact workflow head: `ec4527f0d518f6ba19f57eb42d8a8eb00b98dfca`.

Both exact-head CI events passed:
- Push [42586438396](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42586438396), artifact 13658709148, SHA-256 `02d1673cd82045aaacbae8125c74c7ce7b8896b7abe6def70e81df07331c3d04`.
- Pull request [42586446195](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42586446195), artifact 13659137810, SHA-256 `478fb960cf8a993bb92bf2bbc544f17160b854e0f0f5ac5a28cc0b14a9ccaeaf`.

This verifies the positive-root branch, not a root-convergence study: only one ODE/root-tolerance configuration is used. Frequency-grid refinement, ODE tolerance refinement, independent Python/SciPy comparison and mode identity/continuation remain open.
