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