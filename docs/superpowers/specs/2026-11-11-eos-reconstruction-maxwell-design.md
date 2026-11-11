# EOS reconstruction and genuine Maxwell validation: design

Date: 2026-11-11
Implementation line: PR #1, research/fortran-eos-gap-fill-v1
Single downstream DAG: PR #4, research/sharp-maxwell-cohort-admission-v1
Source-authority pin: research/maxwell-native-fortran-v1 at a09ec30d74b363652d2d1aefe990557c468679dd

## Intent

Resume the pre-existing EOS reconstruction and matching programme. Validate native reconstruct_master and all six smooth crust-core methods against cold first-law identities, then validate the separate one-charge, locally neutral construct_maxwell_transition kernel. The immediate source-backed defect is the Il=0 Q3 convention mismatch documented in the Maxwell audit. The existing C4 JSON plan remains the only machine-readable workflow DAG; this document is a design companion, not a second DAG.

## Goals

1. Use one master function F(ln n_B) for single-phase cold reconstruction; derive P, epsilon, mu_E and c_s^2 consistently.
2. Validate PCHIP/Hermite interval stability and causality over complete polynomial intervals, not only at source knots.
3. Preserve historical matching method tokens and semantics, including method=maxwell as a smooth near-pressure plateau.
4. Validate true Maxwell coexistence in Euler mu_E: common pressure and chemical potential, density/energy jump, and Delta epsilon=mu_t Delta n_B.
5. Preserve CompOSE Q3 convention by Il in both validation and output.
6. Keep raw EOS source rows, original order, pair manifests and legacy paths immutable.
7. Keep composition/provider/C6B blocker and 65-pair root admission gates intact.

## Non-goals and anti-goals

- No Gibbs/global-charge-neutrality construction from one-dimensional EOS tables.
- No reconstruction of missing particle composition or local neutrality inference.
- No provider-specific interpretation of Iphase, Nadd, q_extra or eos.compo without source evidence.
- No Seidov stability verdict, TOV phase-chain implementation or Cowling spectrum promotion in this increment.
- No production root sweep across 65 crust-core pairs while C6/C6B/C7 admission remains blocked.
- No alteration of the historical smooth maxwell token into a physical density jump.

## Approaches considered

1. One unified Maxwell path for all matching: rejected. It conflates the historical smooth plateau with physical phase coexistence and changes established API semantics.
2. Validate reconstruction, smooth matching and true Maxwell separately while sharing thermodynamic invariant helpers/receipts: selected. It follows the existing modules and preserves the physics distinction.
3. Python-only reference implementation replacing native Fortran: rejected. Native Fortran is the implementation authority; Python remains the wrapper and independent end-to-end assertion layer.

## Architecture

- compose_eos_reconstruct.f90: F(log n_B)-based single-phase reconstruction and existing Il=0 Q3 handling.
- compose_eos_matching.f90: historical smooth crust-core methods; maxwell remains smooth_near_pressure_plateau.
- compose_eos_phase_transition.f90: true one-charge Maxwell construction in P(mu_E), one isolated transition only.
- compose_eos_verify.f90: native C1-C11 consistency checks, including C2 Il-aware Q3 relation and phase-boundary derivative masks.
- Python tests: analytic synthetic EOS and C2/serialization integration checks.
- C4 plan JSON on the PR #4 branch: sole machine-readable DAG for future source intake, pair admission and downstream TOV/Cowling.

## Frozen thermodynamic conventions

Use repository units: n_B in fm^-3; P and epsilon in MeV fm^-3; m_n and chemical potentials in MeV. c and hbar remain explicit in dimensional code and conversions.

Single-phase master representation:

    F(x)=f/n_B, x=ln n_B
    Q1=dF/dx
    P=n_B*Q1
    epsilon=n_B*F
    mu_E=F+Q1
    cs2/c2=(Q1+F'')/(F+Q1)

True one-charge, local-neutrality Maxwell representation:

    mu_E=(epsilon+P)/n_B
    n_B=dP/dmu_E
    epsilon=mu_E*n_B-P

At coexistence:

    P_L(mu_t)=P_H(mu_t)=P_t
    mu_L=mu_H=mu_t
    n_H>n_L
    Delta epsilon=mu_t*Delta n_B

Supported cold CompOSE Q3 convention:

- Il=0: mn*(1+Q3)=mu_E-Yq*mn*Q4.
- Il=1 and Q5 approximately zero: mn*(1+Q3)=mu_E.

The Maxwell root search uses Euler mu_E; output Q3 remains in the CompOSE convention. True Maxwell v1 kernel Q3/Euler tolerance 5e-4 is distinct from strict reconstruction/source-admission tolerance 2e-6.

## Matching semantics

| Method | Contract | Exact jump? |
| --- | --- | --- |
| historical maxwell | smooth near-pressure plateau, first-law integral closure | no |
| css | smooth constant-sound-speed bridge | no |
| cubic_F, quintic_F | smooth free-energy interpolation | no |
| cubic_P, quintic_P | smooth pressure/chemical-potential interpolation | no |
| construct_maxwell_transition | two phase coexistence rows at equal P and Euler mu | yes |

## Existing design decisions

- Correct native Il=0 Q3 validation and output serialization with a focused failing regression first.
- Preserve Il=1 analytic Maxwell fixtures and all historical smooth matching tests.
- Keep true Maxwell v1 input limitations explicit: T=0, beta equilibrium, Q2=0, Q5 approximately zero, Q6=Q7, one smooth segment per phase, Nadd=0 and composition-free thermo-only inputs.
- Keep local-neutrality confirmation caller-supplied; the code cannot infer it from an EOS barotrope.
- Any widening to composition, Nadd, internal Iphase boundaries or nonzero Q5 requires a separately reviewed physics/metadata contract.

## Success conditions

- Native build with runtime checks passes.
- Il=0 regression fails before the fix and passes after it.
- Full frontend/pipeline suite passes with 68 tests at the current checkpoint.
- Native quartic polynomial-extrema test passes.
- SHA256SUMS_v4.txt validates every tracked entry after all edits.
- PR #1 and PR #4 stay unmerged research branches.
