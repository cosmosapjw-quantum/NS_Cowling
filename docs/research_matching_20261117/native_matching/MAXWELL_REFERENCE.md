# Maxwell reference-normalization repair

This is a bounded integration repair discovered in the parent's actual native-generated EOS campaign. Twenty SPG/non-APR Maxwell pairs were rejected solely because their CompOSE header reference masses differed. Five additional APR cases were independently blocked by unsupported multiphase retention. Those input-support failures are not addressed by this normalization patch.

## Physical contract and root cause

At zero temperature, the physical source quantities are P=n Q1, epsilon=n M(1+Q7), mu_B=M(1+Q3), and mu_Q=M Q4, where M is the input's neutron rest-energy encoding reference in MeV. The Euler potential for Il=0 is mu_E=mu_B+Yq mu_Q; for Il=1 it is mu_B. Each source already passed native validation in these physical variables.

The former pairwise header equality gate conflated reference encoding with the physical model's rest-energy convention. Removing that gate alone would have corrupted output: copy_source_row copied raw Q7 and Q4 from the high branch into a table normalized by the low branch's M. The transition-row Q4 interpolation had the same normalization assumption.

The corrected native kernel decodes each source with its own reference and re-encodes Q4, Q6 and Q7 using the output reference. Q3 follows the already-validated physical Euler chemical potential, with the Il=0 charge-potential term retained. P, epsilon, mu_B and mu_Q are preserved. No constant is fitted or added to either source energy.

A physical transformation epsilon -> epsilon+n C and mu_E -> mu_E+C remains a different EOS. It may alter or eliminate coexistence. Header-reference normalization does not remove C, nor recover unknown source microphysics.

The parent-owned native JSON output records input/output [mn,mp] references in MeV, physical_energy_shift_mev=0, and mass_normalization_semantics=physical_energy_preserving_Q_reencoding.

## Actual red/green evidence

The analytic native fixture uses x=mu_E-990 MeV,
P_low=8+0.12 x+0.001 x^2 MeV fm^-3, and P_high=P_low+0.01 x, with coefficient units understood. It has an upward density jump at mu_E=990 MeV, P=8 MeV fm^-3, n_low=.12 fm^-3, n_high=.13 fm^-3. Source chemical-potential knots are [950,960,980,1000,1030,1040] MeV. Q4 encodes mu_Q=-40 MeV and Yq=.1, so both Il=0 and Il=1 exercise nonzero charge potential.

Before correction, all six equivalent-reference variants (neutron-only, proton-only, both; Il=0/1) were rejected by the header equality gate. After correction, all six pass and the maximum physical output-energy difference is 2.8421709430404007e-14 MeV fm^-3. Pressure, both chemical potentials, coexistence density and transition chemical potential also agree within the recorded numerical bounds. A physical +100 MeV per baryon high-phase shift still produces NO_COMMON_MU_SUPPORT for both Il conventions.

The public Python tests include independently changed and decimal-rounded reference metadata and validate the complete source/output invariants and the native provenance report. The old metadata-equality assertions in test_maxwell_mass_gate.py and the first four cases of test_mass_normalization_contract.py were replaced because they enforced the incorrect gate. Physical-offset, source-immutability and header-only rewrite controls remain.

Fresh verification:

- Native reference driver, bounds checks and zero/invalid/overflow traps: exit 0, FAILURES=0.
- Seven affected public/native Maxwell test files: 84 passed in 3.22 s, exit 0.
- The first resumed Python run had two test-message regex mismatches; both calculations correctly rejected the shifted input. That log is retained, and the corrected exact native-message assertion passes. No numerical tolerance changed.
- The parent's broader actual-pair rerun and final full regression are separate integration evidence.

## Reproduction

Use bash fortran_eos/build.sh --ensure build/matching_native/maxwell_api, then set COMPOSE_EOS_FORTRAN_EXE to that directory's compose_eos_fortran and PYTHONPATH=python_frontend:scripts. Run the seven test files listed in data/research_matching_20261117/native_matching/maxwell_reference/FINAL_RECEIPT.json. The standalone Fortran driver links compose_eos_types.o and compose_eos_phase_transition.o from the same build.

All production physics remains in Fortran. Python constructs analytic test fixtures, transports tables and compares results. The actual baseline source copy, original red/green logs, rebuilt source/executable identities and final logs are preserved under data/research_matching_20261117/native_matching/maxwell_reference/.
