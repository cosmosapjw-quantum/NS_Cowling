# Independent bounded review: primary-source crust generators

Date: 2026-11-23
Reviewed branch/worktree: research/eos-primary-crust-generation-20261123,
`/home/cosmosapjw/ns_cowling_wt/cowling_c20_native_tov_pchip_20261123`.
Baseline HEAD: `7ad3040d93215db67b23041987e17ea19b6c6715`.

This is an independent read-only implementation/source review, not a project completion verdict.
Production code/tests and the original Dropbox checkout were not modified.
Review outputs are confined to the authorized review directories.

## Finding status

No actionable implementation or source-transcription defect was identified in
the reviewed declared model scopes and exercised input contracts. This statement
does not certify full original BPS mass-search, BBP CLDM, or NV Hartree-Fock
reproduction; each is explicitly outside the implemented identities.

One material scientific admission issue must remain separate from native finite-table
PASS: the raw primary BBP table cannot satisfy stable cold first-law interpolation
through its first two printed knots within the explicitly assumed independent
half-last-digit rounding boxes. The implementation already records this exclusion.
All 137 raw-BBP kernel survivors retain those two knots, so none is a physically
admitted continuous cold barotrope under that source contract. Root's planned
separate source-admission adjudication had not yet been inspected when this note
was written. The raw kernel receipts should remain intact.

## Sources and formulas actually inspected

1. BPS1971 original PDF, DOI 10.1086/151216:
   - Text of printed pp301-304 and scan p304 (Table 2 and Eqs14-18).
   - Scan of Table5 continued p309, all 35 BBP rows.
   - The 12 Table2 isotopes and rounded binding energies agree with the CSV.
   - Electron energy/pressure, bcc lattice energy, fixed-P Gibbs chemical potential,
     and neutron-drip criterion match the stated original-equation scope.
   - Source pre-neutron-drip mass density denotes total energy density/c^2.
2. Negele1974 original PDF, DOI 10.1017/S0074180900099885:
   - Scans of p20 TableII, p23 TableIII, and p24 author fit/equation.
   - All eight fit coefficients agree.
   - Author formula is E_T=m_n+c0+exp(sum from i=1 to7 of c_i*x^(i-1)),
     x=ln(n_b[cm^-3]*10^-35), equivalent to ln(n_b[fm^-3]/10^-4).
   - The TableII row at 5.77e36 cm^-3 reads 1.966 MeV in the 1974 source.
   - The author explicitly cautions that second derivatives amplify energy-fit noise.
     Machine roundoff control does not remove that source approximation.

## Implementation inspection

### BPS

`fortran_eos/src/compose_eos_bps_primary.f90`:
- Lines106-139: dimensions and rest-energy bookkeeping are consistent; nonrelativistic
  electron series coefficients match the analytic integrals and avoid cancellation.
- Lines154-239: fixed-P inversion, global minimum among the supplied 12 nuclei,
  equality of Gibbs chemical potentials and positive density/energy jumps are checked.
- Lines242-359: exact transition endpoints plus one interior per phase are inserted;
  lower support and drip are not extrapolated; phase IDs101-112 and composition masks
  are transported to the common carrier.
- Lines304-314: first-law finite difference removes analytically cancelling rest energy.
- Lines384-398: restricted candidate set, modern constants, rounded historical inputs
  and non-reproduction of the full historical mass catalogue remain explicit.

No sign, factor-of-c, lattice normalization, or first-law defect was found.
The independent 43-row published-table comparison is useful evidence distinct
from internal thermodynamic identities.

### BBP

`fortran_eos/src/compose_eos_bbp_primary.f90`:
- Lines9-16 and152-154: CGS-to-MeV/fm^3 and cm^-3-to-fm^-3 conversions are correct.
- Lines75-90: raw projection retains the 26 declared inner-domain knots; header
  939MeV is an encoding convention and is not added a second time to epsilon.
- Lines236-273: the piecewise-polytropic pressure integral implements
  d(epsilon/n)/dn=P/n^2, using one original energy anchor and stable exprel.
- Lines288-300: rounding-box no-go bounds are conservative with simultaneous
  allowed changes in n, P, and epsilon. For the first interval:
  Delta(epsilon/n) in [0.6088004576851063,1.6494109687876062]MeV;
  allowed pressure integral in [0.2684345407837055,0.2854525535769782]MeV.
  The separating gap is approximately0.32335MeV and is not roundoff.
- Lines307-313: continuous segment sound-speed bounds correctly use both endpoint
  limits and the monotonicity of r=P/(epsilon+P).
- Lines341-387: raw source table versus derived pressure-primary representation
  are distinct identities; the derived model explicitly changes source epsilon.

The no-go is about the finite printed source plus a stated rounding assumption;
it is not a theorem that the original BBP microphysical functional is impossible.

### NV

`fortran_eos/src/compose_eos_nv_primary.f90`:
- Lines99-182 validate exact published numerical source identity and stated
  present-day rest-mass convention.
- Lines185-217 implement coupled Horner evaluation of Q,Q',Q'' and
  P=n*exp(Q)*Q', mu=m_n+c0+exp(Q)*(1+Q'),
  c_s^2/c^2=[exp(Q)*(Q'+Q'^2+Q'')]/mu.
- Lines250-259 omit unavailable composition rather than inventing it.
- Lines275-334 retain the 1974-versus-historical-1973 energy distinction, compare
  with rounded HF eigenvalues as diagnostics, and explicitly report the fit's
  near-drip chemical potential below the free-neutron threshold.

No polynomial-index, logarithm-base, derivative, unit, or fit-source transcription
defect was found. The near-drip threshold issue is a known source-fit limitation,
not a silent repair or a newly introduced Fortran formula error.

## Frontend and orchestration

Inspected:
- `fortran_eos/app/compose_eos_generate.f90`
- `fortran_eos/build.sh`
- `python_frontend/compose_eos/generation_dispatch.py`
- `python_frontend/compose_eos/native_generation.py`
- `python_frontend/compose_eos/primary_model_contracts.py`
- `scripts/prepare_eos_primary_crust_sources_20261123.py`

The production EOS calculations are in Fortran. Python schedules, marshals,
records identities, verifies readback and builds immutable matching views.
The explicit source models do not accept unrelated modern calibration/domain
options. Raw BBP public dispatch requires exactly26 points.

The preparation script checks:
- exact original baseline catalog SHA;
- 32 available original core source records;
- exact SHA of each primary numerical source input;
- generation source hashes and native readback;
- distinct primary-target versus derived-sensitivity identity;
- 896 Cartesian cells, of which672 primary-target and224 derived-sensitivity.

BPS whole-chain views are tried for six bridge methods; exact-Maxwell uses
individual source phases because its existing two-branch kernel requires strict
chemical-potential ordering. Phase labels are not numerically interpolated.
The existing bridge matcher copies retained source labels and uses same-phase
endpoint derivatives.

## Executed independent checks

Command:
```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=python_frontend \
TMPDIR=/tmp/ns_eos_primary_crust_20261123/review \
python3 -m pytest -q -p no:cacheprovider \
  --basetemp=/tmp/ns_eos_primary_crust_20261123/review/pytest_primary \
  python_frontend/tests/test_bps_primary_native.py \
  python_frontend/tests/test_bbp_primary_native.py \
  python_frontend/tests/test_nv_primary_native.py \
  python_frontend/tests/test_primary_crust_dispatch.py
```

Observed result: 73 passed in11.25s.
The tests compile actual modules, perform independent source/reference and
high-precision derivative comparisons, and exercise malformed-source,
domain, sample-count and diagnostic-output failure contracts.
Log: `/tmp/ns_eos_primary_crust_20261123/review/primary_pytest.log`.

I inspected the completed campaign receipts independently:
896 cells =32 cores*4 distinct crust representations*7 methods;
547 native numerical survivors and349 no-admitted-match results.
Input/binary identities are recorded unchanged in the completed summary.
This is the native kernel layer, before source-admission qualification.

### Actual BBP retained-prefix archive check

All137 raw-BBP native survivor archives were opened, archive SHA256 verified,
and their first two retained n/epsilon/P rows compared to the original input
after the stated unit conversion and each output's own mass-header decoding.

Observed maximum relative errors:
- n:0;
- epsilon:2.220446049250313e-16;
- P:0.

No archive hash or prefix-preservation error was found. All137 therefore retain
the exact source obstruction. The six bridge-method counts are22 for the
first-law constrained Hermite pressure bridge and23 for each other bridge;
no raw-BBP exact-Maxwell survivor exists.

### Actual BPS phase-jump archive check

All132 BPS native survivor archives were opened and retained101-112 phase changes
checked. Across retained Gibbs pairs:
- maximum relative pressure difference:0;
- maximum absolute chemical-potential difference:3.410605131648481e-13MeV;
- maximum |delta epsilon - mu*delta n|:4.597017211338539e-17MeV/fm^3.

109 survivors used the primary equilibrium chain;23 used phase000 only.
Retained transition counts were:
3 transitions:4 survivors;
4:12;
5:36;
6:31;
11:26;
0 (phase000):23.

Consequently a surviving BPS bridge must not automatically be described as retaining
the entire outer-crust isotope sequence. Interface trimming remains part of the
construction; only26 inspected survivors retained all11 BPS transitions.
This is a coverage qualification, not a transport defect.

Evidence:
`/tmp/ns_eos_primary_crust_20261123/review/retained_prefix_and_bps_transitions.json`.

## Remaining unverified scope

- Root's final source-qualified admission/summary script was not yet available
  at this review snapshot.
- This review did not re-run the complete896-cell campaign, TOV/Cowling evolution,
  or the full repository test suite; it inspected actual fresh campaign evidence
  and independently executed the73 focused primary tests.
- Full historical BPS candidate mass catalogue, BBP coefficient-level CLDM and
  original NV Hartree-Fock solve remain scientifically unresolved/out of scope.
- Archived rejected endpoint-search outcomes are bounded finite-grid results;
  they do not prove nonexistence of every alternative continuum bridge.
- No claim is made that numeric derivatives of the NV fit remove its physical
  approximation or that BBP pressure-primary energy shifts are source rounding.

Exact reviewed-file hashes:
`/tmp/ns_eos_primary_crust_20261123/review/reviewed_file_sha256.json`.
