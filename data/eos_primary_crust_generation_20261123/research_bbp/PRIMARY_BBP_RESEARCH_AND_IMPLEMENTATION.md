# BBP primary-data recovery and native construction — 2026-11-23

## Scope and identity

The original BBP *model* is Baym, Bethe & Pethick, “Neutron star matter,”
Nuclear Physics A 175 (1971), 225–271, DOI
https://doi.org/10.1016/0375-9474(71)90281-8.
The publisher and INSPIRE agree on this identity.
The full original BBP functional has **not** been reproduced in this change.

A primary-author thermodynamic table was recovered from Baym, Pethick &
Sutherland, “The Ground State of Matter at High Densities: Equation of State and
Stellar Models,” ApJ 170 (1971), 299–317, DOI
https://doi.org/10.1086/151216.
The original paper explicitly attributes the density range 4.3e11–5e14 g/cm^3
to BBP on printed page 307. Table 5 continued, printed page 309 (PDF page 11),
contains the 35 BBP rows transcribed in BPS1971_TABLE5_BBP_CGS.dat.
This is a primary-author table, not the Diener2008 secondary transcription.

Original public source:
https://adsabs.harvard.edu/pdf/1971ApJ...170..299B
PDF SHA256:
91d0f9d04b0cf88b5c88a48c826815c3b55397fb17a2daf088fa52cc53bd215c

The PDF and rendered page are private scratch research inputs:
- /tmp/ns_eos_primary_crust_20261123/sources/bps/BPS1971.pdf
- /tmp/ns_eos_primary_crust_20261123/sources/bbp/BPS1971_p309.png

The whole page was visually inspected; OCR artifacts such as “7.8907229” were
resolved against the image to the printed token “7.890E29”. Numerical table data
alone are preserved in the repository. No full third-party paper is added here.

## Source–equation–code map

| Authority | Meaning | Implementation / status |
|---|---|---|
| BPS1971 p307, paragraph before Eq.28 | rho is total energy density divided by c^2 | rho -> epsilon with exact SI c; implemented |
| BPS1971 p307 Eq.28 | Gamma = (n/P) dP/dn = ((rho+P/c^2)/P) dP/drho | printed Gamma retained as diagnostic; not assumed to be a constraint on the chosen interpolant |
| BPS1971 Table5 p309 | 35 rows of rho, P, n, Z, A, Gamma | all 35 tokens transcribed and native parsing validated |
| Declared existing campaign domain n <= 0.08 fm^-3 | retain inner-crust support | 26 rows, n=0.000267 to 0.07468 fm^-3; no extension beyond last selected knot |
| Cold first law, derived here | d(epsilon/n)/dn = P/n^2 | analytical native pressure-integral constructor |
| Euler relation at T=0 | mu_b=(epsilon+P)/n | native calculation in both representations |
| BBP1971 original functional | compressible liquid-drop nuclear energy, neutron-gas energy, surface prescription and coefficients | not obtained or implemented; exact equation numbers/parameter table are unresolved |

Native module:
fortran_eos/src/compose_eos_bbp_primary.f90

Public interfaces:
- generate_bbp_primary_native(source_path,n_points,eos,output_dir,status)
- generate_bbp_pressure_primary_native(source_path,n_points,eos,output_dir,status)

Both return lowercase status='ok' only after input checks, EOS construction and
diagnostic writes succeed. This is implementation success, not scientific
admission of an EOS or matched core/crust pair.

## Units and semantics

The source rho is in g cm^-3, P in dyn cm^-2, and n in cm^-3.
Internal units are fm^-3, MeV fm^-3 and MeV. Use
c=299792458 m/s and 1 MeV fm^-3=1.602176634e32 J m^-3.
Thus n_internal=1e-39 n_source,
epsilon_internal=(1000 c^2/1.602176634e32) rho,
P_internal=(0.1/1.602176634e32) P_source.

Rest masses are already included in epsilon. The 939 MeV CompOSE header mass
is an encoding constant only, never an added physical rest energy. Z and A alone
do not determine the free-neutron fraction, so composition availability remains
false. The upper nine rows remain in the provenance table and are not relabeled
as inner crust. No historical source or matcher tolerance was modified.

## Representation 1: raw primary table

BBP1971_table5 retains the 26 source n, epsilon and P values after unit conversion.
Requested n_points is recorded but does not cause resampling of this fixed table.
This replaces a missing primary-data problem with an actually evaluable source;
it does not certify cold thermodynamic consistency.

There are four decreasing adjacent mu_b intervals in the selected raw table.
More strongly, the first interval cannot be joined by any stable cold barotrope
inside independent half-last-printed-digit rounding boxes. This is a statement
about that precise published finite table and that rounding assumption, not a
claim that the underlying BBP model is impossible.

For interval endpoints n1<n2, stable P(n) is monotone and the cold first law gives

  delta(epsilon/n) = integral[n1,n2] P(n)/n^2 dn.

Write each published value's closed rounding box as [value-half, value+half].
A conservative necessary bound allowing *simultaneous* changes in rho, P and n is

  Delta_min = epsilon2_low/n2_high - epsilon1_high/n1_low,
  Delta_max = epsilon2_high/n2_low - epsilon1_low/n1_high,

  I_min = P1_low (1/n1_high - 1/n2_low),
  I_max = P2_high (1/n1_low - 1/n2_high).

The two intervals must overlap. All selected density boxes are disjoint.
The first interval has, in MeV,

  [Delta_min, Delta_max] = [0.6088004576851063, 1.6494109687876062]
  [I_min, I_max]         = [0.2684345407837055, 0.2854525535769782].

The lower bound exceeds the upper allowed integral by 0.3233479041081281 MeV.
The gap is many orders above floating-point roundoff. Widening the source's
uncertainty, dropping its first row or changing epsilon could remove a numerical
obstruction but would change the stated source contract; none is done here.

As a second, narrower diagnostic, fixed-central-n/P piecewise-polytropic
integration permits no common epsilon/n integration constant within the rho
rounding boxes: [938.1014305288103,937.1374519576689] MeV is empty.
The standalone Fortran probe and FEASIBILITY.txt preserve both calculations.
The source itself does not state that half-last-digit rounding boxes contain its
true values; that assumption is explicit and only used for this diagnosis.

## Representation 2: explicitly derived pressure-primary EOS

BBP1971_pressure_integral is a different EOS identity. It keeps the tabulated
n_i and P_i, uses only the first source epsilon value as its energy anchor, and
derives all subsequent energies from the cold first law. It is **not** an
original unrounded BBP functional or a reconstruction consistent with all source
rho rounding boxes. Its discrepancy from the raw table remains recorded.

For n_i <= n <= n_(i+1), define

  gamma_i = ln(P_(i+1)/P_i)/ln(n_(i+1)/n_i),
  x = ln(n/n_i), H_i=ln(n_(i+1)/n_i),
  P(n) = P_i exp(gamma_i x),
  exprel(y) = (exp(y)-1)/y, exprel(0)=1.

The exact continuous first-law integral used in Fortran is

  e(n)=epsilon(n)/n
      = epsilon_1/n_1
        + sum_(j<i) (P_j/n_j) H_j exprel((gamma_j-1) H_j)
        + (P_i/n_i) x exprel((gamma_i-1) x).

A short stable series evaluates exprel near zero. At knots epsilon and P are
continuous, although dP/dn has one-sided values. No density extrapolation occurs.
For each smooth segment

  mu_b = e + P/n,
  (c_s/c)^2 = dP/d(epsilon) = gamma_i P/(epsilon+P).

All gamma_i are positive. For r=P/(epsilon+P), differentiation with respect to
ln n yields r'=r(gamma_i-1-gamma_i r). Hence r is monotone or constant on a
segment, and its endpoint values bound the sound speed throughout that segment.
Both sides of every knot are included in the reported bounds.

The output grid is the sorted union of n_points logarithmic nodes and all 26
original knots; n_points=257 produces 281 rows. Near-identical generated nodes
are merged into the exact original source knot within 64 machine epsilons.
This merges only generated sample locations; no original source values are
sorted, dropped, clipped or repaired.

## Diagnostic and rejection contract

Both constructors write:
- primary_source_values.csv: all 35 raw values after stated unit conversion,
  Z/A/Gamma, exact printed-decimal half steps and inner-domain flags.
- primary_source_anchor_comparison.csv: all 26 raw versus derived energy and mu
  values, signed differences and differences in rho half-rounding units.
- primary_source_diagnostics.json: representation, domain, selected/raw counts,
  first-law bound, raw mu reversals, first-anchor convention, derived sound-speed
  bounds and explicit non-claims.

The frontend seals those files with its other generated-source evidence.
Physical matching acceptance is assessed separately by the unchanged matcher
and verifier. A failed raw-source thermodynamic constraint does not automatically
reject the distinct derived construction.

The parser requires exactly 35 numeric rows, six whitespace-separated tokens,
positive finite physical data, scientific columns with the source's four-digit
mantissas, valid integer Z/A, finite positive Gamma, increasing rho/P/n and
disjoint density rounding intervals. Extra columns, rows, malformed records,
nonfinite tokens, overflow exponents, null list tokens and unsafe point counts
fail without returning a generated EOS. Source file comments and blank lines
are allowed; hidden inline fields are not.

## Reuse assessment

compose_eos_spg_native.f90 already eliminates baryon/charge constraints and a
surface/Coulomb radius in a modern GDFM spherical liquid-drop candidate. Its
surface function is explicitly calibrated input under Scurto Eqs.21–25.
compose_eos_tf_models.f90 provides GPPVA RMF kernels. Neither is the historical
BBP nuclear functional or neutron-density-dependent surface prescription.
Reusing those physical terms as “BBP” would be a model substitution. Their
optimizer, common EOS carrier and unit-aware I/O architecture remain reusable
once original BBP equations/coefficients are available.

## Validation completed

python3 -m pytest -q python_frontend/tests/test_bbp_primary_native.py

Result: 14 passed in 0.85 s.
The tests compile the actual Fortran module with -fcheck=all and
-ffpe-trap=invalid,zero,overflow. They check raw unit conversion and 26/35-row
scope, the source-incompatibility certificate, union-grid knot preservation,
an independent analytic integration oracle, positive sampled dmu, continuous
segment sound-speed bounds and malformed-input failure paths.
The test log is BBP_NATIVE_TESTS.log.
This is not a replacement for the parent-run all-core matching matrix.

## Exact literature/data request if full original microphysics is required

No additional user document is needed for the recovered primary table or the
two implemented native representations.

For a faithful coefficient-level BBP functional reproduction, provide either:
1. The full original 47-page paper, Nucl.Phys.A175 (1971),225–271,
   DOI 10.1016/0375-9474(71)90281-8, including all nuclear-energy, neutron-gas,
   surface-energy coefficients, equilibrium equations and any tabulated results;
   exact internal equation/page ranges could not be checked through the abstract.
2. An original-author unrounded n/epsilon/P table with units, energy-zero
   convention, domain and provenance. This would permit validation without
   interpreting the inconsistent published first interval as higher precision.

A legitimate publisher reprint exists as “NEUTRON STAR MATTER” in Selected
Works of Hans A Bethe, DOI 10.1142/9789812795755_0022.
Publisher fulltext requires access; the open abstract and INSPIRE record do not
supply coefficients. INSPIRE API record69606 has no attached documents or URLs.
SciSpace focused queries did not recover the original full text. No paywall was
bypassed, and no unknown equation numbering has been invented.

## Skill context

PhysMath Research Loop 1.3: main instructions, POLICY_INTEGRATION.md and
STAGE_GUIDE_KO.md. Assigned R2/R3/C2/C3 expanded by parent to C5/C6 for the two
explicit native representations. One bounded constructive pass; no child
review/delegation and no commit/push.
