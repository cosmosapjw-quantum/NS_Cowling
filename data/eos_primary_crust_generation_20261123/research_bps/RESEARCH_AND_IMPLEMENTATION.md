# BPS original-source generation, 2026-11-23

## Scope and scientific identity

Implemented `BPS1971_sequence`: the BPS1971 electron + bcc lattice equations,
minimized at fixed pressure over the twelve nuclei and binding energies printed
in that paper's Table 2. This is a primary-equation calculation with a
published-sequence-restricted candidate catalogue. It is not a reproduction of
every candidate of the original UCRL-11980 nuclear mass catalogue, a modern
AME/HFB calculation, or a numerical interpolation of a secondary BPS-labelled
table. Production EOS and verification mathematics are Fortran.

Primary paper:
G. Baym, C. Pethick, P. Sutherland, The Astrophysical Journal 170, 299-317 (1971),
DOI https://doi.org/10.1086/151216 .
Public primary scan: https://adsabs.harvard.edu/pdf/1971ApJ...170..299B .
DOI metadata independently verified at https://inspirehep.net/literature/69636 .

Private downloaded PDF:
`/tmp/ns_eos_primary_crust_20261123/sources/bps/BPS1971.pdf`
SHA256 `91d0f9d04b0cf88b5c88a48c826815c3b55397fb17a2daf088fa52cc53bd215c`.
The PDF is not committed. Tables 2 and 5 were transcribed and checked against
rendered page images, including a higher-resolution Table 5 crop. They contain
numerical facts; no third-party full text is included in this research directory.

## Primary equations and conventions

At zero temperature and local charge neutrality, let n_i=n_b/A and n_e=Z n_i.
W_N=(A-Z)m_n c^2 + Z m_p c^2 - A b is the nuclear energy used by the paper.
The source binding b retains atomic-electron binding energy while electron rest
mass is subtracted from nuclear masses (discussion p302; definition Eq19 p305).
Adding a modern atomic binding correction would change this specific model.

With a=(2/n_i)^(1/3), the original bcc lattice energy per nucleus is
W_L=-1.819620 Z^2 e^2/a (Eq3). Define
x=(hbar c)(3 pi^2 n_e)^(1/3)/(m_e c^2), gamma_e=sqrt(1+x^2), and
K=(m_e c^2)^4/[pi^2 (hbar c)^3]. Then

- epsilon_e=K/8 [x(1+2x^2)gamma_e-asinh(x)].
- P_e=K/24 [x(2x^2-3)gamma_e+3 asinh(x)].
- epsilon=n_i W_N+epsilon_e+n_i W_L.
- P=P_e+n_i W_L/3.
- mu_e=m_e c^2 gamma_e.
- mu_b=[W_N+(4/3)W_L+Z mu_e]/A=(epsilon+P)/n_b.

The minimum mu_b at common P is the stable supplied-candidate state. Every
adjacent phase change solves equality of the two Gibbs energies at the same
pressure. Both one-sided n_b and epsilon are retained. The baryon chemical
potential is continuous; n_b and epsilon jump.

For each homogeneous isotope branch, direct differentiation gives
dP/dln(n_b)=n_e m_e c^2 x^2/[3 sqrt(1+x^2)] + (4/9)n_i W_L.
This defines Gamma=(dP/dln n_b)/P and c_s^2/c^2=(dP/dln n_b)/(epsilon+P).

A small-x convergent series avoids cancellation in the electron integrals.
Every generated state is checked for finite variables, positive pressure and
0<c_s^2/c^2<1. This is a branch stability check; coexistence endpoints explicitly
encode the first-order transitions.

Units are fm, MeV, and MeV/fm^3. Constants are:
hbar*c=197.3269804 MeV fm; m_e*c^2=0.51099895 MeV;
m_n*c^2=939.56542052 MeV; m_p*c^2=938.27208816 MeV;
alpha=1/137.035999084; e^2=alpha*hbar*c.
CGS conversion uses 1 MeV/fm^3=1.602176634e33 dyn/cm^2 and
epsilon/c^2=1.7826619216278976e12 g/cm^3.
The constants are explicit contemporary dimensional constants; the published
b values are historical, rounded to 0.0001 MeV per baryon. Exact historical
byte/numeric reproduction is not asserted.

Neutron drip uses the paper's stated approximation mu_b=m_n*c^2, neglecting the
small continuum-neutron lattice shift. Low support starts at n_b=6.295e-12 fm^-3,
the first Table 5 row of BPS's own calculation (rho about1.044e4 g/cm^3).
The Feynman-Metropolis-Teller part below this density is not extrapolated.
Muon, electron exchange/correlation/screening, finite-nucleus and zero-point
terms of later crust models are absent from this explicitly historical branch.

## Original numerical data and independent validation

`bps1971_table2.csv` contains the 12 source isotopes:
56Fe,62Ni,64Ni,84Se,82Ge,80Zn,78Ni,76Fe,124Mo,122Zr,120Sr,118Kr.
The final jump-percent field is -1 as a documented missing-value marker at drip;
all other entries are direct printed Table 2 facts.

`bps1971_table5_outer.csv` contains 43 rows of the original Table 5
'Present Calculation' section on p308, columns rho,P,n_b,Z,A,Gamma.
rho means total energy density divided by c^2, not a baryonic rest-mass proxy.
n_b converts from cm^-3 to fm^-3 by 1e-39. The table is an independent test input,
not the production generation source.

Native comparisons at each source density and specified isotope produce:
max relative rho error 5.537246063685e-4;
max relative pressure error 4.209205129047e-4;
max relative Gamma error 3.278473934829e-4.
Predeclared comparison limits are rho1e-3, P2e-3, Gamma2e-3 to include the
printed four-digit n,rho,P,Gamma values and explicit constants convention.
These limits are limited to reference-table comparison and do not replace or
relax any project matching/source gate.

Four original pure-isotope rows are inside exact Maxwell coexistence intervals.
Their n_b values are 5e-9,7.924e-7,5e-6,9.976e-5 fm^-3.
The Table 5 isotope label is therefore not the stable fixed-pressure pure phase
at that row's pressure. Each discrepancy was independently classified by
checking that its n_b lies strictly between the computed one-sided limits and
its isotope is one of those two phases. There are four such cases and zero
unclassified composition discrepancies. Raw source numbers are unchanged.
This distinction is required when comparing the paper's tabulated pure-phase
branches with a thermodynamically convex fixed-pressure representation.

For n_points257 the actual output is291 rows, with12 distinct phases and11
coexistence endpoint pairs. Native checks:
finite-difference first-law relative residual1.925193338081499e-11;
pressure inversion relative residual3.9324099532223045e-13;
Euler identity residual4.5474735088646412e-13 MeV;
maximum transition delta-mu2.2737367544323206e-13 MeV.
Neutron drip:
n_b=2.5869979226819602e-4 fm^-3,
rho=4.3242829229315985e11 g/cm^3,
mu_e=26.190774267308036 MeV.
The original paper's last table row is rho4.299e11 and n2.572e-4; this generated
drip boundary is not claimed to reproduce that row exactly. It follows the
explicit original equation with the published rounded b and stated constants.

## Code map and executable contract

Existing `compose_eos_lower_crust_native.f90:isotope_ep` implements Pearson
2018/2019 with additional electron/lattice corrections and modern AME/HFB masses.
It was intentionally left unchanged because it is a different physical model.

New `fortran_eos/src/compose_eos_bps_primary.f90` exports:
`generate_bps_primary_native(source_path,n_points,eos,output_dir,status)`.
The source is Table 2 CSV; output_dir must exist; 16<=n_points<=100000.
Success is lowercase 'ok'. Error paths return a precise status without STOP.
The generated eos_data_t includes complete composition masks, free n/p fractions0,
heavy-cluster baryon fraction1, Ye=Z/A, and isotope-specific phase IDs101-112.

The generator samples a global logarithmic P grid, inserts all exact coexistence
pairs and an interior point of every phase, and checks for an intervening
supplied candidate at every transition and phase interior. Transition positions
are independent of n_points. Only O(N+12) states are stored; sorting starts from
an ordered global grid and appends only34 extra rows. There is no shared mutable
module state; independent calls are reentrant.

Output diagnostics:
`bps_states.csv`, `bps_transitions.csv`, `primary_source_diagnostics.json`.
Root-owned application/build/frontend plumbing integrates these files.

Native standalone verification driver:
`data/eos_primary_crust_generation_20261123/research_bps/bps_primary_probe.f90`.
Pytest transport wrapper:
`python_frontend/tests/test_bps_primary_native.py`.
It compiles the actual Fortran module, invokes all calculations there, checks
printed evidence and JSON, and performs no production EOS mathematics.
16tests passed, including independent43-row verification for16/257/513 samples,
transition-byte equality across16/513 samples, malformed source input,
extra/trailing CSV values, missing/extra rows, wrong isotope order, nonfinite
input, insufficient samples, comment annotations and output I/O failure.
See PYTEST.log and NATIVE_PROBE_257.log for actual captured output.

## Remaining scientific distinction and papers

No BPS paper upload is needed: the original1971 paper and necessary sequence
binding data were accessible. A full repetition of the original Table1 candidate
search additionally needs Myers & Swiatecki, UCRL-11980 (27May1965), including its
mass tables. BPS refers specifically to that report. A repository metadata page
was found at https://escholarship.org/uc/item/65v5x6mx but the full report tables
were not recovered in this bounded phase. This does not block the declared
Table2-restricted model.

Related published mass-model paper: Myers & Swiatecki, Nuclear Physics81 (1966)
1-60. Publisher records expose DOI10.1016/0029-5582(66)90639-0 and a parallel
record DOI10.1016/S0029-5582(66)80001-9. Neither is asserted to be a DOI for the
395-page UCRL-11980 report. Recovering that report would permit a new,
separately-provenanced all-candidate reconstruction; modern nuclear masses are
not used to fill those missing historical data.

The BPS Table5 continued on p309 also supplies n_b, rho, P for the BBP equation
of state. That evidence and the shared private PDF path were handed to the
BBP specialist; this module does not mix BBP inner matter into BPS outer matter.
