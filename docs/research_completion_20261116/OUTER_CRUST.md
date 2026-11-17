# Independent BSk22 outer-crust regeneration

Scope: cold, catalysed, fully ionized, single-nuclide bcc crystal from baryon
density 1e-9 fm^-3 up to neutron drip. This is an independent nuclear-mass
calculation with composition. No CompOSE thermodynamic or composition row enters
the generator. The common BSk22 segment in the thirteen original crust archives
continues beyond neutron drip: its remaining ETF/Strutinsky inner-crust interval
is a separate obligation. No outer-inner or core-crust joining is done here.

## Primary sources and intake

- Pearson et al., MNRAS481,2994(2018), Sec2.1, AppendixB and Table3:
  https://arxiv.org/abs/1903.04981 (paper plus correction).
- Pearson, Goriely, Chamel, PRC83,065810(2011), Eqs16-22:
  https://doi.org/10.1103/PhysRevC.83.065810.
- Author HFB22 mass/radius file:
  https://www.astro.ulb.ac.be/bruslib/nucdata/hfb22-dat.
- AME2016 atomic mass table:
  https://www-nds.iaea.org/amdc/ame2016/mass16.txt.
- Welker et al., PRL119,192502(2017), TableI:
  https://doi.org/10.1103/PhysRevLett.119.192502.

The intake script checks exact SHA256 identities before parsing HFB22 and AME2016.
It reads HFB22 column Mcal as atomic mass excess in MeV and Rch as rms charge
radius in fm. All AME2016 measured entries without the '#' extrapolation flag
override the corresponding HFB22 masses; HFB22 supplies radii. The five Cu masses
in Welker2017 TableI override both. For78Cu the Penning-trap result is used and
the different MRTOF result is recorded explicitly. HFB22 masses have the published
0.01MeV rounding; no fitted missing digits are invented.

The merged catalogue has8481nuclides, including89measured nuclei absent from
the HFB22 radius table. Such candidates receive radius0 only as a lower bound
on their positive finite-size energy. They participate in the competition, and
the solver refuses if any wins. On the supported fixed-isotope equilibrium branch,
pointwise lower energy implies a lower Gibbs minimum; a missing-radius candidate
that loses even at this lower bound cannot be rescued by restoring its radius.
This excludes the measured lower-bound candidates in the executed census, not
all possible unknown nuclei outside the mass inventory.

Exact identities and counts are in data/outer_crust/provenance.json. A previously
acquired public Chamel2020 solver input is HFB27, not HFB22, so it was not silently
substituted for this input. The independent HFB22 source was successfully acquired.

## Equation-to-code contract

Let n_e be electron number density, n_i=n_e/Z and n_b=A n_i. Atomic masses
are converted to electron-binding-corrected atomic energies M'c^2 by Pearson2018
Eq5. The nuclear rest term is n_i(M'c^2-Z m_e c^2). Relativistic electron energy
includes the electron rest energy, so that rest mass is counted exactly once.

The implemented energy is

    epsilon = n_i (M'c^2-Z m_e c^2)
            + epsilon_free,e + epsilon_exchange + epsilon_correlation
            + epsilon_screening + epsilon_bcc + epsilon_finite_size + epsilon_zp.

The first term is independent of density per ion. Free electrons use the exact
zero-temperature relativistic integrals. Exchange uses Pearson2018 B9 at all
relativities, including its sign change; it is not replaced by the high-density
limit. Correlation uses B11 with log10(x_e). Screening uses B7-B8, including the
full b2 numerator 1+0.684 Z^5+0.222 Z^6 and the A-dependent ion plasma term.
The plasma mass follows M' in Pearson2018 B6. This convention is explicit because
an ion-mass comment in the source might otherwise invite an unnoticed substitution.

The bcc Madelung term has coefficient0.895929255682. Finite-size and zero-point
terms follow Pearson2011 Eqs16-22 with the HFB22 rms charge radius and the bcc
phonon average0.5113875. Their energy-density powers at fixed isotope are n_e^2
and n_e^(3/2), so their pressures are respectively epsilon_finite_size and
epsilon_zp/2. Screening pressure differentiates the entire density dependence:
Gamma_pl is proportional to n_e^(-1/6), and its contribution cannot be dropped.
Every pressure term is obtained from P=n_e*d(epsilon)/dn_e-epsilon at fixed A,Z.

All calculations use MeV and fm with hbar*c=197.3269804MeVfm. The neutron rest
energy is939.56542052MeV, proton938.27208816MeV, electron0.51099895MeV and atomic
mass unit931.49410242MeV. These differ from the SPG candidate's common938.9MeV
reference; physical joining requires a separate mass-convention decision.

## Discrete equilibrium and composition

At a common specified pressure the solver finds a separate electron-density root
for every nuclide and minimizes g=(epsilon+P)/n_b over the entire catalogue.
It does not minimize energy across nuclides at common density. Roots use bounded
vector bisection; the actual pressure residual is retained. The selected isotope
has zero free-neutron and free-proton fractions, electron fractionZ/A and cluster
baryon fractionYh=1. Cluster number fraction1/A is a separate field.

If a third nuclide has lower Gibbs energy at an apparent two-nuclide crossing, the solver recursively splits the bracket through that intermediate phase. The initial257-point census exposed this case and its failure log is retained. A16-point/five-nucleus regression exercises the repaired path.

The native table contains both one-sided density limits at each first-order
nuclear transition. A common transition pressure is stored; the density jump is
preserved. Each isotope receives a separate Iphase branch code so derivative
stencils do not cross these jumps. Q4 uses the electron chemical potential from
the complete electron energy/pressure terms (Pearson2018 B24), while Q3 uses g.
The table retains beta-equilibrium Q5=0, zero entropy and zero temperature.

Neutron drip is solved as min g(P)=m_n c^2. Above that point free neutrons are
required and this outer-crust model refuses ordinary state queries. An internal
allow_above_drip option exists only to bracket the drip root. Its continuation is
not exported as a physical inner crust. The default lower boundary is the pressure
of56Fe at n_b=1e-9fm^-3; no atomic atmosphere or physical zero-pressure surface is
claimed.

## Theory and coding research loop

| Step | Concrete output and assessment |
|---|---|
|1Problem|Discrete outer-crust generation identified separately from the post-drip common segment.|
|2Sources|HFB22, AME2016 and Cu primary inputs acquired, with SHA identities and precedence.|
|3Model|Exact source electron, bcc, finite-size and zero-point energy terms transcribed.|
|4Derivation|All pressure terms differentiated at fixed isotope, including screening plasma density dependence.|
|5Limits|Exact relativistic exchange retained; computational bounds and neutron-drip support explicit.|
|6Implementation|Global finite-catalogue minimization at fixed pressure and explicit composition/native export.|
|7Falsification|Independent60-digit mpmath source energy plus automatic differentiation found two failing test cases from an omitted b2 numerator term.|
|8Repair|Full source b2 restored; a coarse-grid third-winner counterexample led to recursive intermediate-layer refinement. Both red logs remain and14targeted tests pass.|
|9Validation|257/513pressure census comparison, refined nuclear transitions, all-candidate compressibility samples, native write/read, and source Table3 diagnostics.|
|10Claim|Regenerates the scoped outer crust; complete original13crust regeneration still requires the post-drip BSk22 segment and each distinct inner-crust implementation.|

The independent energy oracle differentiates the scalar source energy with
mpmath; the production code's analytic pressure is not reused as its expected
answer. Tests also cover measured-mass precedence, common-pressure global finite
catalogue minimum, neutron drip, above-drip refusal, composition closure and
missing-radius refusal. Red/green logs preserve the actual initial source-transcription
failure rather than claiming a clean first implementation.

Run commands:

```sh
python3 scripts/intake_outer_crust_masses.py
PYTHONPATH=python_frontend python3 -m pytest -q python_frontend/tests/test_outer_crust.py
PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. \
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran \
python3 scripts/build_outer_crust.py OUTPUT_DIRECTORY
```

The reproducible campaign receipt is
data/research_completion_20261116/outer_crust/generated/RECEIPT.json. It retains
the full isotope sequence, transition roots, neutron-drip state, source and code
identities, native readback errors and consistency diagnostics. Rounded Table3
pressures/densities are comparison diagnostics; neither that table nor CompOSE
agreement is an unconditional admission gate. Agreement of two finite censuses
does not prove that every arbitrarily thin layer has been excluded.

## Executed result and limitations

The final campaign generated547rows and18nuclear layers. The257- and513-point
pressure censuses give the same layer sequence; refined transition pressures agree
within5.414e-11relative. One coarse interval required adaptive intermediate-layer
refinement. All16native readback arrays are identical in this execution, and native
C1-C6/C10/C11 pass. C7-C9 remain finite-difference diagnostics: their maximum relative
residuals are respectively0.01125,0.00658,1.193e-5. These mesh derivatives are distinct
from the independent60-digit local energy/pressure tests.

The final nucleus is128Sr. Neutron drip occurs at n_b=2.694320064010902e-4fm^-3,
P=4.992059167689819e-4MeVfm^-3 and epsilon=0.2526497904790187MeVfm^-3.
The computed g=939.5654205200001MeV equals the stipulated neutron rest energy.
The full nuclear sequence matches Pearson2018Table3; maximum deviations from its
rounded listed layer endpoints are0.367% in pressure and0.456% in density. These
numbers are source diagnostics and do not establish unpublished-mass precision.

The alternative78Cu mass lowers its energy per baryon by0.000603MeV, while its
minimum nominal gap over the65sampled pressures is0.03549MeV/baryon. This is a
finite-sample sensitivity check. General nuclear-mass uncertainty propagation,
including the quoted0.105MeV uncertainty in79Cu, has not been performed.

The mass CSV is currently loaded from the repository checkout. Editable installation
and repository execution are validated; distributing these data in a non-editable
wheel remains packaging work. The independent oracle dependency mpmath>=1.3 has
been added to the test extra and the separate micromamba native-test environment.
