# GM1 and BBB2 composition regeneration

GM1 now solves uniform nonlinear RMF with the original Glendenning--Moszkowski
Table II parameters and a common nucleon mass of 939 MeV. The provider npe
species set is the default, and an explicit muon option supports a separately
specified npemu calculation. The 939 MeV kernel mass and rounded coupling ratios
are declared reconstruction inputs. The CompOSE header reference masses alone
do not determine the original provider kernel convention.

For S=g_sigma*sigma, in inverse-fm units,
U(S)=S²/(2 C_sigma)+b m S³/3+c S⁴/4 and
dU/dS=n_s,n+n_s,p.
The vector energies are C_omega n²/2+C_rho(n_p-n_n)²/8.
The isovector vertex uses g_rho/2. The root search covers positive effective
mass stationary roots and ranks them by energy. Beta equilibrium and charge
neutrality determine the proton and lepton fractions.

Actual symmetric-matter saturation is n=0.15321152 fm^-3, binding=-16.35818 MeV,
m*/m=0.6999113 and K=300.644 MeV. Published rounded targets are
0.153,-16.3,0.7,300. Across 290 original GM1 core rows from
0.0829067 to 1.0 fm^-3, maximum pressure relative discrepancy is 0.4864%,
and maximum energy-density relative discrepancy is 0.06161%.
The retained comparison includes every density, pressure, energy and proton fraction.

BBB2 uses the Paris+Urbana symmetric/neutron energies in Baldo et al. Table 1,
with E(n,y)=E_SNM(n)+[E_PNM(n)-E_SNM(n)](1-2y)^2, as in Eqs.16--22.
PCHIP of the rounded table values is an explicit numerical representation.
The common neutron/proton mass approximation is retained. The npemu equilibrium
proton fractions agree with all 16 Table 3 targets to maximum absolute 0.000978.

BBB2 pressure differs from the published beta-equilibrium Table 3 by as much as
65.18%, at n=0.16 fm^-3. The original Paris columns, Table 3 identification and
pressure conversion (10^34 dyn cm^-2 = 6.241509074 MeV fm^-3) were independently
rechecked. The energy derivative and Euler identity of this implementation pass.
The original fitted functional coefficients or BHF solver are needed to close the
pressure reproduction discrepancy. This result remains recorded as
PRESSURE_MISMATCH_REQUIRES_ORIGINAL_FIT. The paper's microscopic calculation
extends to 0.8 fm^-3, with its rational extrapolation supplying higher-density
table points through 1.5 fm^-3.

Both generated 81-row grids export n,p,e,mu fractions in eos.compo and pass native
readback, charge/baryon closure and Euler checks. Sampled mechanical/causal
diagnostics cover the generated GM1 0.08--1.0 and BBB2 0.08--1.5 fm^-3 grids.
The BBB2 source-pressure discrepancy is retained alongside those internal checks.
There are 21 passing targeted regression tests. The initial q_extra shape
serialization error and its correction are preserved in the campaign logs.

Evidence is in data/research_generator_20261116/reference_cores, with complete
CompOSE tables, states, compressed receipts and comparisons, checksums and logs.

Sources:
- https://escholarship.org/uc/item/6sw424k1 (GM1 original Table II)
- https://compose.obspm.fr/eos/54 (GM1 provider)
- https://arxiv.org/pdf/astro-ph/9707277 (BBB Tables 1,3 and Eqs.16--22,29)
- https://compose.obspm.fr/eos/69 (BBB2 provider; species description discrepancy
  against the original paper remains documented in the literature intake)

## BBB pressure discrepancy follow-up

The CompOSE BBB2 core has 15 density nodes at 0.16--1.5 fm^-3.
Every provider pressure matches the original paper's Table 3 after unit conversion
to relative 1.07354e-7, consistent with a small conversion-constant difference.
Thus the pressure comparison directly tests the provider target as well as the paper.

Let e(n,y) be baryonic energy per particle excluding the common rest mass,
and let y(n) be the equilibrium proton fraction. Write
epsilon=n[m+e(n,y)]+epsilon_l. Neutrality gives n_e+n_mu=ny, while beta equilibrium
gives e_y=-mu_e and mu_e=mu_mu for populated species. Therefore

d epsilon/dn = m+e+n e_n+n e_y y' + mu_e(y+n y')
             = m+e+n e_n+mu_e y,

and the pressure consistent with the total equilibrium energy is

P_thermo = n² e_n + P_l.

Differentiating only the baryonic equilibrium trajectory and then adding the
ordinary lepton pressure instead yields

P_baryon_path = n² [e_n+e_y y'] + P_l
              = P_thermo - n² mu_e y'.

The original paper's Eq.23 is written as P_nucl=n² d(E/A)/dn at fixed A,
after discussing the baryonic energy in beta equilibrium, and Eq.24 adds P_lep.
This wording permits a composition-path interpretation. Its use inside the
authors' implementation remains an inference, assessed numerically below.

Using CubicSpline as an independent smooth interpolation control, P_baryon_path
matches all 15 provider-density Table 3 pressures within 0.6413%, whereas
P_thermo differs by up to 68.30%. At n=0.16 fm^-3, P_thermo=2.55260,
the composition-path term is 1.03786, and P_baryon_path=1.51474 MeV fm^-3,
compared with the published 1.51669. At n=0.8 fm^-3, corresponding values are
249.15090,36.95611,212.19479 and 212.21131 MeV fm^-3.
The excluded lower endpoint at 0.08 fm^-3 has 8.85% path-pressure discrepancy.
PCHIP's path-pressure discrepancy is larger, reaching 16.07%, demonstrating
additional interpolation sensitivity around the low-density sparse grid.

A 200-realization, fixed-seed perturbation of every Table 1 energy within half
its final displayed digit quantifies finite printed-precision sensitivity.
At n=0.16, the CubicSpline thermodynamic pressure discrepancy spans
48.85%--87.51%, while the path-pressure discrepancy spans -9.27%--9.60%.
At n=0.8, the corresponding ranges are 17.331%--17.495% and
-0.0837%--0.0657%. These sampled ranges distinguish the sizable composition-path
term from rounding effects at those densities. The 0.08 endpoint is poorly
constrained by printed precision and interpolation.

The current default generator retains P_thermo. The follow-up supplies both
labelled pressures and their difference as diagnostics. It supports a specific
source-convention hypothesis; the authors' exact differentiation procedure and
smooth-fit coefficients remain the concrete inputs needed to settle attribution.
A blanket model rejection would discard the useful composition agreement and
the density-resolved information established here.

Reproduce the diagnostic with:

~~~sh
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran \
PYTHONPATH=python_frontend python3 \
data/research_generator_20261116/reference_cores/diagnose_bbb_pressure.py \
build/next_bbb_pressure_check --source build/domain_research_final/core_only/BBB_BHF2
~~~

The retained pressure_followup.json contains both interpolation controls,
all source comparisons and the printed-precision sample bounds.
