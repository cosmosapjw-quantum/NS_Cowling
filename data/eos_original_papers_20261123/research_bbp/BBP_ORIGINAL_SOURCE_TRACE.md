# Original BBP functional: derivation, implementation, and evidence

Date: 2026-11-23. This component is part of branch
`research/eos-original-papers-20261123`, based on
`0fbeeda62f078209cad03798c435933c34e2e1a1`.

## Source authority and scope

G. Baym, H. A. Bethe, C. J. Pethick, *Neutron star matter*, Nuclear Physics A
175 (1971), 225–271, DOI **10.1016/0375-9474(71)90281-8**.
Publisher metadata: <https://www.sciencedirect.com/science/article/pii/0375947471902818>.
The user-provided 47-page PDF is authoritative for equations and printed tables.
Its SHA256 is
`f5d15dbc97c8dae24557b55844c2bfc550984aaa50e514e8d9412e11a05ca495`.
The complete paper was read; decisive equations on printed pp.235–237 and
245–252, and Table7 p.259, were checked against rendered page images.
The PDF itself is not placed in the public source tree.

The previous original-functional source-access blocker is **closed**. No further
paper is necessary to execute the original spherical inner-crust model on the
declared support. Exact reproduction of all rounded published table cells remains
a separate numerical provenance issue; no missing paper is claimed to solve it.

Implemented scope is cold, beta-equilibrated, locally neutral spherical nuclei,
uniform neutron gas and ultrarelativistic electrons, with continuous A and Z.
The generator uses `2.67e-4 <= n_b <= 0.08 fm^-3`. Shell effects, pairing,
nuclear deformation, proton drip, pasta, and the paper's high-density liquid
instability analysis are outside this generator. This is the domain that was
needed by the pending core–crust calculation. A lower energy than the uniform
npe trial is checked; that comparison is not a proof of a global minimum over
all possible inhomogeneous geometries.

## Equation-to-code map

| Source | Physical content | Native implementation |
|---|---|---|
| Eqs3.8,3.11; p235 | Pure neutron W and chemical potential | `neutron_w`, `neutron_dw`, `bbp_original_bulk` |
| Eq3.9; p235 | Proton impurity chemical potential | `bulk_w`, `bbp_original_bulk` |
| Eqs3.15,3.18,3.19; pp236–237 | Kinetic, symmetric-matter and isospin interpolation | `bulk_w` |
| Eq3.21; p238, p249 | w0=16.5MeV, k0=1.43fm^-1, K=143MeV, s=33MeV, sigma=21MeV | Named constants and strict parameter source |
| Eq4.30; p245 | Neutron-gas-dependent surface energy | `bbp_original_point` |
| Eq5.3; p246, thickness convention p248 | Finite-thickness Coulomb correction, d=0.74/k_c | `bbp_original_point` |
| Eq5.9; p247 | Cell Coulomb energy including lattice and finite-size correction | `bbp_original_point` |
| Eq5.13; p248 | Ultrarelativistic electron energy | `bbp_original_point` |
| Eqs6.3–6.4; p249 | Exact elimination of A, surface energy twice Coulomb energy | `a=ws/(2*wc*x**2)` |
| Eqs2.4–2.5; p230 | Total energy and conserved baryon density | Reduced energy in `bbp_original_point` |
| Eqs2.20,6.13–6.18; pp232,250–251 | Pressure and equilibrium conditions | Constrained derivatives; independent direct-pressure test |
| p251, paragraph before Table1 | Deliberately omitted gas derivative in the surface term | `wgs%d=0` only for the published closure |
| Tables1–3; pp251–252, Table7 p259 | Independent published numerical checkpoints | Dedicated tests and benchmark CSV |

All generator physics, gradients, stationarity solves, first-law checks and the
Maxwell construction are native Fortran in
`fortran_eos/src/compose_eos_bbp_original.f90`. Python performs I/O orchestration,
build invocation and independent verification only.

### Conventions

`n_i=k^3/(1.5*pi^2)` and `n_g=k_g^3/(1.5*pi^2)`; the actual neutron Fermi
momentum in the gas is `2^(1/3)*k_g`. `x=Z/A`; `alpha=1-2*x`.
`W` excludes rest energy and is in MeV per nucleon. `n_b,n_i,n_g` are fm^-3;
pressure and total energy density are MeV fm^-3. Nucleon rest energies are
`m_n*c^2=939.56542052MeV`, `m_p*c^2=938.27208816MeV`.
`hbar*c=197.3269804MeV fm` and `e^2=1.43996454784255MeV fm` are explicit.
The paper's rounded neutron kinetic coefficient `19.74MeV fm^2` is retained
throughout its bulk interpolation, including the symmetric kinetic term.
These modern physical constants are declared, not fitted to the table.

The proton exchange term and electron screening are omitted as discussed on
pp246–248. Finite-thickness energy is included, as the authors state on p246.
This source convention yields close agreement with Tables1–3 without tuning
the quoted five nuclear parameters.

The central source formulas are

\[
W_n(k)=19.74k^2-\frac{k^3(40.4-1.088k^3)}{1+2.545k},
\quad
\mu_{p0}=-\frac{k^3(218+277k)}{1+8.57k^2},
\quad
\mu_{n0}=W_n+\frac{k}{3}\partial_kW_n.
\]

The denominator of Eq3.8 is linear in k. The interpolation Eq3.19 is even in
alpha and extends through alpha^6; OCR descriptions suggesting alpha^8 or a
quadratic denominator were not used. The tests enforce `W(k,0)=W_n(k)`,
`W_x(k,0)=mu_p0-mu_n0`, proton–neutron exchange symmetry, the saturation
energy and derivative, and `W_kk(k0,1/2)=K/k0^2`.

Let `u=(n_b-n_g)/(n_i-n_g)` be the nuclear volume fraction. Then

\[
w_s=\frac{\sigma}{w_0}(W_n-W_i)(1-n_g/n_i)^{2/3},\qquad
w_c=\frac35\left(\frac8{9\pi}\right)^{1/3}e^2k
 (1-\tfrac32u^{1/3}+\tfrac12u),
\]
\[
A=\frac{w_s}{2w_cx^2},\quad
k_c=(k^3-k_g^3)^{1/3},\quad d=0.74/k_c,\quad
W_{\rm thick}=-\frac{4\pi}{9}e^2d^2k^3x^2.
\]

After eliminating A, the non-rest energy density is

\[
\varepsilon-m_nc^2n_b =
u n_i[W_i+W_{\rm thick}+(m_p-m_n)c^2x
 +\tfrac32w_sA^{-1/3}]
+(1-u)n_gW_n+\tfrac34\hbar c(3\pi^2n_e)^{1/3}n_e,
\quad n_e=u n_ix.
\]

The independent solver variables are `(k,x,t)` with `t=n_g/n_b` at fixed n_b.
The reduced binding energy is `f=epsilon/n_b-m_n*c^2`. Forward-mode automatic
differentiation has four directions `(k,x,t,n_b)`. A bounded damped Newton
solve finds the three stationary conditions. The gas-free face is solved
separately because its true one-sided gas derivative is singular. Interior
solutions come from five deterministic starting points, and variational
candidates, including the separately solved gas-free face, require reduced
Hessian eigenvalues at least -1e-7 in the declared coordinates. The lowest energy
among accepted gas-free/interior candidates is retained. These are explicit
finite candidate checks, not a claim of a global optimization theorem.
The exact search box is `max(1.05,(1.5*pi^2*n_b*1.000001)^(1/3)) <= k <= 1.5`,
`.025 <= x <= .45`, `1e-9 <= t <= .999999`; the gas-free `t=0` face is separate.
The published closure takes the first converged gas-bearing root. The Maxwell
solve uses log-density seeds `(2.73e-4,3.0e-4)` and constrains
`n_b,min <= n_a < n_b <= 5e-4`. Neither the finite seed search nor the selected
common tangent certifies all spherical branches or their global convex hull.
The `true_g_t=0` entry on the gas-free face is an inactive-coordinate placeholder;
the original surface functional has a divergent one-sided gas derivative there.

## Three scientifically distinct outputs

| Model | Mathematical definition | Intended use |
|---|---|---|
| `BBP1971_original_published` | Original energy with the p251 omission in the surface derivative | Reproduce the authors' numerical closure; retain and report thermodynamic failure |
| `BBP1971_original_variational` | All derivatives of the same energy, pointwise lower energy of the spherical branches | Diagnose consequences of restoring the omitted derivative, including a nonconvex branch switch |
| `BBP1971_variational_maxwell` | Fully variational energy with native common-tangent phase construction | Thermodynamic completion of that approximate functional; explicitly derived, not a claim that BBP published this low-density transition |

For the published closure only `dW_n/dn_g` inside `w_s` is omitted. The neutron
gas bulk derivative and the derivative of the finite-thickness term remain.
The actual energy value is identical between closures at the same state.
The omitted contribution to `df/dn_g` is

\[
C=\frac{u n_i}{n_b}A^{-1/3}\frac{\sigma}{w_0}
(1-n_g/n_i)^{2/3}\frac{dW_n}{dn_g}.
\]

Thus the exact and published gradient differences are `Delta f_t=C*n_b` and
`Delta f_nb=C*t`, with zero differences in the k and x directions. These
identities are explicitly verified by an independent test. Along a smooth
published-closure branch,

\[
P_{\rm derivative}-P_{\rm published}
=n_b^2 C\frac{dn_g}{dn_b},
\]

so exact Euler encoding alone cannot establish the first law. The exported
published pressure is the selected constrained derivative; it is independently
verified against the direct paper expression `P=P_G+P_e+P_L`, with the same
p251 omission in Eq6.13. Its baryon chemical potential also agrees with
`m_n*c^2+mu_n^(G)`. The total derivative along the re-solved branch is tested
separately using `Delta n_b/n_b=1e-4` and no coarse-grid interpolation.

At 257 requested densities the published closure has max selected-gradient
residual `6.388e-11`, but max exact-energy gradient `0.164934MeV` and max relative
same-phase first-law defect **0.06507544**. This is a measured property of the
specified approximation. It is not a no-go theorem about BBP microphysics or
an inference that a generator cannot be built.

The fully variational pointwise path selects four gas-free points on this grid,
then changes to the gas-bearing phase. The raw pressure falls from approximately
`5.68201e-4` to `5.36821e-4MeV fm^-3` across that switch. These raw results are
retained under their own identity. The phase switch and its nonconvexity are
not repaired by adding a numerical pressure slope.

### Native Maxwell completion

The common-tangent solver varies both endpoint densities, retaining the
one-sided stationary phases, and solves `mu_a=mu_b`, `P_a=P_b`. It finds

- `n_a=2.7622175767948561e-4 fm^-3`, gas-free phase101;
- `n_b=3.1566706457246754e-4 fm^-3`, gas-bearing phase102;
- `P=5.4672144575670e-4 MeV fm^-3`;
- `mu=939.8883473106335MeV`.

Before representing common values, the endpoint differences are
`|Delta P|=4.21e-15MeV fm^-3`, `|Delta mu|=9.21e-12MeV`.
The unadjusted endpoint values are saved in `bbp_original_transitions.csv`.
At the output endpoints, common P and mu are represented exactly and epsilon
is their common-tangent value. The discrepancy from the independently solved
energies is numerical roundoff and is inspectable from the transition file.

No mixture-interior density rows are exported. For a requested logarithmic
grid, nodes strictly inside `(n_a,n_b)` are removed, both endpoints are inserted,
and additional one-sided interior anchors are inserted only when necessary to
give each phase at least three points. `rows=requested_points-removed+inserted`
is checked. A request for257 points gives253 rows:6 removed,2 inserted.
For sparse requests3,9,33 the same phase-wise contract is verified.
The native `iphase` values distinguish the two pure phases without assigning
fictitious A/Z to mixed matter.

The maximum same-phase first-law residual is `8.995e-9`, minimum true Hessian
eigenvalue is `2.9474` in the declared scaled coordinates, and the minimum
uniform-minus-crust energy is `0.13277MeV per baryon` on the257-point request.
The finite common-tangent gap lies below the unconvexified source energy at
independent interior probes. No conclusion about a real first-order neutron
drip in nature is attached to this completion of an old approximate surface
functional.

## Source table comparison and unresolved precision

`BBP1971_TABLE12_BENCHMARK.csv` records unmodified Tables1/2 comparison values
and native computed results. Its source is `bbp_original_benchmark.f90`.
For the first eight overlapping rows, pressure differences range from about
`-1.074%` to `+0.0767%`; k,x,A,Z reproduce the source scale without parameter
fitting. At `n_b=.0204`, native `(A,Z,P)=(294.670,64.8784,.0475345)` compared with
printed `(296,65,.0477)`. At `.0475`, native `(547.137,89.0527,.175766)` compared
with `(548,89,.176)`.

The last listed Table2 density `.0789fm^-3` gives native P `.401078`, whereas
the printed reference is `.384`, a **+4.4473%** discrepancy. Native mass density
is `1.33405e14g cm^-3`, compared with the printed `1.30e14`. This exceeds ordinary
last-digit rounding and is retained as an unresolved source-table discrepancy.
It is not hidden inside an acceptance tolerance, and no table coefficient,
density, or model parameter was altered to fit it. The implementable equations
are available; exact recovery of the authors' unrounded numerical execution is
not claimed. The source's two-digit k precision is accounted for when checking
the independent pure-neutron Table7 values.

An independent cross-table check identifies an internal density inconsistency.
From Table1 at the same printed rho, `A=990`, `k=1.21`, `k_g=1.01` and
`n_N=17.8e-6fm^-3`, the defining identity
`n_b=A*n_N+(1-A*n_N/n_i)*n_g` gives `.0769676427403fm^-3`. Taking all four
printed half-last-digit intervals gives the full corner range
`[.0759399063569,.0780056334808]`, which is disjoint from Table2's printed
`.0789 +/- .00005`. A possible `.0769`/`.0789` density typo is an inference;
the printed `.0789` value remains unchanged in every fixture.

An independent fixed-rho native calculation at `rho=1.30e14g cm^-3` instead
finds `n_b=.0768969977621913`, `P=.384202233121209`, `A=987.8898580211`,
`Z=119.643446962`, `x=.121110107559`, `k=1.20837281385`, and
`k_g=1.01008036442`. The pressure lies in the printed `.384 +/- .0005` box,
but A differs from990 by0.213%; this does not establish exact recovery of the
original numerical tables. The independent probe and its output are preserved
as `BBP1971_FIXED_RHO_BENCHMARK.f90` and `.txt`. The new generator therefore has
an explicit source-reproduction qualification, not an unexplained tolerance
change or a silently corrected source density.

## Verification and integration contract

The dedicated suite builds `compose_eos_types.f90`, `compose_eos_io.f90` and the new BBP module
with `gfortran -O2 -fcheck=all -ffpe-trap=invalid,zero,overflow`. No SciPy/native
SLSQP dependency is needed. An initial SLSQP reuse encountered its deliberate
`0/0` initialization of unused multipliers under these traps; the BBP component
instead uses its own small, deterministic Newton solve. This avoids global
floating-point-environment changes and keeps the computation lightweight.

Tests cover the source limits, direct pressure identity, energy derivative
oracles, the exact omitted-term identity, curvature and phase comparisons,
common tangent and one-sided first-law checks, raw failure preservation,
composition conservation, sparse grids, invalid/nonfinite/extreme inputs,
strict parameter-file records, and parallel byte determinism. In particular,
an overflowing decimal exponent is rejected before Fortran real conversion.
There is no `STOP` or `ERROR STOP` in the production module.

The independent integration review found and corrected a CompOSE convention
error: q5 is the conserved lepton-number chemical potential, so `q5=0` in this
neutrino-free beta-equilibrated model. q4 remains `-mu_e/(m_n*c^2)`. Every tested
generation now passes a real `write_compose` followed by `read_compose` with
explicit q5-zero and negative-q4 assertions. The gas-free-face Hessian gate and
the finite-search/global-claim qualifications were also added in response to
that review. Final dedicated generator tests: **27 passed in7.55s**.

Each generator returns `status='ok'` for successfully produced native data.
This is **not** the matching acceptance decision. Export full `eos.compo`:
`Yn=(1-u)n_g/n_b`, `Yh=u n_i/n_b`, `Ye=u n_ix/n_b`, `Yp=Ymu=0`, A and Z.
All presence flags and C3/C4 inventories are explicit. All three diagnostic
files must be integrity-sealed by the shared frontend:

- `bbp_original_diagnostics.json`;
- `bbp_original_states.csv`;
- `bbp_original_transitions.csv` (header-only for the two raw closures).

The component is stateless and reentrant. Independent model and pair tasks can
run concurrently. Within one density sweep, continuation is deterministic and
serial; its three-variable solve is too small to benefit from nested threads.
The matching batch, source admission, final publication and backups are owned
by the parent integration task and are not claimed complete here.

## R/C execution trace and handoff

| Stages | Concrete component evidence |
|---|---|
| R1–R3 | Fixed original-functional target, uploaded primary identity, verified equations and conventions, separate p251 closure |
| R4–R5 | Three explicit model constructions, reduced constrained energy, omitted-term formula, common-tangent construction |
| R6–R7 | Independent pressure and chemical-potential derivation, finite-difference gradient oracle, source tables and native numerical probes |
| R8–R9 | Published closure is executable but thermodynamically defective; raw variational branch nonconvex; derived Maxwell representation preserves one-sided phases |
| R10 | This source trace, strict parameters, Fortran benchmark and full test log |
| C1–C4 | Source-based targets, equation-to-routine mapping, bounded native solver, deterministic coarse-grained parallel execution contract |
| C5–C7 | New Fortran module, native diagnostics/composition, strict-trap tests, source and thermodynamic checks |
| C8 | Native derivative errors, curvature, endpoint residuals, sparse-grid behavior, parallel byte identity; per-model257-point generation takes about1–2s on the active remote workstation |
| C9 | Independent integration review corrected q5, gas-free-face curvature admission and scope overstatements; all9 printed anchors are preserved |
| C10 | Exact wrappers/IDs/data contract handed to parent; publication and all-pair survival results remain parent-owned |

No existing original-table model, pressure-integral model, source tolerance,
matching method, shared catalog or app entry point was changed by this component.
