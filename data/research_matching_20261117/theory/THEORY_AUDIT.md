# Core–crust matching theory audit, 2026-11-17

Scope: one bounded independent R2–R9 theory/audit pass against checkout base
478c7b5de1514ac512f91156f36551478f36f688. Production repairs are owned by
matching_native_20261117 and matching_pipeline_20261117. This lane changes only
this report/evidence and new matching_theory_* tests. It does not claim a new
global no-go theorem or microscopic validation of a phenomenological bridge.

Skill context actually read: PhysMath Research Loop 1.3.0 main, POLICY_INTEGRATION,
SUBAGENT_MODE, THEORY_LOOP R2–R9 from
skill://plugins_6acdfbba921881919e0c17843706aef4/physmath-research-loop.
No recursive delegation, no Git mutation, one compute thread.

## 1. Physical variables and exact identities (derived)

T=0; one conserved baryon number after the selected composition constraints have
been imposed. Let n>0 [fm^-3], epsilon [MeV fm^-3] include rest energy,
F=epsilon/n [MeV], P [MeV fm^-3], mu=(epsilon+P)/n [MeV].
Here F is energy/free energy per baryon, not a free-energy density.
Let s=c_s^2/c^2=dP/d epsilon. CompOSE's numerical mass mn means m_n c^2 [MeV];
no c=1 assumption about a dimensional speed is needed. hbar and k_B do not enter
these cold barotropic differential identities.

Extensivity and the first law give

    d epsilon = mu dn,
    P = n mu - epsilon = n^2 dF/dn,
    dP = n dmu,     dmu/dP = 1/n.

Euler consistency at tabulated nodes alone does not imply the differential first
law between nodes. Assigning mu=(epsilon+P)/n can repair an algebraic residual,
but it cannot certify a separately interpolated P(n), epsilon(n).

Mechanical stability of a smooth branch requires dP/dn>0. For mu>0 the causal
contract is 0<s<=1. At a sharp first-order equilibrium interface, P and mu are
continuous, n jumps upward and

    Delta epsilon = mu_t Delta n > 0.

The equal-pressure density plateau is a pair of one-sided states. Values inside
the density gap are not additional homogeneous equilibrium states.

For a smooth P(mu) interpolation with retained endpoint densities,

    n_1 < Delta P/Delta mu < n_2

is a necessary strictly-stable matching condition. Causality further implies
g(mu)=n(mu)/mu is nondecreasing, hence

    n_2/mu_2 >= n_1/mu_1,
    (n_1/2mu_1)(mu_2^2-mu_1^2) <= Delta P
      <= (n_2/2mu_2)(mu_2^2-mu_1^2).

These are endpoint necessity tests, not complete certificates of the selected
interpolant. A window-search heuristic failing its extra margin requirements
does not establish that no thermodynamically admissible bridge exists.

## 2. What is actually exposed (implementation read)

python_frontend/compose_eos/api.py match_crust_core accepts five method strings:
css, cubic_F, quintic_F, cubic_P, quintic_P. Its historical “six-method” docstring
is stale. The separate construct_maxwell_transition API is available with
explicit local-charge-neutrality confirmation. method=maxwell in the common
matcher is deliberately unavailable. No direct-join or polytrope method is
exposed by this matching API.

| Method | Actual mathematical object | Physical meaning |
| --- | --- | --- |
| cubic_F / quintic_F | Polynomial or piecewise Hermite q(x)=P/n, x=ln n; F from its integral | Phenomenological smooth cold barotrope, with selected endpoint continuity |
| cubic_P / quintic_P | Polynomial or piecewise Hermite n(mu); P from its integral | Phenomenological smooth cold barotrope, with selected endpoint continuity |
| css | One constant-s bridge; continuous core join, crust-side upward density jump | An explicitly constructed additional phase/interface, not a root of the two original source P(mu) branches |
| construct_maxwell_transition | Equality root of two source P(mu) Hermite curves | Single-conserved-charge, locally neutral equilibrium of the supplied phases |

The native historical first_law_constrained_hermite_pressure_bridge is reachable
below the public Python whitelist but is not a sixth public frontend method.
A directly concatenated table requires either compatible continuous endpoint
states, an explicit admissible first-order interface, or disclosure of its
unresolved join. Monotone P and epsilon alone do not establish that contract.

A phenomenological smooth bridge need not maximize pressure among the two
original microscopic models at every mu: it explicitly replaces that interval
by a new effective barotrope. It must instead satisfy its own first law,
endpoint, stability, causality and source-retention contracts. Calling it
microscopic equilibrium would require further evidence.

## 3. Concrete defect: valid decreasing q is rejected

At baseline compose_eos_matching.f90 bridge_profile rejects dy_tot=y1-y0<=0
before distinguishing canonical q(x) from n(mu). Increasing n(mu) is physically
required. Increasing q(x)=P/n is not.

An explicit exact counterexample is

    x=ln[n/(0.04 fm^-3)],
    q(x)=(10-x) MeV,
    F(x)=(940+10x-x^2/2) MeV.

It satisfies P=nq, epsilon=nF, mu=F+q and

    dP/d epsilon = (q-1 MeV)/(F+q).

On n in [0.005,0.12] fm^-3, q,F,mu are positive and 0<s<1, even though dq/dx=-1
MeV. Direct differentiation proves both first-law identities. Wolfram Language
independently returned zero for d epsilon/dn-mu and n^2 dF/dn-P and returned the
same sound-speed ratio and dq/dn=-1/n (with the implicit 1 MeV coefficient).

Native reproduction uses 41 geometrically spaced crust rows from 0.005 to 0.04
and 61 linearly spaced core rows from 0.045 to 0.12. Both are samples of this same
physical EOS, with M=939 MeV normalization and exact Euler Q3.

Observed baseline public-API result:

| Method | Result |
| --- | --- |
| cubic_F | rejected: bridge endpoints must be finite and strictly increasing |
| quintic_F | same rejection |
| cubic_P | success, selected source indices (40,0), s range 0.00935349–0.00948951 |
| quintic_P | success, same indices, s range 0.00933935–0.00947370 |
| css | rejects because its required positive crust-side jump is absent |

The CSS rejection is appropriate for its declared topology. The F-method
rejection is a construction restriction incorrectly coupled to an unphysical
monotonic-q requirement. It is a false rejection, not an accepted acausal EOS.

Smallest repair: allow the canonical polynomial candidate for decreasing q,
then certify the actual physical inequalities over each whole interval.
Do not simply remove all shape gates from the fallback splines. A candidate
outside the implemented family should retain an explicit unsupported-construction
outcome, distinct from physical impossibility.

Regression: fortran_eos/tests/matching_theory_soft_barotrope.py.
Baseline command:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q
fortran_eos/tests/matching_theory_soft_barotrope.py.
Observed: 2 failed, 1 passed in 0.26 s; preserved in decreasing_q_red.log.
Production correction and green run are owned by the native/integration lanes.

## 4. Whole-interval repair contracts (derived; independently CAS checked)

For an F-method segment write x=x_1+L u, 0<=u<=1, L>0, q=q(u) and

    F(u)=F_1+L integral_0^u q(v) dv,
    mu=F+q,
    A=q+q_u/L,
    s=A/mu.

Sufficient directly checkable segment conditions are q>0, F>0, A>0 and

    mu-A = F-q_u/L >= 0.

All are polynomials when q is polynomial. Their extrema are found by evaluating
endpoints plus real derivative roots in the segment. This proves positivity of
P, epsilon, mu and their required monotonic derivatives, and 0<s<=1. It permits
q_u<0 when q+q_u/L remains positive. Fixed 401-point sampling cannot replace
this certificate.

For a P-method write mu=mu_1+D u, D>0, n=n(u),
P=P_1+D integral n du and epsilon=mu n-P. Sufficient conditions are n>0,
n_u>0 and

    H=mu n_u-D n >= 0,     s=D n/(mu n_u).

The endpoints additionally supply positive P and epsilon; their derivatives
are then positive. H_u=mu n_uu, an identity independently confirmed in Wolfram.
A piecewise interpolant needs these checks on every actual spline segment,
including short knot intervals; the global diagnostic grid can miss them.

Endpoint derivatives must use the same physical Euler mu as the bridge and
the output. A provider Q3 residual accepted within the declared 5e-4 tolerance
must not silently become an inconsistent core endpoint in just some methods.

## 5. Maxwell topology controls and limits

The actual native Maxwell kernel derives source mu from Euler, uses cubic
Hermite P(mu) with prescribed slopes n, checks complete-interval stability and
causality, and isolates all stationary points of the cubic pressure difference.
These are materially stronger than scanning a fixed dense grid.

If D(mu)=P_low-P_high has exactly one isolated zero in the common support and
Delta n=n_high-n_low>0 there, then D'(mu_t)=-Delta n<0. Continuity plus the
absence of another zero implies P_low>P_high to the left and P_high>P_low to the
right wherever both branches are defined. Thus the selected local topology
already has the correct two-phase pressure envelope on that common support.
No unsupported comparison outside either branch's data domain follows.

Independent exact-cubic native controls:

    x=mu-990 MeV,
    P_low=8+0.12x+0.001x^2 [MeV fm^-3],
    P_high=P_low+10^-6(x+20)x(x-20).

The units of each coefficient supply the stated pressure dimension. Source mu
knots [950,960,1030,1040] MeV place all three roots 970,990,1010 in a single
shared interval [960,1030]. Both phases remain positive, stable and causal.
Native output correctly rejects MULTIPLE_COEXISTENCE_CROSSINGS. Replacing the
perturbation by -0.01x gives one root with reversed density ordering and correctly
rejects NONPOSITIVE_DENSITY_JUMP. Disjoint source mu domains correctly reject
NO_COMMON_MU_SUPPORT.

fortran_eos/tests/matching_theory_maxwell_topology.py: 3 passed in 0.17 s;
preserved maxwell_topology.log. A support gap is an input-domain limitation.
It is not evidence against a crossing in unprovided source data.

## 6. Reference mass versus physical energy offsets

With M=m_n c^2, physical F=M(1+Q6), epsilon=nM(1+Q7).
Changing only the encoding reference to M' requires

    1+Q6'=(M/M')(1+Q6),
    1+Q7'=(M/M')(1+Q7),

and the corresponding chemical-potential re-encoding. This leaves physical
P, epsilon and mu unchanged. For Il=0 at fixed Yq, the one-dimensional Euler
chemical potential is mu_E=mu_B+Yq mu_q; source Q3 alone is not mu_E.

By contrast F -> F+C changes epsilon -> epsilon+nC and mu -> mu+C while leaving
P unchanged. It preserves the differential first law but changes the physical
energy density that gravitates and can change relative phase coexistence.
It is a model adjustment, not a normalization gauge. Fortin et al.'s explicit
core energy correction belongs to this latter category. No such physical shift
is authorized merely to make a matching root appear.

## 7. Primary literature and provenance

1. Fortin et al., Phys. Rev. C 94, 035804 (2016),
   https://arxiv.org/abs/1604.01944 ; original manuscript read via Exa.
   Section III.2, Eqs. (6)–(11): first-law integration and the physical changes
   caused by endpoint chemical-potential/density corrections.
   Section III.3, Eqs. (12)–(17): continuity, convexity and causal matching bounds.
2. Zdunik, Fortin and Haensel, A&A 599, A119 (2017),
   https://www.aanda.org/articles/aa/full_html/2017/03/aa29975-16/aa29975-16.html .
   Sections 1–2, Eqs. (6)–(7): cold equilibrium mu and enthalpy integral relation;
   interpreted with epsilon=rho c^2 here to avoid the source's notation ambiguity.
3. Canullan-Pascual et al., arXiv:2502.02373,
   https://arxiv.org/html/2502.02373 .
   Section 2, Eqs. (3)–(10): constant-s P(mu) interpolation and endpoint density
   gaps. The extracted Eq. (1) has a dimensional typo; this audit derives
   epsilon=mu n-P directly rather than copying that expression.
4. CompOSE Reference Manual, arXiv:2203.03209,
   https://arxiv.org/abs/2203.03209 ; retrieved through Parallel primary-source
   excerpts. Thermodynamic definitions include rest energies and define Q1–Q7.
   Quick Guides, https://arxiv.org/html/2311.04715v1, confirms table encoding.
5. SciSpace independently recovered the 2016 journal record:
   https://scispace.com/papers/neutron-star-radii-and-crusts-uncertainties-and-unified-3ln9628sxi .
   Used for bibliographic cross-check, not as authority for equations.

Evidence classes: formulas in §§1,3–6 derived; referenced physical contracts
literature-supported; decreasing-q rejection and three topology outcomes
implementation-verified on the baseline native executable. Whole-interval
smooth-bridge repair remains a native-lane implementation obligation.

## 8. Integration readback

The native lane subsequently implemented polynomial-extrema checks and enabled
decreasing-q canonical candidates. This lane read the actual updated
certify_profile, certify_spline and recursive polynomial_roots_unit code.
The local segment coefficients, antiderivative base accumulation, canonical
scale factors and shifted mu multiplication agree with §4. No algebraic defect
was found in that narrow integration review. This is floating-point polynomial
extrema checking, not rigorous interval arithmetic. Fresh executable results
must be recorded after the coordinated full build.

Identity qualification: the checkout base SHA fixes the starting project
identity. Production agents were changing separate owned paths concurrently;
the initial auto-built executable was not accompanied by a build-input manifest.
Thus “baseline” above means the observed pre-correction behavior of the
increasing-q gate, not a claim of byte identity of every build input with the
base tree. EVIDENCE.json records explicit artifact/executable readback digests.

A degenerate limit of the same counterexample family is q=10 MeV and
F=940 MeV+(10 MeV) ln[n/(0.04 fm^-3)]. It is cold, stable and causal with
s=q/(F+q). A pre-repair CLI probe on the same grids rejected cubic_F,
quintic_F and even cubic_P as “no convex matching window found”.
find_window divides by q_2-q_1 in its heuristic thf, so this limit requires
an explicit zero-denominator guard in a floating-point-trap build.
This is an observed construction-domain restriction, not physical nonconvexity.
A truthful unsupported-flat-profile outcome is a sufficient minimal repair;
the current audit does not require extending the interpolation family.

## 9. Final targeted verification

After native owner freeze, this independent lane ran its six tests against
build/matching_native/api/compose_eos_fortran via COMPOSE_EOS_FORTRAN_EXE.
Result: **6 passed in 0.26 s** (theory_six_green.log).
Executable SHA256:
920a4512e56d028f3f7e3e03680d5e32b915d297b3936963c94467b9ebde716f.
Matching source SHA256:
72ff79e1709f51fa479996fd051a3735ca539fea64e213c28fd3fabc884391f2.
EVIDENCE.json binds the command, readback hashes and observed result.
The decreasing-q defect is closed for the tested cubic_F/quintic_F cases;
the three Maxwell topology controls remain green.
The native owner also reports trap-tested exact-flat canonical support; that
flat-case native-driver result belongs to its receipt, not to these six tests.

CSS root uniqueness follows directly by writing its pressure gain as
Delta P(s)=integral_[mu1,mu2] n2 (mu/mu2)^(1/s) dmu.
For 0<s<=1 and mu1<mu2, every interior integrand increases strictly with s,
so the gain is strictly increasing and the scalar root is unique when
bracketed. The separate n_bridge>n_crust check enforces the declared positive
density-jump topology. The Seidov diagnostic concerns onset stability of a
small new-phase core; it does not certify an entire stellar sequence.

This closes this lane's bounded primary scientific pass and the narrow review
of the repaired solver mathematics. Production corpus integration remains
owned by the parent, without a further theory-review cycle.
