# Native cold crust-core matching: current implementation contract

**Reviewed against source: 2026-11-23.** The historical filename is retained for existing links. This document describes the current Fortran kernels and Python frontend; historical six-token behavior is recorded only in section 10. It is not a survivor census or a claim that every source pair admits a physical match.

The implementation authorities are [compose_eos_matching.f90](../fortran_eos/src/compose_eos_matching.f90), [compose_eos_phase_transition.f90](../fortran_eos/src/compose_eos_phase_transition.f90), [compose_eos_verify.f90](../fortran_eos/src/compose_eos_verify.f90), the [native CLI](../fortran_eos/app/compose_eos_fortran.f90), and [api.py](../python_frontend/compose_eos/api.py). Thermodynamic calculations, root solves and admission checks execute in Fortran. Python supplies source selection, serialization, provenance and orchestration.

## 1. Seven methods and their dispatch

| Campaign method | Python entry point / native operation | Construction |
|---|---|---|
| `cubic_F` | `match_crust_core` / `match-core` | Continuous free-energy bridge, using a quadratic or piecewise cubic profile for `Q1` |
| `quintic_F` | `match_crust_core` / `match-core` | Continuous free-energy bridge with endpoint derivative constraints and a polynomial/spline construction |
| `cubic_P` | `match_crust_core` / `match-core` | Continuous pressure-potential bridge, using a quadratic or piecewise cubic profile for density as a function of chemical potential |
| `quintic_P` | `match_crust_core` / `match-core` | Continuous pressure-potential bridge with endpoint derivative constraints and a polynomial/spline construction |
| `first_law_constrained_hermite_pressure_bridge` | `match_crust_core` / `match-core` | Continuous near-pressure-plateau Hermite bridge; the first law fixes its integral |
| `css` | `match_crust_core` / `match-core` | One constant-sound-speed phase, a positive crust-side density jump, and a continuous join to the selected core endpoint |
| `maxwell` | `construct_maxwell_transition` / `maxwell-transition` | Genuine coexistence of two supplied single-phase barotropes at common pressure and chemical potential |

`match_crust_core(method="maxwell")` is rejected deliberately. A campaign must dispatch that method to `construct_maxwell_transition`; it must not substitute the continuous pressure bridge. The full first-law method name is a supported public token. A Gibbs construction is not among these seven methods: a one-dimensional cold barotrope does not supply the independent charge-chemical-potential surface required for a general two-charge Gibbs construction.

## 2. Variables, source conventions and input contract

The input slices are cold, one-dimensional barotropes. Densities are in `fm^-3`, pressure and total energy density in `MeV fm^-3`, and chemical potentials and energies per baryon in `MeV`. The CompOSE header values `mn` and `mp` are rest-mass energies, i.e. masses multiplied by `c^2`. Define the dimensionless sound-speed ratio

\[
 s=\frac{c_s^2}{c^2}=\frac{dP}{d\epsilon}.
\]

For the smooth/CSS kernel, the source energy per baryon and pressure are decoded using each source's own reference:

\[
 F=(1+Q_6)m_n,\qquad P=n_BQ_1,\qquad
 \epsilon_{\rm source}=n_B(1+Q_7)m_n.
\]

At exactly zero temperature, `Q6=Q7` gives `F=epsilon/n_B`. The constructor's allowed input discrepancy is stated below; it is not an instruction to add a rest-mass correction. The effective chemical potential used in the cold one-variable construction is

\[
 \mu_E=F+Q_1.
\]

The separately decoded source convention is

\[
 \mu_{\rm source}^{\rm eff}=
 \begin{cases}
 (1+Q_3)m_n,& I_l=1,\\
 (1+Q_3)m_n+Y_qQ_4m_n,& I_l=0.
 \end{cases}
\]

Thus, for `Il=0`, the stored `Q3` chemical potential excludes `Yq*mu_q`. The construction must undo that convention when comparing with Euler's relation, and restore it on output. Matching a fixed physical `Yq` slice does not permit relabeling it as `Yq=0`. The cold identities used by the bridge are

\[
 dF=\frac{P}{n_B}\,d\ln n_B,\qquad
 dP=n_B\,d\mu_E,\qquad
 \epsilon=\mu_En_B-P.
\]

The `match-core` input checks require at least three density points, finite arrays, positive strictly increasing densities, valid positive reference masses, `Il` in `{0,1}`, equal `Il` for both inputs, and equal physical charge-fraction values to `1e-12`. Integer axis indices alone do not identify a physical slice. The explicit cold-state checks are

\[
 |T|\le10^{-12}\ {\rm MeV},\quad
 \max|Q_2|\le10^{-10},\quad
 \max|Q_6-Q_7|\le2\times10^{-6}.
\]

Both inputs must carry beta-equilibrium metadata; for `Il=1`, the smooth/CSS kernel also enforces `max(abs(Q5)) <= 1e-8`. Every supplied source row is checked against the provider Euler relation before endpoint search, with relative tolerance `5e-4`. A metadata flag cannot override a failed chemical-potential check. Source-input rejection therefore remains distinct from an exhausted endpoint search.

The interface supports `core_compo_mode="strict_thermo_only"` and integer `n_bridge >= 8`. `literature_beta_eq` composition reconstruction is unavailable. If either input has an `eos.compo` file or usable abundances, `drop_composition=True` is required. The separate Maxwell kernel has its own input contract, including single-phase inputs, `Nadd=0`, explicit local-neutrality confirmation and its own source-interpolation checks; see section 6.

### Reference reencoding and Euler closure are separate operations

Changing the output reference from `m_old` to `m_ref` changes the dimensionless encoding, not the physical energy:

\[
 1+Q_{j,\rm ref}=\frac{m_{\rm old}}{m_{\rm ref}}(1+Q_{j,\rm old}),
 \quad j=3,6,7,
 \qquad
 Q_{j,\rm ref}=\frac{m_{\rm old}}{m_{\rm ref}}Q_{j,\rm old},
 \quad j=4,5.
\]

Density, pressure, physical chemical potentials and total energy density remain unchanged under that operation. Adding `mn_core-mn_crust` to `F`, or adding `n_B*(mn_core-mn_crust)` to `epsilon`, would change the physical EOS and is not performed. The smooth/CSS output uses the core's reference masses. `dm_mass_gauge_meV` and frontend `dm_crust` report the applied physical energy shift, **zero**; `normalization_reference_delta_mev` separately reports the header difference. See [the mass-normalization correction](MASS_NORMALIZATION_CORRECTION_20261115.md).

Output assembly also imposes the cold closure: it derives `Q3` from `F+Q1`, with the `Il=0` charge term restored, and sets `Q7=Q6`. This is separate from reference reencoding. It need not preserve a provider's small allowed `Q3` or `Q6-Q7` tabulation residual bitwise. `source_euler_relative_residual` reports the all-input Euler discrepancy, and `source_euler_residual_meV` reports the retained crust's absolute discrepancy. For source rows already satisfying the exact cold identities, the retained `n_B`, `P`, `F`, `epsilon` and physical charge chemical potential are preserved.

## 3. Endpoint search and its finite scope

The native search returns the first admitted construction in a deterministic order. Historical choices remain a **preference pass**, followed by a fallback over every remaining source endpoint pair. A failed endpoint derivative or failed bridge construction causes a retry of the next candidate; it does not immediately reject the source pair.

The preferred pass uses the historical core onset (`n_B >= 0.03 fm^-3` at `Abar=0` when composition is usable, otherwise `n_B >= 0.045 fm^-3`), upper-crust index bound `max(12, n_crust-148)`, endpoint margins `n2>1.1*n1`, `P2>1.1*P1`, `mu2>mu1+0.2 MeV`, `F2>F1+0.2 MeV`, and the historical normalized integral fractions in `(0.12,0.88)`. These preferences do not exclude a pair from the fallback.

Using zero-based indices, the finite endpoint domain is

\[
 i_1=N_{\rm crust}-1,\ldots,2,\qquad
 i_2=0,\ldots,N_{\rm core}-3.
\]

This leaves at least three source points on each endpoint's retained derivative side. The endpoint derivative uses up to six points, stops at a source phase-label boundary, and requires at least three same-phase points. Its one-sided PCHIP estimate is

\[
 s=\frac{P}{\epsilon}
 \frac{d\ln P/d\ln n_B}{d\ln\epsilon/d\ln n_B}.
\]

The current source-endpoint gate requires `0<s<1`; it does not clip an invalid derivative. These derivatives are cached within the matching call and are independent of other threads or processes.

After inexpensive retained-source checks, a candidate must satisfy `n2>n1`, `P2>P1`, `mu2>mu1` and `F2>F1`. Causality implies `d(n_B/mu_E)/dmu_E >= 0`, hence the necessary average-pressure bound

\[
 \frac{n_1(\mu_1+\mu_2)}{2\mu_1}
 \le \frac{P_2-P_1}{\mu_2-\mu_1}
 \le \frac{n_2(\mu_1+\mu_2)}{2\mu_2}.
\]

The cheap comparison allows only a machine-roundoff allowance of `64*epsilon` times the comparison scale. A surviving candidate is passed to the actual method constructor and its interval checks. The endpoint first-law residual must satisfy

\[
 \frac{|F_{\rm bridge}(1)-F_2|}{\max(1\ {\rm MeV},|F_2-F_1|)}\le5\times10^{-9}.
\]

The native receipt exposes `window_pairs_total`, `window_pairs_examined`, `window_candidates_constructed`, `window_candidates_rejected`, `window_preferred_candidates`, `window_fallback_used` and `window_search_exhausted`. The CLI includes search counters on matching rejection as well as success; frontend reports expose them under `window_search` when available. Inputs rejected before search do not acquire an exhaustion claim.

Exhaustion means that this implementation examined its finite source endpoint grid under the selected method and declared input domain. It does not prove that no different interpolation family, additional source data, different physically authorized support or interior endpoint choice could work. The code does not enumerate every possible bridge function or every successful endpoint pair. Internal polynomial/spline branch selection is also part of the method prescription, rather than a universal existence search.

## 4. The five continuous bridge constructions

Write `x=ln(n_B)`, `Delta x=ln(n2/n1)`, `Delta mu=mu2-mu1`, and `0<=u<=1`. All methods in this section preserve continuous density, pressure, energy and effective chemical potential at the joins. Their method names identify a construction family; they do not guarantee that a successful fallback is one global polynomial of that degree.

### Free-energy families: `cubic_F` and `quintic_F`

The profile is `q(u)=Q1=P/n_B`, with

\[
 n_B(u)=n_1e^{u\Delta x},\qquad
 F(u)=F_1+\Delta x\int_0^u q(v)\,dv,
\]
\[
 \int_0^1q(u)\,du=\frac{F_2-F_1}{\Delta x},\qquad
 P=n_Bq,\quad\mu_E=F+q,\quad
 s=\frac{q+q'/\Delta x}{F+q}.
\]

For the endpoint-derivative-constrained quintic family,

\[
 q'(0)=(\mu_1s_1-q_1)\Delta x,\qquad
 q'(1)=(\mu_2s_2-q_2)\Delta x.
\]

The cubic family first attempts a quadratic profile, corresponding to a cubic potential. Its prescribed fallback is a C1 cubic profile with one interior knot, whose integral is fixed by a scalar solve. The cubic family does not promise equality of the source and bridge sound speeds at both endpoints.

The quintic family first attempts a quartic profile obtained from a quintic potential, with exact supplied endpoint slopes. Its fallback uses a C1 piecewise cubic profile with five knots: the two endpoints and three interior points `ua, um, ub`. The historical receipt name is `adaptive_2knot_c2_hermite`; the actual five-knot array in `make_c2_knots` is authoritative. Integrating the C1 profile gives a C2 potential. The reported `spline_mode` distinguishes the polynomial and fallback cases.

### Pressure-potential families: `cubic_P` and `quintic_P`

Here `mu_E(u)=mu1+u*Delta mu` and the interpolated profile is `n_B(u)`:

\[
 P(u)=P_1+\Delta\mu\int_0^u n_B(v)\,dv,\qquad
 F=\mu_E-\frac{P}{n_B},\qquad
 s=\frac{n_B\Delta\mu}{\mu_E n_B'}.
\]

The required integral is `(P2-P1)/Delta mu`. For the quintic family, the endpoint slopes are

\[
 n_B'(0)=\frac{n_1\Delta\mu}{\mu_1s_1},\qquad
 n_B'(1)=\frac{n_2\Delta\mu}{\mu_2s_2}.
\]

These methods use the same quadratic/quartic and piecewise cubic profile machinery described above, now for density rather than `Q1`. Integrating the density profile enforces `dP=n_B dmu_E`. Smooth endpoint values alone do not establish mechanical stability or causality; the constructed profile must also pass the interval inequalities in section 7.

### `first_law_constrained_hermite_pressure_bridge`

This method uses a C1 cubic Hermite pressure spline in `u`, with internal knots `0.22` and `0.78` and internal pressures

\[
 P_a=P_t-0.015(P_2-P_1),\qquad
 P_b=P_t+0.015(P_2-P_1).
\]

The nominal endpoint pressure derivatives are `dP/du=n_B*mu_E*s*Delta x`; the pressure-spline shape limiter may reduce them. The scalar pressure parameter `P_t` is fixed by

\[
 F_2-F_1=\Delta x\int_0^1\frac{P(u)}{n_1e^{u\Delta x}}\,du.
\]

Then

\[
 F(u)=F_1+\Delta x\int_0^u\frac{P(v)}{n_1e^{v\Delta x}}\,dv,\quad
 \mu_E=F+\frac{P}{n_B},\quad
 s=\frac{P'(u)}{\Delta x\,n_B\mu_E}.
\]

The implementation uses the same 16-point Gauss-Legendre rule separately on each smooth pressure-spline piece for both the endpoint root and partial integrals. Quadrature nodes are cached within each construction. The former order-changing whole-interval quadrature is not the current algorithm. This is a continuous near-plateau bridge with zero density jump, `is_exact_maxwell=false` and `seidov_applicable=false`. Its full method name must not be shortened to `maxwell`.

## 5. CSS: one speed and a positive crust-side jump

Let the retained crust endpoint be `(n_c, P_t, mu_t)` and the selected core endpoint `(n2, P2, mu2)`. The constructor solves for one constant `s=(c_s/c)^2` and a bridge-side transition density `n_b`:

\[
 n_b=n_2\left(\frac{\mu_t}{\mu_2}\right)^{1/s},
\]
\[
 P_2-P_t=n_2\mu_t\frac{s}{1+s}
 \left[\frac{\mu_2}{\mu_t}-
 \left(\frac{\mu_t}{\mu_2}\right)^{1/s}\right].
\]

The current bisection bracket is `s in [1e-10,1]`. The solver rejects an unbracketed root, a nonpositive inferred density jump or a nonpositive energy jump. Across the crust/CSS interface,

\[
 P_-=P_+=P_t,\quad\mu_-^{E}=\mu_+^{E}=\mu_t,\quad
 \Delta n=n_b-n_c>0,\quad\Delta\epsilon=\mu_t\Delta n>0.
\]

Within the bridge,

\[
 n_B(\mu)=n_b\left(\frac{\mu}{\mu_t}\right)^{1/s},\qquad
 P(\mu)=P_t+n_b\mu_t\frac{s}{1+s}
 \left[\left(\frac{\mu}{\mu_t}\right)^{1+1/s}-1\right],
\]
\[
 \epsilon(\mu)=\mu n_B(\mu)-P(\mu).
\]

The bridge joins the core endpoint continuously in `n_B`, `P`, `epsilon` and `mu_E`. Its single sound speed need not match the source core derivative. There is no second internal CSS speed.

The output contains two adjacent transition rows with equal `P_t` and `mu_t`, but different `n_B` and `epsilon`. The CSS receipt sets `is_exact_maxwell=true` for this sharp crust-to-constructed-CSS interface. That field does not turn CSS into a coexistence solve between the two supplied source branches. The separate `maxwell-transition` operation has the latter meaning.

CSS also reports the Seidov small-core onset quantities

\[
 \Delta\epsilon_{\rm crit}=\tfrac12\epsilon_-+\tfrac32P_t,
 \qquad \mathcal S=\frac{\Delta\epsilon}{\Delta\epsilon_{\rm crit}},
\]

with `seidov_small_core_stable = (Delta epsilon < Delta epsilon_crit)`. This Boolean is a diagnostic, not a complete stellar stable-branch calculation or an unconditional matching rejection. The receipt separately identifies the transition's low/high indices, continuous core-join index, bridge speed and source endpoint speeds.

## 6. Genuine Maxwell coexistence

`construct_maxwell_transition` dispatches to the separate native phase-transition kernel. It requires caller-confirmed local charge neutrality and works with the supplied cold single-phase branches in `P(mu_E)`. Each branch uses cubic Hermite pressure interpolation whose knot slopes are the source densities, so `n_B=dP/dmu_E` and `epsilon=mu_E*n_B-P`.

The constructor seeks equal pressure at equal effective chemical potential with a positive high-minus-low density and energy jump. Source-branch interpolation, coexistence-root and phase-selection gates are separate from the `match-core` endpoint heuristic. Phase input validation requires at least three source knots, one constant `Iphase`, no retained abundance inventory and `Nadd=0`. A provider-specific auxiliary column may be omitted only by an explicit projection that preserves the immutable source, its auxiliary values and the primary thermodynamic columns; it must not be silently interpreted as a generic phase variable.

The Maxwell source interval certificate checks positivity, mechanical stability and causality for that particular cubic `P(mu_E)` representation. It is not a mandatory representation test for every `F(ln n_B)` bridge. A failure of that interpolant does not, by itself, prove that the underlying nuclear model admits no physical EOS.

See [the separate Maxwell interface document](FORTRAN_EOS_MAXWELL_TRANSITION_V1.md) for its dedicated API and output structure; the current code linked above remains authoritative for exact guards and status categories. No inference about a complete globally minimizing multi-phase equation of state follows from a two-branch finite-support coexistence result.

## 7. Validation: bridge, retained table and source coverage

### Continuous bridge admission

The polynomial and spline constructors check physical inequalities over their complete polynomial pieces. In the free-energy representation, with `x=ln n_B`, the checked quantities include

\[
 q>0,\quad F>0,\quad\mu_E=F+q>0,\quad
 D=q+q_x>0,\quad\mu_E-D\ge0.
\]

In the pressure-potential representation, they include positive pressure, density and energy, and

\[
 \frac{dn_B}{d\mu_E}>0,\qquad
 \mu_E\frac{dn_B}{d\mu_E}-n_B\ge0.
\]

Polynomial minima are evaluated at endpoints and all stationary roots found by recursive derivative-root isolation, rather than only at output nodes. For the first-law pressure spline, a local cubic `P(t)` with `s_x=Delta x*Delta u>0` obeys

\[
 \epsilon_t=s_x(\epsilon+P),\quad
 A=P+\frac{P_t}{s_x}+\frac{P_{tt}}{s_x^2}+\frac{P_{ttt}}{s_x^3},
 \quad\epsilon=(\epsilon_0+A(0))e^{s_xt}-A(t).
\]

Its causal margin `G=epsilon+P-P_t/s_x` is an exponential minus a quadratic. The implementation partitions the derivative by the possible single zero of `G''`, isolates stationary points of `G`, and checks them and the interval endpoints. It uses stable positive quadrature moments to evaluate the energy instead of subtracting large primitives. Pressure and its derivative receive their polynomial interval checks as well.

These are continuous-piece checks implemented in floating-point arithmetic with finite root/iteration tolerances, not formal interval-arithmetic proofs. The generated analytic bridge is additionally checked at 401 diagnostic points, including endpoints. The allowed bridge range is `0<s<=1`; the reported minimum/maximum sound speeds are sampled summaries, not necessarily exact global extrema. The finite samples do not replace the continuous-piece admission tests.

### Retained-source and assembled-table checks

Before constructing each candidate, cached crust-prefix/core-suffix tests check positive finite `P`, `F`, `epsilon` and Euler `mu`, nondecreasing pressure, increasing energy, nondecreasing `mu` within roundoff, and the adjacent causal secant bound

\[
 \Delta P\le(1+10^{-6})\Delta\epsilon.
\]

These are necessary finite-table conditions. They do not certify a chosen interpolation between every retained source knot. The supplied source domain is not repaired or truncated at failed physical intervals during this search.

The Python smooth/CSS API independently invokes the native consistency verifier on the assembled table. Its hard result includes C1-C6, C10-C11, pointwise physical-domain checks, resolved derivative-domain failures and thermodynamic-jump consistency. C7-C9 remain finite-difference diagnostics. The verifier records derivative masking and sample counts; phase boundaries and unresolved large log-density steps require attention to that coverage. Positive finite pressure/energy/effective-potential domains cannot be removed by derivative masking.

`match_crust_core` retains its tuple return even when that independent consistency result fails: `info["diagnostic_only"]` is then **true**. Such output is not an admitted matched EOS. A campaign must require the independent hard gate and its declared source-support/representation criteria before classifying a survivor. The genuine-Maxwell Python API instead raises when its assembled output fails the independent consistency verifier.

### Whole-selected support versus a finite component

Retaining an `Iphase` integer is categorical provenance, not authorization to interpolate across an inherited phase boundary. A source census must preserve the full selected source view and distinguish `WHOLE_SELECTED_RETAINED_CHAIN` from `FINITE_COMPONENT_ONLY`. A component candidate may be prescribed by source phase labels or an explicit authoritative source selection; deleting intervals because they failed a physical test cannot establish whole-source admission.

The [2026-11-23 source catalog](../data/eos_matching_reaudit_20261123/SOURCE_CATALOG.json) keeps 101 source records, including blocked sources. Its `matching_views` use source-declared contiguous phase runs and explicit auxiliary projections; legacy diagnostically selected good intervals are not matching authority. These are input/coverage records, not precomputed physical-survivor claims. The complete all-pair campaign's own receipt, source hashes and output checks determine its result.

## 8. Assembly, phase provenance and export

For a continuous bridge, the output rows are the retained crust prefix through `i1`, exactly `n_bridge` interior samples and the retained core suffix beginning at `i2`. Endpoints are not duplicated. CSS additionally inserts the explicit high-density transition row immediately after the crust-side transition row; the equal-pressure jump pair must not be deduplicated or smoothed away.

Abundance arrays are absent as physical data: their numeric placeholders are zero and presence/closure masks are false. The output nevertheless has `composition_file_present=true`, `composition_available=false` and a **phase-only `eos.compo` sidecar** (`Npairs=Nquad=0`). Retained source `Iphase` labels are copied verbatim. The bridge receives the smallest nonnegative label absent from both source label sets. Equality of labels from different providers does not assert equality of their physical phases.

The selected physical `T`, `Yq` and `Il` conventions are retained; the output explicitly represents a `1x1` selected thermodynamic slice. Matched output has `Nadd=0`. Original auxiliary quantities must remain accessible through the source/projection provenance, rather than being invented on the bridge.

The physical charge chemical potential `mu_q=mn*Q4` is retained on source portions. Across the bridge it is joined using `6u^5-15u^4+10u^3`; this is not a composition or charge-neutrality solve. Output has `Q2=Q5=0`, `Q7=Q6` and

\[
 Q_3=\frac{F+Q_1}{m_{n,\rm out}}-1
       -\mathbf 1_{I_l=0}Y_qQ_4.
\]

Export must preserve phase labels and any adjacent first-order transition rows. A successful CompOSE round trip does not establish that a downstream smooth TOV/Cowling consumer supports the discontinuity; that consumer's interface contract and tests are separate.

## 9. Interfaces and bounded regression evidence

For a smooth/CSS method, prepared and source-authorized inputs can be passed as follows:

```python
from compose_eos import match_crust_core, write_compose

matched, info = match_crust_core(
    crust, core,
    method="quintic_P", n_bridge=40,
    core_compo_mode="strict_thermo_only",
    drop_composition=True,
)
if info["diagnostic_only"]:
    raise RuntimeError(info["consistency"]["overall"])
# Apply the campaign's retained-source representation and coverage gates
# before reporting a survivor. This export alone is not that classification.
write_compose("output/matched_compose", matched)
```

A genuine-Maxwell call is separate. The local-neutrality argument is an explicit physical premise that must be supported for the supplied phases:

```python
from compose_eos import construct_maxwell_transition

matched, report = construct_maxwell_transition(
    phase_low, phase_high,
    assume_local_charge_neutrality=True,
    drop_composition=True,
    return_report=True,
)
```

Native CLI equivalents, from the repository root:

```sh
bash fortran_eos/build.sh build/native_matching
mkdir -p output/smooth_internal output/maxwell_internal
build/native_matching/compose_eos_fortran match-core \
  CRUST_INTERNAL_DIR CORE_INTERNAL_DIR output/smooth_internal quintic_P 40 true
build/native_matching/compose_eos_fortran maxwell-transition \
  LOW_PHASE_INTERNAL_DIR HIGH_PHASE_INTERNAL_DIR output/maxwell_internal true
```

Each internal directory contains `eos_data.csv` and `eos_data.nml`. `match-core`'s last Boolean means explicit composition-drop consent; `maxwell-transition`'s means local-neutrality confirmation and requires already prepared thermo-only inputs. They are different contracts.

The 2026-11-23 native matching regression source is pinned by SHA256 `4cc5f32c858a7d9f952f6f7c326c9465667b609ec02caaa264a828aa9ecf186b`. [test_matching_reaudit_20261123.py](../python_frontend/tests/test_matching_reaudit_20261123.py) builds two Fortran test programs under checked/FPE-trap and optimized configurations; the observed result was **4 passed**. Controls include endpoint retry, support outside the old last-149-row cutoff, sub-0.2-MeV endpoint differences, low-density core support, invalid retained energy, an acausal retained core tail, and an inter-node acausal pressure-bridge counterexample with a separate valid control.

The four suites driven by [run_matching_native_adversarial.py](../fortran_eos/tests/run_matching_native_adversarial.py) also passed at that source hash: thermodynamic adversarial cases, interval adversarial cases, slice conventions and extrema controls. Reproduction commands are:

```sh
PYTHONPATH=python_frontend OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 -m pytest -q python_frontend/tests/test_matching_reaudit_20261123.py
python3 fortran_eos/tests/run_matching_native_adversarial.py \
  --build-dir /tmp/ns_matching_controls_build \
  --log-dir /tmp/ns_matching_controls_logs
```

The complete frontend/pipeline CI also exposed 14 failures sharing the analytic Maxwell fixture in [test_maxwell_mass_gate.py](../python_frontend/tests/test_maxwell_mass_gate.py): 12 reference-reencoding cases and two explicit-composition-drop cases. The native coexistence solve succeeded, but the former six-knot grid left **zero** resolved derivative points in its eight-row output after correct gap/transition stencil masking. C10 and C11 therefore failed the existing coverage contract. This was not a change to the EOS curves or an energy-reencoding defect: the corrected mask also excludes one-sided three-knot endpoint stencils when they cross a masked interval.

The fixture now samples the same quadratic `P(mu)` and its exact derivative on 31 points from 950 to 1040 MeV. The unchanged transition at 990 MeV remains off-grid. Both `Il` controls have 31 resolved derivative points in the 33-row output and pass the consistency gate. The former six-point grid is preserved as two explicit negative-coverage tests: coexistence succeeds, `derivative_points=0`, C10/C11 fail, and the public API rejects that uncertified output. The focused mass-reference plus explicit-drop suites passed **25 tests in 0.62 s** with the pinned native executable. Only test-fixture sampling changed; no production tolerance, mask or physical equation was relaxed. Evidence is recorded in `/tmp/ns_eos_reaudit_20261123/matching/CI_MAXWELL_FIXTURE_REPAIR.json` and `CI_MAXWELL_FIXTURE_GREEN.log`.

An additional dependency check found a test-side identity error in [test_mass_normalization_contract.py](../python_frontend/tests/test_mass_normalization_contract.py): it matched output rows to both source phases by density alone. In the refined analytic fixture, the discarded low-phase state at `mu=995 MeV` and the retained high-phase transition at `mu=990 MeV` both have `n_B=0.13 fm^-3`. The regression now restricts comparison to the retained chemical-potential domain and the corresponding output branch, and explicitly retains that density-alias counterexample. The three affected test modules passed **45 tests in 2.31 s**; `CI_MAXWELL_FIXTURE_ALL_GREEN.log` records the run. No production change was needed.

The observed run receipts are `MATCHING_FINAL_RECEIPT.json`, `PYTEST_NATIVE_REGRESSIONS_FINAL.log` and `legacy_controls_final/record.json` under `/tmp/ns_eos_reaudit_20261123/matching/`; final campaign packaging must preserve those receipts. This evidence establishes the listed finite regression outcomes. It is not an all-source pass, proof of interpolation uniqueness, microscopic equilibrium calculation, stellar stable-branch result or blanket guarantee about arbitrary grids.

## 10. Historical records, retained for provenance

An earlier version of this document described a six-token interface that included `maxwell` as a continuous near-plateau alias, two-segment continuous CSS, an additive crust energy shift, zeroed phase labels, a last-approximately-150-row endpoint filter and sample-only bridge checks. Those descriptions are superseded by the current implementation above. Historical fixtures and recorded pass counts do not validate those behaviors as current physical contracts.

The previous document recorded branch `research/fortran-eos-gap-fill-v1` at commit `e2b7ef319b9bb73dc85023a757dc1449c8f22e3c`, with successful Actions runs [42378176837](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42378176837) and [42378180891](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42378180891). It reported 44 frontend/provider tests, 58 manifest entries and source packaging at that commit. The recorded [source artifact](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42378176837/artifacts/13576661259) had ZIP SHA256 `28793bbbfa0708a8b24d0cc2e128b0833e9eea16e3d30ea547172b172b62474b`. These links and values are preserved historical evidence; availability and CI status were not revalidated by this documentation edit, and the artifact cannot substitute for the current exact-source receipt.

The [2026-11-15 mass-normalization correction](MASS_NORMALIZATION_CORRECTION_20261115.md) explicitly identifies the physical error in the earlier additive energy shift. The older [CSS decision note](CSS_SINGLE_SPEED_CRUST_JUMP_CONTRACT_20261115.md) records the change to one speed and a crust-side jump, but its historical numerical example includes the subsequently rejected energy-shift convention and must not be used as a current numerical acceptance fixture. Use the current Fortran implementation, source hashes and independently observed regression/campaign receipts for current claims.
