# EOS matching adversarial repair — 2026-11-23

Baseline: e51dadd132125f0c5389da8c92c50b84eeb524e8.
Parent owns full source census, CLI/API exposure, integration and publication.
This worker edited only compose_eos_matching.f90 and three dedicated regression files.

## Confirmed defects and fixes

1. A first heuristic window was returned before its endpoint derivative or bridge
was certified. Any later failure immediately rejected the entire pair/method.
Minimal real-data reproduction: normalized GPPVA_FSU2H -> GM_GM1, label the last
crust row as a new one-row phase. Full input rejects at the endpoint derivative;
removing only that unretained last row admits i1=298, i2=79, overall consistency PASS.
Fix: native candidate iteration, cached one-sided derivatives, actual bridge
certification before selection, and continuation after each rejected construction.

2. Last149 crust rows, j>=12, n_core>=.045 (or .03 dissolved), 10% n/P margins,
.2 MeV F/mu margins and normalized-mean fractions (.12,.88) excluded legitimate
windows. Four native counterexamples failed before the repair and passed afterward.
Fix: preserve these rules only as preferred ordering; then examine every remaining
source endpoint pair with at least three retained points on its derivative side.
Necessary n/P/F/mu ordering and causal integral bounds replace arbitrary exclusions.
For a causal monotone n(mu), d(n/mu)/dmu>=0, so its pressure-average density obeys
n1*(mu1+mu2)/(2mu1) <= Delta P/Delta mu <= n2*(mu1+mu2)/(2mu2).
No interpolation uniqueness or continuum window completeness is claimed.

3. The first-law constrained cubic P(log n) bridge had only a401-point causal check.
A real native counterexample passed with max sampled cs2=.99999491972783217 but
at u=.1021 had cs2=1.0000050816809485. The old run exits1 in RED_FIRSTLAW.log.
Fix: continuous piecewise interval check. With local t in[0,1], span s and cubic P,
A=P+P_t/s+P_tt/s^2+P_ttt/s^3 and epsilon=(epsilon0+A0)exp(s*t)-A.
The causal margin epsilon+P-P_t/s is exponential minus quadratic; its second
derivative has at most one zero. Native code partitions there, isolates every
stationary point and evaluates the margin at endpoints and all stationary points.
P and P_t positivity use existing polynomial-extrema code. Positive integral
moments avoid cancellation in the exponential primitive. Physics reviewer separately
verified the analytic expression with Wolfram. Tests reject the counterexample and
admit its shifted-energy control (maxcs2=.81132482410733664).

4. A causal bridge could retain an acausal core tail. A Gamma3 synthetic source
reproduced native status ok before the guard (RED_TAIL.log, exit1). Retained prefix
and suffix must now have positive finite energy, nondecreasing pressure and Euler
chemical potential, increasing energy, and causal adjacent secants. DeltaP/Deltaepsilon
uses the preexisting verifier C10 tolerance1e-6. These necessary table constraints
do not certify a source interpolant or replace full assembled output verification.

5. match_report.method was32 characters, truncating the43-character first-law method.
It is now64; spline_mode is96. The exact method identity is tested. Parent exposes
this previously hidden native method in the frontend and serializes search counters.

## Verification

Dedicated pytest:4 passed in6.34s, native checked/trap and O2 executions.
Native search tests: rejected endpoint retry, beyond-last149 coverage, sub.2MeV
span, below.045 core, exact method identity, exhaustive failure counts, negative
retained energy, acausal retained core.
Existing native adversarial suites cover normalization, decreasing/constant Q,
first-law derivative, interior acausality, slice conventions and polynomial extrema.
See legacy_controls_final/record.json and MATCHING_FINAL_RECEIPT.json for final hashes.

## Runtime and scope

No saved mutable kernel state. Endpoint derivatives are cached per native call.
Gauss nodes are computed once per first-law construction and reused across output
queries, rather than rebuilt per point. Candidate search uses interval certificates
and endpoint/midpoint checks; the existing401-point construction check remains in
final assembly. Parent runner must separately require full matched-output consistency.
The exhaustive counters count finite supplied endpoint pairs; only the first certified
construction per requested pair/method is returned. Failed search is a finite native
construction outcome, not a universal physical nonexistence theorem.
