# Density-resolved source diagnostics

The source-preserving constructors now attach a physical_domain record to each
core/crust receipt. Native whole-table C1–C11 diagnostics remain available for
comparison. The new record retains every source row and reports pressure, total
energy density, baryon chemical potential, Euler applicability/residual and
composition masks. Adjacent-row intervals report two dimensionless estimates
of c_s²/c², dP/dε and (dP/dn)/[(ε+P)/n] using midpoint specific enthalpy.

Mechanical/causal findings, Euler consistency and composition availability are
distinct reported quantities. Phase changes and coarse sampling split connected
interpolation ranges. COARSE_SAMPLING only identifies spacing above the existing
log-step 0.25 threshold. It does not assert a missing physical density interval.
Its secant_diagnostic_pass and actual secant values remain available.

The first acausal APR source secant is bracketed by n=0.88 and 0.90 fm^-3, with
dP/dε=1.02523235. Source samples through 0.88 pass these diagnostics. The physical
crossing is unresolved within/near the bracket until the model or an accepted,
converged interpolant is evaluated. The corresponding KDE0v1 bracket is
[1.425315,1.442943] fm^-3. Both are density-local findings. Attribution to an
intrinsic model limitation, table error or convention mismatch needs independent
reproduction; the diagnostic does not assign that attribution.

The five core and thirteen crust source tables were built and read back natively.
The detailed receipts and full row/interval diagnostics are retained in
data/research_generator_20261116/source_domain_receipts.json.gz.
source_domain_summary.json provides ranges and input/source identities.
All 13 crust outputs contain 4172 rows in total.

TW uses pinned hashes of all five CompOSE source files. Its first 217 physical
and composition rows are compared against the archived DD2 outer segment before
extraction. The receipt keeps UNRESOLVED_BPS_BSk22 and the bibliographic conflict.
The source-faithful extraction is supported by the table identity/equality.

Tests cover acausal tails with usable lower-density samples, coarse-grid
diagnostics without physical exclusion, phase boundaries, missing lepton masks,
the Il=0 Euler charge term, nonfinite rows, immutable inputs, invalid density
grids, all 18 component readbacks and tampered TW source bytes.
The final targeted suite has 23 passing tests after independent review.

Independent review reproduced two diagnostic defects and repaired them. Cold
thermodynamic support now requires the per-row Q6=Q7 identity at the native
2e-6 tolerance. Baryon chemical-potential decreases are recorded independently
from compressibility, because charged/composition-changing paths require their
appropriate constrained chemical potential. Both counterexamples are retained
in domain_logs/review_red.log and the repaired run in review_green.log.

The monotonic-pressure comparison allows 64 times float64 machine epsilon
scaled by max(|P_left|, |P_right|, float64 tiny). Each interval records this
absolute roundoff threshold. This resolves the APR equal-pressure endpoint
difference of -2.398e-14 MeV fm^-3 with a threshold of 8.341e-14, while a
relative pressure decrease at pressures around 1e-20 remains flagged. The
physical sound-speed tolerance remains unchanged. APR jump recheck and the
roundoff regression RED/GREEN logs are retained alongside other evidence.
