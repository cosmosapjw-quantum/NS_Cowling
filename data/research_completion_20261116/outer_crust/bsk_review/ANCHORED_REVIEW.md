# Independent review: boundary-anchored BSk22 analytic variant

Status: PASS_SCOPED_BOUNDARY_ANCHORED_ANALYTIC_VARIANT. This is a distinct,
first-law-consistent barotrope retaining the published pressure as a function of
baryon density and the separately fitted composition. It is not a new microscopic
ETF/Strutinsky stationary solution and does not independently determine neutron
drip from that composition.

The construction is mathematically consistent. Solve P_source(n0)=P_outer and
set e0=mu_outer-P_outer/n0. Then e(n)=e0+integral(P_source(u)/u^2 du) satisfies
P=n^2 de/dn and mu=e+P/n. Equal pressure and chemical potential are satisfied at
the anchor, with positive density jump. Since dmu/dP=1/n, n0>n_outer gives
the inner branch the smaller chemical potential immediately above that pressure.
A numerical comparison at relative pressure offsets 1e-5, 1e-4 and 1e-3 confirms
this ordering against the formal metastable continuation of the outer branch.
This is local ordering of two declared approximations, not global microscopic
phase competition.

An independent 60-digit implementation of the published pressure formula, its
root and its quadrature agrees with the module: anchor density absolute error
1.20e-18 fm^-3, energy error at most1.14e-13 MeV/baryon, pressure error at most
2.01e-14 relative, and sound-speed/chi error at most1.81e-14 relative. Nine
independent quadrature cases include both endpoints.

A 40-digit interval calculation over512 contiguous source-density intervals
bounds dP/dn between0.51959 and1.91432 MeV. Positive P' implies g(n)>=g0.
Also dg=dP/n<=dP/n0 implies g(n)<=g0+(Pmax-P0)/n0. These bounds enclose
c_s^2/c^2 between0.00054916 and0.00203745 throughout the reconstructed interval.
This enclosure concerns the explicit analytic formulas and reconstruction, not
calibration of the microscopic BSk22 functional.

The source composition remains a finite-error fit. At the anchor its free-neutron
fraction is0.00491945 and its independently fitted neutron chemical potential
exceeds the reconstructed baryon chemical potential by0.0237608 MeV. Across1025
density points the latter residual ranges from-0.021495 to+0.025630 MeV. The
energy per baryon changes from the original energy fit by+0.006174 to+0.102600 MeV.
These differences are preserved as diagnostics. The construction resolves the
barotropic boundary/first-law mismatch; it does not resolve microscopic chemical
equilibrium of the composition.

A metadata defect was found and repaired by the owner. Copying the source state
dictionary initially retained the old dlogP/dlogrho and the old energy-derivative
sound-speed field after epsilon had changed. Current fields now use the
reconstructed energy; the old source fields have explicit source labels. No
reviewer edits were made to the implementation.

Nonfinite and out-of-range state densities, nonfinite anchor inputs, negative
pressure, inconsistent Euler energy, invalid quadrature tolerance, and a missing
positive-density pressure bracket are tested as refusals. The runner's native
table tests are separate from this review's independent quadrature and interval
checks.

Reproduce:

```sh
python3 data/research_completion_20261116/outer_crust/bsk_review/review_bsk22_anchor.py
```

ANCHORED_REVIEW.json records reviewed source hashes, high-precision cases, interval
bounds, explicit anchor inputs, chemical residuals, local ordering and refusal
tests. anchored_review.log retains the execution output. REVIEW.json and
REVIEW.md separately audit the unchanged original analytic representations.
