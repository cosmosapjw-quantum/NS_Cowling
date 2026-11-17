# B2 native core support audit

Status: COMPLETE_BOUNDED_B2. All three actual native inputs contain 385 states on [0.02, 1.0] fm^-3. The higher Skyrme metadata limits are not generated support.

## Numerical evidence

| Model | species | continuous cubic P(mu) cs2/c2 | native probe cs2/c2 | max relative first-law residual |
|---|---|---|---|---|
| GM_GM1 | npe | [0.0006333208331554208, 0.756492775710528] | [0.0006333325785907743, 0.7564909824951508] | 2.032619e-11 |
| RG_SLY4 | npemu | [0.002977075640087629, 0.8318776264233295] | [0.002977067802630777, 0.8318619126432452] | 1.197242e-11 |
| RG_KDE0v1 | npemu | [0.0034266044293448502, 0.7490182824772836] | [0.0034266097528477036, 0.7490058188855996] | 1.786193e-11 |

Every one of the 384 intervals/model was evaluated at endpoints and every real stationary candidate of the relevant polynomial or rational derivative. This certifies the declared interpolant numerically; it does not prove the full underlying functional is stable at every real density. Native geometric midpoints and one-sided endpoint derivative stencils provide additional functional evidence. No densities outside [0.02,1.0] were sampled by the native probe.

## Thermodynamic construction

With x=(mu-mu0)/Delta_mu, P(x)=P0+Delta_mu*n0*x+[3*Delta_P-Delta_mu*(2*n0+n1)]*x^2+[-2*Delta_P+Delta_mu*(n0+n1)]*x^3. Define n=dP/dmu and epsilon=mu*n-P. Then d epsilon/dn=mu, dP/dn=n/(dn/dmu), and cs2/c2=n/(mu*dn/dmu). Native source states and energies are retained at knots to floating point roundoff.

## Interpretation and exclusions

- No TOV, maximum-mass, radial-stability or Cowling claim.
- No global inhomogeneous ground-state or crust-transition-density determination.
- No spin, gradient, finite-wavelength or pairing stability claim.
- No source-provider equivalence claim; unresolved GM1 calibration/provider residual is preserved.
- No support beyond actual n=1.0 fm^-3 certified; metadata upper bounds are computational capability only.
- Positive beta-path compressibility is not full microscopic stability.
- Native interior derivative evidence and polynomial certificate are distinct from a rigorous continuum theorem for the functional.
- GM1 is intentionally npe; excluded muons are not added silently.
- The native frozen diagnostic fixes proton fraction while leptons equilibrate at combined susceptibility; it is not a fully frozen-species buoyancy closure.

The upper endpoint 1.0 fm^-3 is an actual-input boundary, not a measured causal crossing. Subsaturation homogeneous states can be metastable relative to nonuniform matter; a beta-path mechanical check does not determine the crust transition. Skyrme muon onset is a continuous lepton threshold, not a first-order nuclear transition. GM1 retains its intentional npe species convention.

Full exact hashes, residuals, onset locations, native probe identities and A-anchor differences are in REPORT.json. Reproduce with:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python scripts/audit_core_support_20261117.py
