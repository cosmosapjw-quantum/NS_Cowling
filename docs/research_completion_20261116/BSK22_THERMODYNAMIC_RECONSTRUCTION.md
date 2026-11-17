# Optional boundary-anchored BSk22 thermodynamic reconstruction

This separately named approximation preserves the published BSk22 pressure as a function of baryon density and its published composition fit, then reconstructs the total energy by the first law. It leaves bsk22_analytic.py and both original fit representations unchanged. The [author's EOS page](https://www.ioffe.ru/astro/NSG/BSk/nsesubr.html) recommends calculating one thermodynamic function from selected fitted functions to enforce consistency. The choice of the independently regenerated outer-crust endpoint as the integration anchor is an explicit construction here, not a new source fit or an ETF/Strutinsky calculation.

Let P_s(n) denote the original published pressure fit evaluated using the source's original rho(n) mapping. Determine n_0 from P_s(n_0)=P_drip, require n_0>n_outer, and set

e_0 = mu_drip − P_drip/n_0,
e(n) = e_0 + integral[n_0,n] P_s(u)/u² du,
epsilon(n)=n e(n),   mu_B(n)=e(n)+P_s(n)/n.

Here e includes the neutron rest energy. Integrating over x=log(n/n_0) gives integrand P_s(n_0 exp x)/(n_0 exp x), which is smooth on this segment. The code evaluates this with adaptive quadrature and checks it independently with96-point Gauss–Legendre quadrature. Pressure is preserved as a function of n; it is not reevaluated using the reconstructed rho in the old pressure-versus-rho fit.

The analytic pressure derivative is inherited from the source P_s(n), whereas sound speed and the logarithmic pressure derivative with respect to the **new** energy density are recalculated:

Gamma = n/P_s dP_s/dn,
cs²/c² = Gamma P_s/(epsilon+P_s),
dlogP/dlogrho_new = Gamma epsilon/(epsilon+P_s).

Independent review caught inherited derivative metadata that still referred to the original energy fit. It was repaired, and the source values now carry explicit source_ names. Tests compare these fields against finite differences of the reconstructed functions.

For the independently regenerated outer endpoint, the resulting inner endpoint is n_0=2.71440983347e-4fm^-3, compared with n_outer=2.69432006401e-4fm^-3: a positive0.745634% density jump. Pressure residual is3.36e-18MeVfm^-3 and the baryon chemical-potential residual is zero to reported precision. Since dmu/dP=1/n on both thermodynamically consistent branches, 1/n_inner−1/n_outer=−27.4695fm³ gives the correct local ordering: the inner branch has lower chemical potential immediately above the matching pressure. This establishes local ordering of two declared approximations, not a global microscopic minimum.

The analytic inner composition at n_0 has Ynf=.00491945. Therefore this boundary is not a newly solved neutron-drip condition. The free-neutron fraction need not vanish at the imposed connection of separate approximations. The reconstructed chemical potential differs from the independently fitted composition chemical potential by -0.02149 to+0.02556MeV across the sampled segment. That residual remains explicit and prevents a claim of an exact common variational solution.

The total energy per baryon shifts by+0.006174 to+0.102600MeV from the original published energy fit; the upper endpoint shift is+0.006174MeV. These changes are fully reported. No source coefficients are fitted to CompOSE or tuned to the upper-crust models.

Actual129-row and257-row tables were generated, emitted through the native CompOSE writer, and read back with composition preserved. Both native hard-check results were PASS. On common density points, the two direct analytic evaluations agreed exactly. First-law finite differences had maximum relative discrepancy3.77e-9, sound-speed differences2.58e-8, and the independent Gauss–Legendre/adaptive energy integrals differed by at most1.14e-13MeV. The combined original-fit and reconstruction tests passed23 cases.

API: BoundaryAnchoredBSk22 in python_frontend/compose_eos/bsk22_thermodynamic_reconstruction.py; explicit outer_density, outer_pressure, outer_mu_b and optional outer_epsilon arguments; .state(n), .eos_data(grid) and .anchor_diagnostics(). The supported interval is n_0 through0.002fm^-3. Reproducer: scripts/build_bsk22_reconstruction.py. Evidence: data/research_completion_20261116/bsk22_reconstruction/generated/RECEIPT.json and both native table directories.

This result admits an optional thermodynamically consistent, boundary-anchored analytic representation. It does not admit microscopic ETFSI regeneration, exact original energy-fit reproduction, a new neutron-drip solution, global phase completeness, or a core–crust join.


An independent review subsequently checked the reconstructed function against60-digit pressure quadrature: energy differences were at most1.14e-13MeV/baryon, with sound-speed and logarithmic-derivative relative errors at most1.81e-14. Its interval-enclosure calculation over the entire declared support bounded dP/dn between0.51959 and1.91432MeV and cs²/c² between0.00054916 and0.00203745. Local Gibbs ordering also passed at three positive pressure offsets. A denser1025-point sample found composition-fit chemical-potential residuals from -0.021495 to+0.025630MeV; this supersedes the smaller maximum of the257-point reporting sample. The independent verdict is PASS_SCOPED for this named approximation, retained in outer_crust/bsk_review/ANCHORED_REVIEW.json.
