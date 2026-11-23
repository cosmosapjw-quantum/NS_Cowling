# Supplementary interpolated BPS + NV crust

This cold scalar barotrope retains the reconstructed BPS outer crust and the Negele–Vautherin (1973) ground-state energy fit of Table 4. A cubic interpolation between them is an explicit assumption. Both source normalizations and endpoint slopes are held fixed. The source supports that enter this composite are 6.02914×10⁻¹² through 2.58691171607×10⁻⁴ fm⁻³ for BPS and 6.30839140553×10⁻⁴ through 7.89×10⁻² fm⁻³ for NV. The historical outer crust includes 11 isobaric density jumps. In the native phase field, 201–212 are the BPS nuclear layers, 213 marks interpolated samples, and 214 denotes the NV author fit. The last two segment boundaries are continuous in pressure and baryon chemical potential. Composition is given for BPS alone. The compulsory Q4/Q5 columns in the thermodynamic-only extension follow the source's scalar equilibrium convention. [BPS](https://doi.org/10.1086/151216) · [NV1973](https://doi.org/10.1016/0375-9474(73)90349-7)

## Conditions on a direct join

At common pressure the difference μ_NV−μ_BPS decreases from +0.313396 to +0.166270 MeV over P = 3.69281137×10⁻⁴ through 4.90875255×10⁻⁴ MeV fm⁻³. Its derivative is 1/n_NV−1/n_BPS < 0 throughout this common support, since even the smallest NV density exceeds the largest BPS density. Thus the unshifted EOSs have no equal-P, equal-μ Maxwell join.

A stable interpolation with fixed endpoints must also have convex ε(n). Let n₀ be the BPS drip density and let H(n) = ε₀−ε_NV(n)+μ_NV(n)(n−n₀). Convexity requires H > 0 for a bridge with positive compressibility. H′ = (n−n₀) μ′_NV > 0 on the interval between the first two NV Table 3 densities. For the published fit, μ′_NV = exp(f) [f′+f′²+f″]/n. An exact rational Bernstein bound for the bracket on 1 ≤ log(n/10⁻⁴) ≤ 1.4 gives a lower bound of 0.453154. Yet H(4.00×10⁻⁴ fm⁻³) = −2.55905×10⁻⁵ MeV fm⁻³. Every ending density in the first two-cell interval consequently fails this condition.

## Chosen cubic and checks

Fix the left endpoint to BPS neutron drip. In ascending density order, test every row in the native 257-row NV fit. For each candidate let h = n₁−n₀, t = (n−n₀)/h, s = (ε₁−ε₀)/h. Use ε(t) = ε₀+h[μ₀t+a t²+b t³], with a = 3(s−μ₀)−(μ₁−μ₀) and b = (μ₁−μ₀)−2(s−μ₀). The endpoint values of ε and μ = dε/dn fix all coefficients. P = nμ−ε, so both pressures agree with their sources. The source sound speeds need not equal the interpolant's one-sided endpoint speeds.

The second derivative (2a+6bt)/h is linear in t. Positive endpoint curvatures prove positive compressibility throughout, and n₁ max(ε″)/μ₀ < 1 suffices for causality. The first accepted native row is **38**. The retained cubic has 46.2426 ≤ ε″ ≤ 4559.14 MeV fm³ and 1.2720×10⁻⁵ ≤ cₛ²/c² ≤ 0.00306108 throughout the bridge. The interval increases density by a factor of 2.43858. It uses 63 of the 906 output rows and occupies 3.827% of the full logarithmic density span. It replaces the first 37 NV fit rows and spans the first three tabulated NV-cell densities. The algorithm selects the first admissible row on the specified 257-row grid. Other grids or interpolation families can change that endpoint.

`generated/join_candidates.csv` records every rejected grid candidate and the final admissible one. `bridge_verification.json` records phase support, immutable inputs, endpoint residuals and independent checks. The Fortran reader recovered all 906 native rows. Thermodynamics and the full original BPS composition sidecar prefix were independently compared to the two input files. The complete output passes the native thermodynamic consistency gate. Composition coverage remains partial.

## Reproduction

From the repository root, with an existing BPS reconstruction and the 257-point NV1973 generation in the indicated directories, run the following in a POSIX shell. The generator writes to a fresh existing directory. The independent Python verifier checks the copy distributed in the campaign.

```sh
d=$(mktemp -d)
gfortran -std=f2018 -O2 -ffree-line-length-none -fcheck=all -J"$d" -I"$d" \
  fortran_eos/src/compose_eos_types.f90 fortran_eos/src/compose_eos_io.f90 \
  fortran_eos/tools/bps_nv1973_bridge.f90 -o "$d/bridge"
mkdir "$d/generated"
"$d/bridge" \
  data/eos_reconstructed_matching_20261124/generated_sources/BPS_MTW1964_MS1966/generated \
  data/eos_original_papers_20261123/generated_sources/NV1973_author_fit/generated "$d/generated"
diff -rq "$d/generated" data/eos_reconstructed_matching_20261124/supplementary_bps_nv/generated
python data/eos_reconstructed_matching_20261124/supplementary_bps_nv/verify_bridge.py
```
