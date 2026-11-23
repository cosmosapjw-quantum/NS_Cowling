# Stellar solutions with reconstructed BPS, a bridge, genuine NV rows and a core

Each input is a **reconstructed version** of the outer crust, with historical mass and lattice conventions documented in the preceding release. The 110 public matches span 22 core source keys and five matching methods. Candidate IDs distinguish different source representations of the same named EOS. The tests here use cold, beta-equilibrated, spherical perfect-fluid stars in general relativity. They evaluate the equilibrium quadrupolar tide. A dynamical Cowling f-mode on these eleven-interface stars requires perturbation and conversion conditions at every interface and a verified mode identification. Those eigenfrequencies remain unresolved.

## Native solver

The solver reads the original CompOSE directory, including its phase labels. Pressures and energy densities are in MeV/fm³, and interior integration uses geometrized units with lengths in km. Interpolation is linear in log energy versus log pressure inside each continuous table interval. The active portion of each interval must satisfy positive compressibility and the existing sound-speed ceiling of 1.00001 in units of c². This interpolation supplies the barotrope used by the stellar equations. The independent endpoint test in `../causal_endpoint` also uses the tabulated baryon density.

Write x = ln(Pc/P), u = r², v = m/r and y = rH′/H. With g = v + 4πuP and A = P(1−2v)/[(ε+P)g], du/dx = 2uA and dv/dx = (4πuε−v)A. The code advances the corresponding Hinderer equation for y using the analytic derivative of the current log-linear cell. The regular central expansion starts at x ≤ 10⁻⁸ with m/r = (4π/3)εc u to leading order and y(0) = 2. Classical RK4 clips every accepted step to the current cell endpoint. The reported maximum step is in x.

Two bounding pressures with different phase labels may coincide within 32 units of binary64 relative machine epsilon to accommodate the CompOSE P/n round trip. A positive energy-density jump is required there. The mass and radius continue across each interface. For an outward crossing the code adds Δy = −4πr³(εinside−εoutside)/(m+4πr³Pt). Widely separated changes of phase label remain ordinary positive-width cells. An unlabelled exact plateau is rejected. At exactly the transition pressure, a vanishing-radius inner core is handled by starting in the outer phase. Sub-step comparisons with smooth and one-jump radial DVODE solutions, central interfaces and consecutive jumps are in `regression.txt`. [Takátsy and Kovács, 2020, Eq. 11](https://doi.org/10.1103/PhysRevD.102.028501) · [Hinderer, 2008, corrected manuscript](https://arxiv.org/abs/0711.2420) · [Hinderer, 2009, erratum](https://doi.org/10.1088/0004-637X/697/1/964)

The grid ends at Pcut = 5.628660207790336 × 10⁻¹⁵ MeV/fm³ with εcut = 5.609588603804460 × 10⁻⁹ MeV/fm³. The reported R, k₂ and Λ use this finite endpoint and should be interpreted as cutoff proxies. Extending them to a physical P = 0 surface needs an outer boundary model. The output separately gives the magnitude 4πR³εcut/(m+4πR³Pcut) of a hypothetical local energy-density drop by εcut. This quantity measures one local correction to y. It does not bound the radius contribution of an arbitrary unmodelled exterior. M is the enclosed gravitational mass at the same cutoff.

## Public outcomes

Four requests per candidate were run at maximum log-pressure steps of 0.04, 0.02 and 0.01. All converged stars cross 11 labelled discontinuities. For the mass target, the native solver samples 32 geometric intervals from 5 to min(300,Ptable,max) MeV/fm³. It bisects the first **sampled upward mass crossing** to an absolute tolerance of 10⁻⁸ solar masses. This survey can miss a narrow branch between consecutive samples. Radial mode stability and the maximum mass require separate calculations.

| Request | Public candidates | Converged | No failing endpoint interval anywhere in table | Only above central pressure | Star reaches beyond first failing interval |
|---|---:|---:|---:|---:|---:|
| Pc = 25 MeV/fm³ | 110 | 105 | 80 | 20 | 5 |
| Pc = 50 MeV/fm³ | 110 | 105 | 80 | 15 | 10 |
| Pc = 100 MeV/fm³ | 110 | 100 | 75 | 15 | 10 |
| M = 1.4 solar masses | 110 | 100 | 75 | 15 | 10 |

At 1.4 solar masses each of five methods on core_002 (APR_APR) and core_003 (BBB_BHF2) lacks a sampled crossing. Other input representations named BBB_BHF2 do reach the target. Five methods each on core_028 (XMLSLZ_DDME-X) and core_031 (XMLSLZ_TM1) traverse a source interval failing the cold, causal endpoint bound. Their stars are numerical diagnostics. For core_001 (2020-MVCD-D1MS-full), core_029 (XMLSLZ_NL3) and core_030 (XMLSLZ_PK1), the first offending bracket is above the computed central pressure.

In the 75 target-mass configurations on tables passing the independent endpoint test, Rcut ranges from 11.14211 to 14.57250 km and Λcut from 210.482 to 1245.861 at the finest mesh. Both lower extrema are from core_016 (BBB_BHF2), and the upper extrema from core_010 (NL3_PROVIDER_PRECISION). Across the five matching methods at fixed core, the largest radius spread in this subset is 0.0005233 km, and the largest Λ spread is 1.384, both for core_006 (GM_GM1). These spreads measure sensitivity to the five matching prescriptions at the supplied source resolution. The local surface-jump scale ranges from 6.25 × 10⁻¹¹ to 1.40 × 10⁻¹⁰ over the same 75 cases.

Between the two finest meshes, the maximum relative changes across 100 numerical 1.4-solar-mass solutions are 2.49 × 10⁻⁹ in Rcut and 8.24 × 10⁻⁸ in Λcut. Across all mesh levels the maximum absolute mass residual at the target is 9.92 × 10⁻⁹ solar masses. These tests hold the interpolation, EOS data and surface cutoff fixed. The full per-request output, including failed requests, is in `public_stars.csv`. The per-core method spread and the observed maximum step-to-step changes are in `matching_method_spread.csv` and `mesh_convergence.csv`.

## Replay from repository root

The 110 inputs are selected by `qualified`, `has_BPS`, `has_synthetic_bridge` and `has_genuine_NV` in COMPOSITE_RETENTION.csv. The additional protected inputs are excluded from the public run. The public source archive has SHA-256 `18a7dabd7bb52963daa82ed597d7e45bc5290060710774014a8287b039e3e898`. On a clean checkout, first unpack the previous release's `Crust_matching_evidence_20261124.zip` member directory named by `src` below. It contains both the nested public EOS archive and the retention index.

```sh
src=data/eos_reconstructed_matching_20261124/remaining_matches/BPS_MTW1964_MS1966_NV1973_Hermite
receipt=data/eos_native_upgrade_20261124
run=$(mktemp -d)
unzip -q "$src/PUBLIC_MATCHED_EOS.zip" -d "$run/tables"
bash fortran_eos/build.sh --ensure "$run/build"
bin="$run/build/compose_eos_fortran"
"$bin" tov-love layered "$run/tables/cell_000001" 50 .01
"$bin" tov-love layered-target "$run/tables/cell_000001" 1.4 5 300 .01
python scripts/run_layered_tov_campaign.py --binary "$bin" --tables "$run/tables" \
  --index "$src/COMPOSITE_RETENTION.csv" --causality "$receipt/causal_endpoint/FULL_CRUST_110.csv" \
  --workers 3 --output "$run/public_stars.csv"
gfortran -O2 -std=f2018 fortran_eos/tests/summarize_layered_mesh.f90 -o "$run/summary"
awk -F, '$10=="ok" {print $1,$5,$6,$11,$12,$15,$17}' "$run/public_stars.csv" | \
  LC_ALL=C sort -k1,1 -k2,2 -k3,3gr | "$run/summary"
awk -F, '$10=="ok" && $6==".01" && $5=="M1.4" {print $2,$9,$12,$15}' "$run/public_stars.csv" | \
  LC_ALL=C sort -k1,1 | "$run/summary" method-spread
```

Both native entry points accept an optional final `--first-failure Pleft Pright` or `--audited-clean` that represents an independently audited table. With no annotation the native result records `unaudited`. A central pressure equal to Pleft only touches the lower boundary. One equal to Pright has traversed the entire interval. The frontend dispatches independent native processes and writes their diagnostics. A discrepancy in byte-level replay across compiler builds should also be checked against the reported numerical residuals.
