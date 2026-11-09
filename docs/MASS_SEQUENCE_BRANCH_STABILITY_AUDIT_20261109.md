# Mass-sequence topology and radial-stability audit

Date: 2026-11-09  
Scope: GPPVA_FSU2 × BBB_BHF2 six-method pilot, frozen v5 EOS-cell TOV backend  
Evidence state: CI artifact inspected; 24-point branch topology is an observed numerical diagnostic, not a radial-stability result. A 96-point refinement is the next exact-head experiment.

## 1. Why the first mass-branch gate is not yet a stability gate

The initial mass-sequence extractor deliberately retains only the first connected rising segment of M(Pc), and stops at the first resolved decline. This is a conservative turning-point proxy for target-mass searches. It is not a proof that every point in the retained segment is radially stable, nor does it prove that every later rising segment is unstable.

The radial stability criterion is the sign of the lowest eigenvalue of the relativistic radial-pulsation problem. With the convention exp(i omega t), a nodeless fundamental radial mode satisfies:
- omega_0^2 > 0: radially stable;
- omega_0^2 = 0: marginal;
- omega_0^2 < 0: radially unstable.

Bardeen–Thorne–Meltzer (BTM) turning-point counting is a useful sequence-level diagnostic under its assumptions, but the direct radial Sturm–Liouville spectrum is the decisive test. Alford, Harris & Sachdeva explicitly compare the two and show why mass-radius extrema alone can misclassify exotic low-compactness branches. Essick gives further examples of stability re-appearance in specially structured low-pressure EOS sequences. Therefore neither “all later rises are stable” nor “all later rises remain unstable” is a safe general rule.

## 2. Exact pilot evidence at the preceding 24-point resolution

Actions run: https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42610037961  
Tested HEAD: d1d01e6b114d97cd496640037aafaf7ff7be126e  
Artifact: 13668724637, SHA-256 fb62570658e44f9847d4a08f73174e36f14c091b6960b3b847206232a73b1534

The branch test was performed on six methods for the same GPPVA_FSU2 × BBB_BHF2 pair. The v5 TOV baseline M and R at the frozen reference Pc reproduced to near machine precision for all six methods. The helper tests and both mock-response test groups passed, but the physical pilot validator failed because no complete six-method cell existed at 1.0 Msun.

For the css row, representative native TOV samples were:

| Scan index | Pc (Pa) | M (Msun) | R (km) |
|---:|---:|---:|---:|
| 0 | 7.3959e30 | 0.671717 | 636.227 |
| 3 | 2.7549e31 | 0.743214 | 1231.731 |
| 4 | 4.2706e31 | 0.684088 | 2777.204 |
| 5 | 6.6201e31 | 0.083621 | 51.056 |
| 14 | 3.4217e33 | 0.681138 | 11.645 |
| 18 | 1.9758e34 | 1.455607 | 11.138 |
| 23 | 1.7685e35 | 1.914837 | 9.320 |

Across all six matching prescriptions, the first sampled maximum lies in the narrow interval 0.743199–0.743409 Msun, with radius approximately 1231–1235 km. A later high-density rising interval reaches approximately 1.9147–1.91485 Msun at the scan endpoint. The very large excursion in radius and mass between scan indices 3–6 is common to all six methods at this coarse sampling.

These values are TOV sequence observations from a 24-point logarithmic Pc mesh. They do not show that the first peak is the physical maximum mass, do not demonstrate that the later high-density interval is stable, and do not establish that the sharp transition in the sampled M/R sequence is a physical EOS feature rather than insufficient resolution or a numerical issue. The fact that six bridges show similar topology suggests a common upstream core/low-density EOS contribution, but this remains a hypothesis.

## 3. Radial pulsation formulation and thermodynamic closure

For reference, in geometric units G=c=1 and metric
ds^2 = -exp(2 nu) dt^2 + exp(2 lambda) dr^2 + r^2 dOmega^2,
the standard first-order radial pulsation system can be written in terms of the dimensionless Lagrangian radial displacement xi and Lagrangian pressure perturbation DeltaP:

d xi/dr = (nu' - 3/r) xi - DeltaP/(r Gamma_1 P),

d DeltaP/dr =
{ exp(2 lambda) [omega^2 exp(-2 nu) - 8 pi P]
  + nu' (4/r + nu') } (epsilon+P) r xi
- { nu' + 4 pi (epsilon+P) r exp(2 lambda) } DeltaP.

For a regular centre, normalize xi(0)=1 and impose
DeltaP(0) = -3 Gamma_1(0) P_c xi(0).
For a true zero-pressure stellar surface, impose DeltaP(R)=0. These boundary conditions must not be copied uncritically to the current v5 endpoint: the current background is terminated at a finite positive pressure cutoff. The radial problem therefore needs either a documented low-density continuation to P=0 or a separately derived and independently tested finite-cutoff boundary condition.

The equilibrium barotropic index is
Gamma_eq = (epsilon+P)/P * dP/d epsilon.
The perturbative adiabatic index is instead
Gamma_1 = (epsilon+P)/P * (partial P/partial epsilon)_(s,Y_i).
They coincide only under the adopted instantaneous-equilibrium/barotropic perturbation closure. For cold beta-equilibrium matter, the frozen-composition perturbation can have a different Gamma_1. The current projected EOS table does not by itself provide all composition derivatives needed for a microphysical frozen-composition stability calculation. A result computed with Gamma_1=Gamma_eq must therefore be labelled “barotropic-equilibrium radial stability”, not universal physical stability.

At a genuine sharp first-order interface the radial variables require their own junction law. In the usual limiting cases, [DeltaP]=0 and:
- slow conversion: [xi]=0;
- rapid conversion: [xi - DeltaP/(r P')]=0.
More general reactive junction conditions require composition variables and their gradients. These conditions are not to be imposed at the historical smooth-bridge matching interval. The native sharp-Maxwell PR #2/#3 track remains separate.

## 4. New diagnostic implemented in the mass-dependent track

The code now partitions each contiguous valid TOV scan block into monotonic M(Pc) candidate segments:
- rising;
- falling;
- flat/degenerate;
- unresolved single-sample segments.

Local maxima/minima are shared endpoints of adjacent segments. Failed TOV samples and omitted scan indices are hard barriers. Each row carries the explicit label
radial_stability_status = unresolved_requires_radial_mode_spectrum.

The original initial-rising-segment selector is retained only as a conservative target-mass search proxy. Its output metadata and log labels are being renamed so its first local peak is not called Mmax and no unverified segment is called radially stable. Candidate segments are written to candidate_mass_segments.csv and are not used for Cowling spectrum roots.

The pilot scan is increased from 24 to 96 logarithmic Pc samples for the same six EOS/method rows. This experiment asks whether the observed topology is resolved with materially improved sampling. It does not itself assign radial stability.

## 5. Acceptance gates before secondary-branch spectra

1. Inspect the 96-point TOV samples for all six methods, including solver failures, radius, mass, central pressure, and EOS-cell/event locations.
2. Demonstrate convergence of every candidate maximum/minimum location and M/R values under a second grid refinement, for example 96 to 192 points, before interpreting the topology.
3. Reproduce selected TOV points with an independent implementation of the same raw native-knot P(mu) cubic-Hermite interpolation; a generic P(epsilon) spline is not an equivalent comparator.
4. Implement the radial eigenproblem with an explicitly selected Gamma_1 closure, centre series, and physically admissible surface condition. Verify the radial solver against published reference tables such as Kokkotas & Ruoff and/or Gondek et al.
5. Track the nodeless omega_0^2 along Pc and identify zero crossings. Radial-stability labels are assigned only after eigenvalue and mode-node checks converge.
6. If a true sharp phase boundary is present, implement and test its perturbative junction law separately; do not use smooth-bridge continuity or infer slow/rapid conversion from the equilibrium EOS alone.
7. Only after a radially admissible mass branch is established should mass-target solving and f/p1 mode scans be attempted on that branch.

## 6. Current disposition

- The 24-point six-model pilot: baseline M/R reproduction passes; no 1.0 Msun complete pair; zero admitted f/p1 spectra.
- The forward Cowling tangent mock and fixed-mass response mock pass their scoped CI tests but are not physical-star calculations.
- The 96-point scan and candidate-segment output are implemented in the next code increment; exact-head Actions are pending.
- No value in this note is a published maximum mass, a certified radial eigenfrequency, or a new physical f/p1 result.

## References

1. S. Chandrasekhar, “The Dynamical Instability of Gaseous Masses Approaching the Schwarzschild Limit in General Relativity,” ApJ 140 (1964) 417–433.
2. J. M. Bardeen, K. S. Thorne & D. W. Meltzer, “A catalogue of methods for studying the normal modes of radial pulsation of general-relativistic stellar models,” ApJ 145 (1966) 505.
3. K. D. Kokkotas & J. Ruoff, “Radial oscillations of relativistic stars,” A&A (2001), DOI 10.1051/0004-6361:20000216.
4. D. Gondek, J. L. Zdunik & P. Haensel, “Radial pulsations of proto-neutron stars,” A&A 325 (1997) 217; see also arXiv:astro-ph/9901167 for relativistic radial-mode equations and mode evolution.
5. M. G. Alford, S. P. Harris & P. S. Sachdeva, “On the Stability of Strange Dwarf Hybrid Stars,” ApJ 847, 109 (2017), arXiv:1705.09880.
6. R. Essick, “Exotic Stable Branches with Efficient TOV Sequences,” arXiv:2405.05395 (2024).
7. P. B. Rau & G. G. Salaben, “Non-equilibrium effects on stability of hybrid stars with first-order phase transitions,” arXiv:2309.08540v2.


## Exact 96-point pilot: observed mass-sequence topology (2026-11-09)

Exact-head run: https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42611801783  
HEAD: dd5a33e521ecde554945c82b79cfd534ed91b306  
Artifact: 13670407745; ZIP SHA-256 ecafc110e8b9d0196dfb4bb341da2ee54273f5d61767cca7c36e502a40737b7d

Execution facts:
- All six model rows completed 96/96 valid TOV integrations; no failed scan samples.
- All six frozen-v5 reference M/R reproductions passed at near machine precision.
- Candidate-segment tests: 10 passed.
- Fixed-mass/Cowling mock tests: 6 passed.
- 48 requested f/p1 rows were unavailable (0 valid); the physical validator failed at the required 1.0 Msun complete six-method cell. This is an expected scientific hold, not a code/test failure.

The 96-point first sampled local maximum is now consistently near:
- M = 0.7650195–0.7651069 Msun;
- Pc = 0.226758–0.227326 MeV fm^-3;
- R = 1879.8–1889.0 km.

By comparison, the 24-point run placed the first sampled maximum at M = 0.74320–0.74343 Msun, Pc = 0.17195–0.17236 MeV fm^-3, and R = 1232–1235 km. The sizable movement between 24 and 96 samples proves that the 24-point maximum was not a converged turning-point estimate. The 96-point estimate is still grid-resolved only, not a converged continuum extremum.

For five methods (css, cubic_F, cubic_P, quintic_F, quintic_P), the candidate pieces are rising/falling/rising/falling. The maxwell compatibility method (which is not an exact Maxwell transition and has no density jump) displays rising/falling/rising/falling/rising/falling. Its two intermediate turning candidates are:
- local maximum near Pc = 1.11640 MeV fm^-3, M = 0.216841 Msun, R = 17.856 km;
- local minimum near Pc = 1.89763 MeV fm^-3, M = 0.202879 Msun, R = 17.315 km.

These lie inside its smooth near-pressure-plateau bridge, which is documented as P_transition,left = 0.25468 to P_transition,right = 3.308 MeV fm^-3, n_left = 0.054 to n_right = 0.2 fm^-3, with cs2_min = 9.0712e-4, and zero density jump. The other five bridges do not show this extra pair of turning candidates at 96-point resolution. This is an observed method-dependent sequence feature, not yet a radial-stability or mode-spectrum result.

All five non-maxwell methods have a later rising candidate piece from approximately M = 0.0713 to 1.9196 Msun over Pc = 0.312 to 892.7 MeV fm^-3, followed by a small final decline to about 1.9149 Msun at the scan endpoint. The maxwell row has the additional local fold noted above before its later high-density rise. No candidate segment is labelled stable; the JSON/CSV explicitly retain radial_stability_status=unresolved_requires_radial_mode_spectrum.

A notable coincidence is that the common first maximum sits shortly below the smooth bridge's left boundary at 0.25468 MeV fm^-3. This suggests a common crust/bridge-entry response, but the 96-point grid still moves the first maximum significantly relative to the 24-point grid. The hypothesis must be tested with 192-point refinement and a separate radial eigenvalue calculation before physical attribution.


## 7. Next refinement

The pilot workflow is changed from 96 to 192 logarithmic central-pressure samples for the same six EOS/method rows. Acceptance requires:
1. each model completes all 192 native TOV integrations;
2. first-fold Pc, M and R changes are compared with the 96-point values, not merely the discrete segment count;
3. all local candidate turns, especially the method-specific maxwell bridge folds, are recorded without stability labels;
4. the exact same EOS SHA, v5 source build key, r0, radius limit, step cap and ODE tolerance are retained;
5. no Cowling spectra are admitted until radial omega_0^2 is computed with a declared perturbative Gamma_1 closure and a valid surface condition.


## Exact 192-point pilot: observed mass-sequence topology (2026-11-09)

Exact-head run: https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42615970377  
HEAD: 40177800d4484a19149485581e0f2d5bd69594a9  
Artifact: 13673520039; ZIP SHA-256 e460c295af38a70444e6139fb82aa037b3be4718e29d172db7be31d358abbff5

Execution:
- All 6 model rows completed 192/192 valid TOV integrations (1152/1152 total; no failed samples).
- All frozen-v5 M/R reference reproductions passed at near machine precision.
- Candidate-segment unit tests: 10 passed.
- Fixed-mass and Cowling tangent mock tests: 6 passed.
- Valid f/p1 spectrum rows: 0/48; complete six-method mass cells: 0. The scientific validator again fails correctly at the required 1.0 Msun paired cell.

First sampled maximum across methods:
- M = 0.7652434–0.7652747 Msun;
- Pc = 0.2248758–0.2254394 MeV fm^-3;
- R = 1849.995–1858.853 km.

Relative to the 96-point results, the first-maximum mass changes by only about 2.2e-4 in relative terms, while Pc shifts by roughly 0.8% and R by about 1.6%. This is promising but not sufficient to call the structural turning point converged, especially because the first subsequent minimum remains poorly resolved in radius: for css, the sampled minimum changes from about 418 km at 96 points to about 195 km at 192 points, while its mass remains about 0.0714–0.0716 Msun. The transition region has a sharp response in M(Pc), and radius convergence lags mass convergence.

The method-specific maxwell bridge fold persists and becomes better localized:
- intermediate local maximum: M = 0.2179489 Msun, Pc = 1.0415144 MeV fm^-3, R = 17.8614 km;
- intermediate local minimum: M = 0.2028220 Msun, Pc = 1.9619452 MeV fm^-3, R = 17.2259 km.

At 96 points the corresponding candidates were M = 0.2168408 Msun at Pc = 1.1164034, R = 17.8562 km and M = 0.2028787 Msun at Pc = 1.8976322, R = 17.3147 km. The additional pair therefore persists under doubling the grid, although its Pc locations still move by several percent. It remains a kinematic feature of the smooth near-pressure-plateau compatibility bridge (not an exact Maxwell jump), not a radial-stability result.

For the five other methods, the candidate direction pattern remains rising/falling/rising/falling. For the maxwell compatibility row it remains rising/falling/rising/falling/rising/falling. The late high-density sampled maximum remains near 1.9195–1.9196 Msun at Pc about 893.7 MeV fm^-3 and R about 9.53 km, followed by a small decline at the upper scan edge; this is not certified Mmax because the high-pressure endpoint is not a demonstrated EOS-domain/sequence boundary.

The 192-point result shows that:
1. the first maximum mass is becoming stable across grid refinement, but its radius and the following minimum are not;
2. the maxwell-labelled smooth bridge has a reproducible additional M(Pc) fold absent in the other five prescriptions at these resolutions;
3. no segment can be assigned radial stability without the nodeless radial eigenvalue omega_0^2.

The pilot is therefore refined once more to 384 central-pressure points. At the same time, an independent SciPy DOP853 TOV calculation is added. It reconstructs the exact native-knot cubic-Hermite P(mu) interpolation (with raw n=dP/dmu knot derivatives) and integrates the SI TOV equations directly in radius. Five selected Pc values cover the first fold, the early falling interval, both maxwell bridge folds and a high-density reference. The cross-check target is relative M and R agreement <= 1e-6 against the frozen v5 cell_solve TOV-only backend. This is a structure/TOV benchmark, not a radial-mode or Cowling-spectrum validation.


## Exact 384-point pilot and SciPy TOV comparison (2026-11-09)

Actions run: https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42616463263  
HEAD: 6fb68cf1204b9656aa526f437a380e3c2e340a58  
Artifact: 13671829050; ZIP SHA-256 a8bd61678ec74e1dea5451718f11386b44bb40162cad6ab29e10f1fe6268befc

### Native TOV sequence

- Six selected EOS/method rows completed 384/384 integrations each: 2304/2304 valid TOV points, no solver failures.
- All six frozen-v5 reference M/R anchors reproduced to near machine precision.
- Candidate-segment tests: 10 passed.
- Fixed-mass/Cowling mock tests: 6 passed.
- Native-Hermite EOS tests: 4 passed.
- New mass-dependent f/p1 spectra: 0/48; no complete six-method mass cell; every candidate segment remains radially unresolved.

The first sampled maximum across the six methods is:
- mass: 0.7652866–0.7652926 Msun;
- central pressure: 0.2239479–0.2245094 MeV fm^-3;
- radius: 1835.56–1844.27 km.

From 192 to 384 samples the first maximum mass changes by only about 2.6e-5 relative, while its radius decreases by about 0.75% and the pressure location by about 0.4%. The mass estimate is now highly stable under this grid refinement; radius is not yet within a 0.1% structural criterion. The 24-point estimate was severely under-resolved, and the 96-point estimate had a much larger radius.

The first post-maximum minimum remains poorly resolved:
- css: M=0.0710791 Msun, Pc=0.3153889 MeV fm^-3, R=322.49 km;
- cubic_F: M=0.0711298 Msun, Pc=0.3153347, R=363.82 km;
- cubic_P: M=0.0711171 Msun, Pc=0.3153452, R=357.44 km;
- maxwell compatibility: M=0.0711879 Msun, Pc=0.3160936, R=258.24 km;
- quintic_F: M=0.0710918 Msun, Pc=0.3154177, R=341.17 km;
- quintic_P: M=0.0710763 Msun, Pc=0.3155437, R=317.60 km.

This minimum radius is non-monotonic across 96→192→384 and is still not converged. No interpolated M/R fold is used for mode roots.

The legacy maxwell compatibility bridge retains an intermediate max/min pair:
- max: M=0.2180954 Msun, Pc=1.0606304 MeV fm^-3, R=17.8517 km;
- min: M=0.2028434 Msun, Pc=1.994654 MeV fm^-3, R=17.1793 km.
Their masses are nearly stable compared with 192 points, while Pc locations move by roughly 1.8% and 1.7%. The metadata still identifies this method as a smooth near-pressure-plateau construction with zero density jump; it is not a true Maxwell phase transition and no Seidov criterion applies. Five other methods retain four M(Pc) candidate segments; this method retains six.

The late high-density sampled maximum is around M=1.91955 Msun, Pc=918.07 MeV fm^-3 and R=9.504 km, followed by a small decline at the upper grid edge. It is not a certified global Mmax or radially stable branch.

### Independent SciPy TOV check

The 384-point CI also ran:
- 4 native-Hermite EOS tests: PASS;
- five selected direct-r SciPy DOP853 TOV comparisons against frozen v5 production TOV;
- the first comparison used SciPy rtol=2e-10, max_step=2000 m and native Fortran settings rtol=1e-9, step_cap=0.05.

The raw comparison reached maximum relative mass difference 3.60164e-7 and radius difference 1.46572e-6. The predeclared acceptance was 1e-6 for both M and R, so the radius criterion narrowly failed. No tolerance was widened. The failure did not invalidate the 2304 successful native integrations or the six frozen-anchor reproductions; it means the independent structural comparison was not yet accepted.

The first version of the comparison script raised before writing its casewise receipt. This is a diagnostic/representation defect, not evidence that the two TOV formulations disagree physically. The next commit repairs this by always serializing per-case coarse/tight SciPy and default/tight native results before checking acceptance. It tightens the SciPy direct-r integration and also refines the native q-step. The same 1e-6 criterion remains in force, now applied to tight SciPy versus tight native results. The next exact-head run must inspect the saved per-case receipt even if acceptance fails again.


### Exact-head 384-point and independent TOV receipt: run 42617251360

The exact-head run used commit `57511496a57e1a1cee6bafc5254de63f4f5ab514` and artifact `13673098470` (SHA-256 `0d376752ba610819d3ca59fceedbf50ea5eebc5d03b3f440e3b87135546b3964`). The 384-point pilot completed all 2304/2304 native TOV integrations, all six frozen-v5 M/R anchor reproductions, and the 10 mass-sequence, 6 fixed-mass/Cowling mock, and 4 native-Hermite tests. There were still 0/48 admitted f/p1 spectra and zero complete paired target-mass cells; the physical validator was skipped after the preceding step failed.

The saved independent SciPy receipt has status `PASS_REFINED_SAME_EOS_INDEPENDENT_TOV`. Across five selected central pressures, the maximum tight SciPy-versus-tight-native relative differences are:
- mass: `1.6232670240598032e-8`;
- radius: `1.1432510017123261e-7`.

Both are below the unchanged `1e-6` acceptance threshold. Native baseline-to-refinement differences are at most `1.299910193106295e-10` in mass and `7.057766721929021e-11` in radius. The coarse SciPy comparison remains deliberately diagnostic and is not the acceptance comparator. This closes the **structural TOV cross-check only** for the sampled points and shared finite-pressure surface convention.

The overall Actions run nevertheless concluded `failure`, not green. The numerical script completed and wrote the passing receipt, then its CLI summary raised `KeyError: 'max_relative_mass_difference'`: the printer referenced obsolete v1 field names, while the persisted v2 receipt correctly stores nested fields `max_relative_difference_scipy_tight_vs_native_tight`, `max_relative_difference_coarse_vs_native_reference`, and `max_relative_native_refinement`. This is an implementation/reporting defect after successful numerical acceptance, not a TOV discrepancy. A regression test now exercises the CLI against the v2 receipt schema. The six-method validator remains intentionally unrelaxed.

Scientific boundary unchanged: the first rising-proxy endpoint is only about 0.7653 Msun, all candidate segments have `radial_stability_status=unresolved_requires_radial_mode_spectrum`, and no 1.0 Msun six-method cell or new Cowling spectrum exists. Do not infer stability from the high-density rising interval or from the structural TOV agreement.


### Surface-boundary continuation diagnostic from the raw projected EOS (2026-11-09)

Source: `oscillation_release_v1/data/eos_projected/GPPVA_FSU2__BBB_BHF2__css.csv`, Git blob `5826c044e18054a8709dacde4bc648559cfb5bb9`, 323 rows. The first five raw knots were checked in all six projected method files (`css`, `cubic_F`, `cubic_P`, `maxwell`, `quintic_F`, `quintic_P)); the first-five tuples ((n_b,P,\epsilon,\mu_b)) are exactly identical across all six methods. Thus this particular outer-crust continuation can be shared across the six-method pilot, subject to the remaining EOS-source/provenance checks.

At the finite-pressure cutoff:
- (P_{\rm cut}=9.000000079422\times10^{-12}\;\mathrm{MeV\,fm^{-3}});
- (\epsilon_{\rm cut}=4.354134734530827\times10^{-7}\;\mathrm{MeV\,fm^{-3}});
- (n_{b,\rm cut}=4.67962\times10^{-10}\;\mathrm{fm^{-3}});
- (\mu_{b,\rm cut}=930.4654511545\;\mathrm{MeV});
- right-cell Hermite (c_s^2=3.1507096132423\times10^{-5});
- (\Gamma_{\rm eq}=(\epsilon+P)c_s^2/P=1.52432195647).

The right-sided local (\Gamma_{\rm eq}) values at the first five knots are approximately (1.52432,1.50950,1.49434,1.47885,1.46155). A log-log least-squares fit (\log P=a+\gamma\log\epsilon) across successively wider low-density windows is not a stable global prescription:
| First knots | fitted (\gamma) | RMS residual in (\log P) |
|---:|---:|---:|
| 5 | 1.58002 | 0.00241 |
| 10 | 1.54586 | 0.01140 |
| 20 | 1.46579 | 0.05136 |
| 40 | 1.39937 | 0.08822 |
| 80 | 1.36705 | 0.10260 |
| 160 | 1.32284 | 0.13340 |

A candidate *local tangent* zero-pressure continuation is:
[
P=K\epsilon^{\gamma_s},\qquad
\gamma_s=\frac{\epsilon_{\rm cut}}{P_{\rm cut}}c_{s,\rm cut}^2,
\qquad K=\frac{P_{\rm cut}}{\epsilon_{\rm cut}^{\gamma_s}}.
]
It matches (P) and (dP/d\epsilon=c_s^2) at the cutoff and has (P,\epsilon\to0) together. This is a mathematical C1 closure candidate, not a physically sourced crust EOS. A thermodynamically admissible implementation must additionally continue (n_b) using (d\epsilon=(\epsilon+P)dn_b/n_b), and set (\mu_b=(\epsilon+P)/n_b), preserving (dP= n_b d\mu_b). Do not just extend (P(\epsilon)) while leaving (n_b,\mu_b) inconsistent.

The continuation must be sensitivity-tested against alternative smooth continuations / fit windows, TOV radius and mass changes, and radial eigenvalues. The present data support a C1 tangent baseline but do not justify selecting a unique exterior EOS. First-five-knot identity across the six methods indicates that the outer continuation itself should not create a spurious matching-method effect if shared identically.

The first radial-stability implementation should use (\Gamma_1=\Gamma_{\rm eq}) only as an explicitly labelled *barotropic-equilibrium* closure. Frozen-composition (\Gamma_1=(\epsilon+P)(\partial P/\partial\epsilon)_{s,Y_i}/P) requires composition derivatives that are not present in this projected table.
