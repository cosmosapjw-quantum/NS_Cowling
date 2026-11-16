# Independent review: published BSk22 lower-inner representation

Status: PASS_SCOPED_PUBLISHED_ANALYTIC_VARIANT. This admits the equation implementation
and its separately identified thermodynamic representations. It does not admit an
independent microscopic ETF/Strutinsky minimization or chemical equilibrium of the
separately fitted composition.

The reviewed API is bsk22_analytic.state(n_b, thermodynamic_mode=...), restricted
to [2.69403e-4, .002] fm^-3. The author reference is bskfit18.f, SHA256
4bd3b716f04e40c69165fa83b5cd8ecb03c30aa1f2b45342fb070cef09a3a99c, version
2023-02-13. Its header and INCRST implementation explicitly set free protons to
zero below proton drip. The present density range is entirely below the BSk22
proton-drip switch .0592 fm^-3.

Source review confirms Z_cell=Z_cluster=40 in this range, A_cell=40/Ye,
A_cluster=A_cell*(1-Yn_free), Yh=1-Yn_free and Yp_free=0. These imply the
repository's baryon and charge closures. A 4097-point full-domain sweep gives
maximum closure residual 5.56e-17. Invalid densities NaN, infinity, nonpositive
values and one representable step outside each boundary are refused.

An independently transcribed scalar source energy was differentiated with mpmath
at 70 decimal digits in 33 cases. Maximum relative differences from the module
were 1.78e-15 for energy-derived pressure and 5.22e-15 for sound speed squared.
The paper's exact exponent 7/6 is retained; it differs from the author's rounded
Fortran literal 1.16667, an identified source representation difference.

A separate 40-digit outward interval-arithmetic evaluation covers the entire
declared density interval in 512 contiguous subintervals. For the energy-derived
branch it encloses dP/dn in [0.9122301,1.4422189] MeV and d epsilon/dn in
[939.2517125,941.4344157] MeV. It encloses c_s^2/c^2 in
[0.0009696495,0.0015354978]. Thus the analytic branch is compressible, stable and
causal on this interval. This statement concerns these published analytic
coefficient formulas; it is not a proof of calibration of an underlying ETF model
or a general proof concerning every floating-point implementation.

The separately fitted source pressure differs from n^2 de/dn by
-8.2646% to +12.0750% in the 4097-point sweep. Merely making the Euler identity
true by defining mu_B=(epsilon+P)/n does not remove this first-law mismatch.
Choosing pressure=n^2 de/dn repairs the barotropic first law, but leaves a
different mismatch: exported mu_B minus the independently fitted composition
mu_n ranges from -0.2905045 to +0.0623491 MeV. At the nominal source drip,
mu_B(energy derivative)=939.2928173 MeV, whereas the composition fit sets
mu_n=939.5654 MeV and Yn_free=0. This is a finite-error fit-consistency issue,
not evidence that the numerical derivative is broken.

The raw outer and inner endpoints are not a common Gibbs/pressure matching pair.
For the energy-derived branch, matching the outer endpoint chemical potential
alone leaves pressure larger by 1.81814e-5 MeV fm^-3; matching its pressure alone
leaves chemical potential lower by 0.056326 MeV. These are checks of these
unmodified representations, not an impossibility theorem for other models or
for a separately constructed thermodynamic reconstruction. No automatic join is
admitted by this review. An explicitly named boundary-anchored analytic variant
is being reviewed separately.

One implementation defect was found in the runner: it accessed outer['Ye'],
where the independently generated outer receipt uses outer['Y_e']. The observed
KeyError was reported and the owner repaired it. The numerical module was not
edited by this reviewer.

Reproduce the numerical and interval evidence:

```sh
python3 data/research_completion_20261116/outer_crust/bsk_review/review_bsk22.py
```

REVIEW.json contains exact reviewed code/source hashes, every high-precision
comparison case, interval enclosures, invalid-input results, composition checks,
and both unmatched boundary states. review.log retains the actual execution.
