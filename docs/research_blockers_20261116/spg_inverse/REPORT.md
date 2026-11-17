# SPG: public-composition inverse recovery and new variational states

Baseline: `a49aaf9f82ddfda484e0f965aefdfe73aa7918f2`. Work date: 2026-11-16.

**Result.** The missing surface coefficients can be reconstructed conditionally from public CompOSE composition and the published CLD energy, without fitting final pressure or energy density. All 760 eligible proton-free rows have one detected mechanical root under the full dense-phase interpretation of A. Their virial-inferred surface values closely follow the two-parameter Scurto formula, including held-out rows. The resulting coefficients produce 15 new stationary native-Fortran equilibria. This closes the absence of an executable public-data surface prescription; it does not authenticate the original nuclear-mass fit or establish equivalence to every source equilibrium.

The source composition does not satisfy the full beta-stationarity equation of the published energy to its printed precision. Most of its residual is proportional to `f*3 sigma'(y)/(R ni)`. The final thesis equations do not contain the corresponding `(1-f)` multiplier; this is an observed numerical pattern, not proof of an author-code defect. Production physics was not changed.

## 1. Source, conventions, and domain

Literature: Scurto, Pais & Gulminelli, [arXiv:2402.15548v1](https://arxiv.org/html/2402.15548v1), Eqs. (21)–(25). The source-input lane also inspected the original TeX and Scurto thesis Appendix A, especially A46–A48. Its independent read-only formula check agreed with the mechanical and beta derivatives below. That check did not find a published `(1-f)` surface-derivative prescription or an explicit exported-A definition.

Public inputs are the pinned `eos_generator/eos_data/SPG_M{1..5}/eos.nb`, `eos.compo`, and `eos.thermo`; exact hashes are in `data/research_blockers_20261116/spg_inverse/summary.json`. Bulk coefficients are the exact archived rows bound in `data/research_generator_20261116/spg_author_candidates.json` and implemented by `compose_eos_gdfm_native.f90`. The source-input lane checked all 16 numbers for every model against that archive. This does not establish that the original internal solver used no additional unpublished digits.

At T=0, use spherical, sharp-interface, uniform-electron CLD cells, no muons and no gas protons. `ni`, `ng`, `nb` have units fm^-3, R has units fm, energies and chemical potentials MeV, energy density MeV fm^-3, and sigma/sigma' MeV fm^-2. Here sigma' means d sigma/d y; y is dimensionless. The common baryon rest energy is 938.9 MeV; hbar*c=197.3269804 MeV fm is retained in the bulk and lepton calculations, and e^2=1.43996448 MeV fm denotes the Coulomb coupling including the SI 1/(4 pi epsilon0). No neutron-proton mass splitting is inserted into this source model.

Raw inner rows have placeholder zero neutron/proton entries and a group slot that equals Ye, not a standard cluster-number abundance. They are not silently promoted to physical gas abundances. Under the explicitly conditional full-A map,

    y = Z/A,  t = Ye/y,  f = t nb/ni,
    ng = nb(1-t)/(1-f),  R = [3 A/(4 pi ni)]^(1/3).

The source imported outer/lower segment (raw group slot=1) is outside this inversion; no SPG equilibrium claim is made for it. Source M3 row366, nb=0.07943282347 fm^-3, explicitly has gas-proton fraction 0.00003 and is excluded from this proton-free branch. Its raw row is preserved in `summary.json`; it is not classified as a failed inverse solve. Eligible counts M1..M5 are 148,160,148,158,146.

## 2. Constructive derivation

Let epsI=eps_b(ni,y), epsg=eps_b(ng,0), and

    E = f epsI + (1-f) epsg + eps_e(f ni y) + Ccell + Scell,
    Ccell = 2 pi e^2 f Phi(f) R^2 (ni y)^2,
    Scell = 3 f sigma(y)/R,
    Phi = (2-3 f^(1/3)+f)/5,
    ng = (nb-f ni)/(1-f),  ne=f ni y.

All derivatives below hold nb fixed. The independent internal variables are (f,ni,y,R). R variation gives Scell=2 Ccell. Write C=Ccell/f. With muI=mu_nI+y(mu_pI-mu_nI), density and volume variation give

    (1/f) dE/dni = muI-mu_ng+mu_e y+2 C/ni = H,
    dE/df = epsI-epsg+mu_ng(ng-ni)+mu_e ni y
            + C[3+f Phi'/Phi].

Since ni muI-epsI=PI and ng mu_ng-epsg=Pg, subtraction eliminates chemical potentials, electrons and sigma:

    ni H - dE/df = PI-Pg-C[1+f Phi'/Phi] = M.

Thus solve M(ni)=0 from the public composition and known bulk model, then obtain

    sigma = 2 C R/3 = e^2 Phi A ni y^2

for the full-A convention. This is the mechanical-plus-virial inverse; it does not use public P or epsilon as fitting targets. H is an unfitted chemical-transfer diagnostic.

Composition variation supplies a second unfitted equilibrium equation:

    B = (1/(f ni)) dE/dy
      = mu_pI-mu_nI+mu_e+2 C/(ni y)+3 sigma'/(R ni).

The coefficient of sigma' is full `3/(R ni)`: the explicit f in the surface energy cancels the divisor f ni. A `(1-f)` multiplier is not obtained by this constrained variation. Phi is independent of y at fixed f, so its derivative belongs in M, not an extra term in B. Combining H-yB gives the neutron relation

    N = mu_nI-mu_ng-3 y sigma'/(R ni).

At a stationary cell M=H=B=N=0, although these four identities are not four independent constraints. Given a sigma value at a mechanical root, beta stationarity can instead infer

    sigma'_local = -(R ni/3)[mu_pI-mu_nI+mu_e+2 C/(ni y)].

For D=y^-3+(1-y)^-3 and sigma=sigma0(b+16)/(b+D), one locally obtains

    b_local = -D-D' sigma/sigma'_local,
    sigma0_local = sigma(b_local+D)/(b_local+16).

Local recovery is ill-conditioned when sigma' is small, at y=1/2 where D'=0, and for small y where b is subdominant to D. A source beta defect also contaminates it directly. Therefore final candidate coefficients are instead obtained from multiple mechanical/virial values using the exact linear relation

    1/sigma = alpha D + gamma,
    b=gamma/alpha,  sigma0=1/[alpha(b+16)].

Identifiability requires at least two distinct D values and nonzero, physically consistent alpha and b+16. At one y only the value sigma(y) is identified. No universal uniqueness theorem for M is asserted: the bounded scan tests a declared dense-phase density interval.

## 3. Actual native computation

`fortran_eos/tests/spg_inverse_probe.f90` implements both A conventions, the mechanical root search, surface inversion, linear fit, residuals, finite-nucleus energy minimization, synthetic checks, and forward states. Python scripts only prepare inputs, call native physics, and summarize results.

The search scans 160 cells over ni in `[max(0.03,nb*(1+1e-6)),0.25]` fm^-3 and bisects every detected sign change for55steps. Multiple roots are counted; the highest-density root is retained. Failure rows are kept with status2 and NaN derived values, never merged into valid-row residual statistics. A sign scan is not a proof that arbitrarily narrow pairs of roots do not exist.

Training uses nb<=0.04 fm^-3 and raw row index modulo5 !=0. Every fifth row plus every row above0.04 is held out. The original printed values are used directly. There is no pressure, epsilon, Q3 or Q4 optimization.

| Model | sigma0 [MeV fm^-2] | b | train/heldout | held-out max relative surface residual |
|---|---:|---:|---:|---:|
| M1 | 0.88358941 | 10.911369 | 97/51 | 2.35e-4 |
| M2 | 0.97186520 | 26.488074 | 100/60 | 6.20e-4 |
| M3 | 0.89299226 | 11.393456 | 95/53 | 1.50e-4 |
| M4 | 0.91922898 | 21.784891 | 95/63 | 1.00e-3 |
| M5 | 0.92029539 | 19.499992 | 97/49 | 3.12e-4 |

Digits in this table identify the experiment, not an accuracy claim on original author parameters. The train RMS relative surface residual is1.98e-5–2.51e-5. The scaled two-column linear-design determinant is0.152–0.269, so the sampled training data are not close to the single-y rank degeneracy. Parameter covariance under original hidden implementation errors is not determined.

The alternative gas-excess map was actually evaluated:

    ng=nb[1-Ye A/Z],
    f=(nb-ng)/(ni-ng), R=[3A/(4pi(ni-ng))]^(1/3),
    y=nb Ye/(f ni).

It produces25,32,26,36,24 missing mechanical roots within the stated density interval, respectively, and much larger remaining-row held-out surface discrepancies. M4 has two rows with multiple roots. All these rows and status codes are retained under `excess/`. This strongly favors the full-A interpretation within the tested model; it is numerical semantic evidence, not provider metadata authentication or a no-go theorem for every possible excess definition.

## 4. Residual mechanism and precision

The full published-energy beta residual reaches0.597,0.784,0.610,0.715,0.655 MeV on held-out M1..M5 states. Chemical-transfer residuals reach0.0771,0.0897,0.0783,0.0854,0.0770 MeV. The independently solved mechanical residual is near machine precision. These are separate outcomes: the good two-parameter surface fit does not erase the beta residual.

All eight corners of the printed rounding cell (Ye ±5e-6, A/Z ±0.0005) were propagated at the first, middle, and last eligible row of each model. Chemical residual changes are at most7.36e-4MeV and beta changes at most5.40e-3MeV in those15tested neighborhoods, well below their central residuals. Surface residuals at those representatives are comparable to their rounding sensitivity, including increasing uncertainty near large f. This is a measured15-state sensitivity, not a rigorous enclosure for every row or for unpublished bulk inputs.

The data show

    B approximately f Sbeta, H approximately y f Sbeta,
    Sbeta=3 sigma'_fit/(R ni), N approximately0.

Replacing B by B-f Sbeta reduces the all-row beta RMS to0.00428–0.00465MeV, with maximum0.00727–0.00810MeV. This corresponds algebraically to replacing the full surface term by `(1-f) Sbeta`. It is only a diagnostic of the source-table pattern. The neutron-equilibrium RMS is0.00048–0.00063MeV, while source -938.9Q4 and the free-electron chemical potential differ by at most0.00445MeV in the separate read-only check, compatible with printed Ye precision. No term was deleted from production to fit this pattern. Rounding, source bulk precision and source solver conventions remain distinct from an authenticated source-code defect.

Synthetic native-forward states with the deliberately known (sigma0,b)=(1.1,25) at nb=0.008 and0.03 for each model recover sigma0 within4.21e-8 and b within1.93e-6. The full beta residual is at most7.78e-8MeV; the hypothetical `(1-f)` variant leaves -0.453 to -0.0732MeV on those same stationary states. This tests the inverse/derivative implementation under a shared bulk routine; it is not an independent EOS-physics validation.

## 5. Constructive alternatives and independent mass data

A new native solve used the fixed inferred coefficients at three held-out source densities near0.004,0.02 and0.06fm^-3 per model. All15states are stationary, have accepted native curvature and insertion checks, and have maximum projected gradient1.14e-8MeV. Maximum relative differences from their public source compositions are Ye2.46%, A1.87%, Z2.68%. Thus the inferred prescription can drive a distinct variational EOS calculation while preserving the actual differences. These15states are not a density-wide EOS release.

An independent2408-nucleus AME2016 cohort was evaluated without refitting the inferred surface. The declared candidate functional is

    E/A = eps_b(n,Z/A)/n - 938.9
          +4piR^2 sigma(Z/A)/A +3 e^2 Z^2/(5 R A),
    R=[3A/(4pi n)]^(1/3), n in[0.04,0.24]fm^-3.

All12040model/nucleus minima are interior; |d(E/A)/dn|<=9.67e-12MeV fm^3. RMS total binding-energy residuals are53.61,34.11,56.20,50.64,35.78MeV for M1..M5; per-baryon RMS values are0.550,0.358,0.576,0.519,0.373MeV. No shell, pairing, curvature, exchange, electron-binding, or source-specific finite-nucleus prescription is silently supplied. This is an independent data prediction under a specified functional, and it does not recover the missing source fit recipe. `mass_summary.json` binds the exact2408-row AME input to its archived hash.

## 6. Reproduction and next use

From the repository root, with OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1:

```sh
bash fortran_eos/tests/spg_inverse_build.sh /tmp/spg_inverse_build
python3 fortran_eos/tests/spg_inverse_run.py /tmp/spg_inverse_build/spg_inverse_probe
python3 fortran_eos/tests/spg_inverse_run.py /tmp/spg_inverse_build/spg_inverse_probe excess
python3 fortran_eos/tests/spg_inverse_mass_check.py /tmp/spg_inverse_build/spg_inverse_probe
/tmp/spg_inverse_build/spg_inverse_probe synthetic > data/research_blockers_20261116/spg_inverse/synthetic_final.log 2> data/research_blockers_20261116/spg_inverse/synthetic_stderr.log
python3 fortran_eos/tests/spg_inverse_forward_check.py /tmp/spg_inverse_build/spg_inverse_probe
python3 fortran_eos/tests/spg_inverse_evidence.py
```

Compiler flags, source/executable hashes, actual branch counts and implementation checks are in `diagnostic_evidence.json` and `summary.json`. The compiler reports the existing SLSQP executable-stack requirement; runtime array-temporary warnings originate in the unchanged native solver. The synthetic log also preserves an IEEE_INVALID_FLAG report from the forward search; the accepted final states and their measured residuals are finite. This flag was not promoted to a scientific failure or suppressed.

`inferred_surface_parameters.json` is the concrete public-data research prescription. It is suitable as an explicitly labeled input to further variational construction. The remaining source-equivalence obligation is to resolve the beta surface-derivative pattern and original finite-nucleus calibration, not to force CompOSE equality. No production module, build script or default parameter selection was modified in this workstream.
