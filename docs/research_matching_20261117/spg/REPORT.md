# Conditional SPG generator closure — 2026-11-17

Baseline: `478c7b5de1514ac512f91156f36551478f36f688`. Selected work: R3–R8/C2–C8.

The public native generator now accepts the explicit, immutable variant
`public_composition_inferred_20261116`. All five models have actually generated
129 stationary rows from their independently solved neutron drip to their
fixed-density cluster/uniform energy crossing: 645 rows, native composition
readback, all native consistency hard gates passed. This is a usable conditional
SPG EOS reconstruction. The default remains `AME2408_candidate`; the new variant
is never described as the original author mass calibration.

## Input and scientific scope

The source JSON is
`data/research_blockers_20261116/spg_inverse/inferred_surface_parameters.json`,
SHA256 `7d6974d11c5986eb0b10c04aa6c085b3f1508711735edd9f895abca0b7e0cfcd`.
Its bytes are checked before request construction. The public receipt carries
classification, full-A assumption, original-calibration unresolved status and
source beta-stationarity discrepancy. Dispatch now also checks that the surface
input did not change during execution. All 15 pinned CompOSE source files were
rehashed and matched.

The previous 760-row mechanical/virial inversion and full-surface beta derivative
finite differences are reused from `docs/research_blockers_20261116/spg_inverse/REPORT.md`
and `docs/research_coexistence_20261116/SPG.md`. No source search, refit, or repetition
of that 760-row derivative campaign was performed. Those data establish a precise
conditional surface prescription and a source-state beta residual, not hidden
input identity. In particular, the source residual pattern approximating an
extra `(1-f)` surface-derivative factor remains an observation about the source
states. No such factor was inserted into the generator.

The new forward calculation retains the original spherical CLD energy, uniform
electrons, common nucleon rest energy 938.9 MeV, hbar*c=197.3269804 MeV fm,
e²=1.43996448 MeV fm, and no muons. Density is fm^-3; energy density and pressure
are MeV fm^-3. The inverse input was identified on proton-free source states;
the forward solver permits gas protons where its insertion condition requires
them. At the generated endpoint M1/M3/M5 have positive free-proton fractions
2.4044e-4, 2.0683e-4 and 7.0933e-8 respectively. M2/M4 remain proton-free in all
129 rows. This distinction is explicitly retained.

## Actual EOS and variational checks

| Model | drip n_b [fm^-3] | crossing n_b [fm^-3] | max KKT [MeV] | min interior scaled curvature |
|---|---:|---:|---:|---:|
| M1 | .00023317369094 | .07653917944431 | 1.806e-8 | .027865 |
| M2 | .00024162193293 | .08921915113926 | 1.931e-8 | .089286 |
| M3 | .00023554102995 | .08115016114712 | 1.964e-8 | .030822 |
| M4 | .00024756028502 | .09882555854321 | 1.997e-8 | .085235 |
| M5 | .00023611274269 | .07303679645062 | 1.976e-8 | .035018 |

Every requested row is in `data/research_matching_20261117/spg/native_129_repaired/SPG_M*/generated/`.
The native writer independently re-reads density, thermodynamics, phase, every
composition value and all seven masks; the frontend then reads the exported
tables again. Baryon/charge closure is within 1.11e-16/5.55e-17. Muon insertion
margins remain positive, with the minimum 10.8315 MeV.

The drip routine's zero curvature field is a placeholder, not a measured zero
eigenvalue. This was addressed by a separate explicit-energy Hessian calculation
in the vacuum tangent variables `(n_i,y)`, each scaled by 0.1. Three finite steps
per model give minimum eigenvalues 2.0692–2.3536, all positive. Raw results are
`drip_curvature.csv`. Interior curvature and drip KKT retain their distinct scopes.

The fixed-density crossing has |delta epsilon| <=6.36e-12 MeV fm^-3. Its pressure
jumps are .00445–.00807 MeV fm^-3 and chemical-potential jumps .06060–.08167 MeV.
It is the explicitly selected source-style energy crossing; it is not a Maxwell
common tangent. No cross-phase interpolation or globally optimal phase diagram
is certified by this campaign.

## Independent thermodynamics and sensitive derivatives

`spg_finish_validation.f90` independently assembles the explicit bulk+electron+
Coulomb+surface energy and eliminates the radius through its virial relation.
It reconstructs `(t,n_i,y,v)` from exported composition and radius, where
`t=Yh`, `y=Z/A`, `n_i=3A/(4pi R³)`, `v=y_g/y`,
`f=t n_b/n_i`, `n_g=n_b(1-t)/(1-f)`. It shares only the native nuclear-bulk and
free-lepton evaluations with production; surface evaluation, coordinate mapping,
energy assembly and finite differences are independent implementations.

All 645 explicit energies agree within 8.53e-14 MeV fm^-3. Five step sizes
1e-3,3e-4,1e-4,3e-5,1e-5 test fixed-coordinate
`mu_b=partial epsilon/partial n_b` and `P=n_b² partial(epsilon/n_b)/partial n_b`.
At relative step 1e-4, maximum errors are 3.07e-8 MeV and 6.75e-10 MeV fm^-3.
The smaller steps show the expected roundoff increase; no tolerance was widened.

Forty selected states include drip, the next row, interior quartiles, the last
two rows, maximal-KKT and closest gas-proton-threshold points. Their four relevant
constrained derivatives were evaluated directly from the independent energy.
At neutron/proton vacuum boundaries the kinetic term contains a density^(5/3)
contribution, so a one-sided derivative has leading finite-step error of order
h^(2/3), rather than an analytic fourth-order remainder. Observed maximum t/v
errors decrease from .001938/.003102 to .0000900/.0001414 MeV as h decreases from
1e-3 to1e-5. The t convergence order between the last two steps is approximately
2/3. Smooth n_i/y derivatives reach roughly 1e-7/1e-8 MeV before roundoff grows.
These finite-step errors are reported, not mislabeled as KKT failures.

At the same 40 points, 800 neighboring states were separately reoptimized
(five steps, four neighbors each). Maximum |d epsilon/dn_b−mu_b| decreases from
3.05e-5 to2.27e-7 MeV. The pressure derivative residual reaches2.99e-9 MeV fm^-3
at h=1e-4 and grows to2.23e-8 at h=1e-5. Sound-speed estimates remain positive
and subluminal. The largest change between the last two sound-speed estimates
is1.563e-5 in c_s²/c², at neutron drip; no finer sound-speed accuracy is claimed.
Raw steps and errors are retained in `validation_refined/{audit,thermo}.csv`.

## Reproduced and repaired implementation failure

The first actual M5 campaign stopped at `CLUSTER_BRANCH_LOWER_UNCONVERGED_ATTEMPT`.
At n_b=.073062500000000058, a finite cluster had E/A−m=8.00578802360003 MeV and
KKT5.51e-9 MeV. A lower failed trial had t=1.0127478667e-8,
f=1.012746791e-8 and E/A−m=8.005766724990087 MeV. The independent uniform beta
solution is8.005766724989599 MeV. They differ by4.89e-13 MeV.

The previous uniform-limit test required t<=1.0001e-8, so this optimizer stopped
just above the artificial coordinate floor and was misclassified. The minimal
repair requires both negligible volume fraction f<1e-7 and agreement with an
independent uniform bulk energy within1e-9 MeV. It does not accept the failed
trial as a stationary cell. A materially finite phase or a dilute energy-distinct
trial still fails this new classifier. The true finite-cluster lower-attempt
HOLD remains. No bulk, surface, derivative, KKT or curvature formula changed.

The exact crossing regression returned exit3 before the repair and exit0 after
it, with three classification controls and |delta epsilon|=4.31e-12 MeV fm^-3.
The original failed M5 run and diagnostic attempt coordinates are preserved.
A separate build attempt returned75 because another task changed unrelated
sources during compilation; it published no binary. The completed repaired run
used a frozen source snapshot, recorded in its input manifest, and unchanged
physical dependency hashes. Existing runtime array-temporary/IEEE-invalid search
warnings are retained; accepted states and diagnostics are finite.

An independent review also reproduced missing integrity coverage for optional
`SUPPORT.json`. The optional-file hash and sealed-missing checks now include it.
The reviewer supplied three red tests; the final targeted generator suite,
including the six new inferred-input tests, ran89 tests successfully.

## Reproduction and handoff

From the repository root with OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1:

```sh
bash fortran_eos/tests/spg_finish_build.sh /tmp/spg_finish_build
/tmp/spg_finish_build/spg_finish_regression
python3 fortran_eos/tests/spg_finish_run.py /tmp/spg_finish_build /tmp/spg_finish_run
python3 fortran_eos/tests/spg_finish_validate.py /tmp/spg_finish_build /tmp/spg_finish_run /tmp/spg_finish_validation
/tmp/spg_finish_build/spg_finish_drip_curvature < /tmp/spg_finish_validation/inputs.txt
PYTHONPATH=python_frontend pytest -q python_frontend/tests/test_spg_inferred_generation.py
```

The public Python entry is `run_generation('SPG_M5', fresh_output,
variant='public_composition_inferred_20261116', points=129)`; the CLI variant
choices are derived from the same catalog. The root integration task owns its
final public CLI smoke and core/crust matching campaign. The direct native full
campaign and independent diagnostics above are already complete.

`data/research_matching_20261117/spg/EVIDENCE.json` binds source/input identities,
actual return codes, all raw evidence and exact selected states. The repaired
run receipt SHA256 is
`955dec02748f124a72fab3d0010e379c702e1475aacec25ca3cfba7d3b451ea4`;
the frozen production SPG source SHA256 is
`ac02b84295210037de44d3a34d0d90f32c9f1e6e3001638ff0bc609b96ec2dbc`.

Status: implementation-verified and numerically checked conditional finite EOS
construction. Original nuclear-mass calibration/source-state equivalence remains
unresolved as before. It does not block use under the explicit candidate identity.
No all-geometry global minimum, original author calibration or complete joined
crust is claimed. Production edits in this workstream are frozen.
