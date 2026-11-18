# COWLING-MAXWELL-NATIVE-03: source-protected core + inner/outer crust and native f/p1

Date: 2026-11-19 KST. Work branch: `research/cowling-maxwell-native-core-crust-20261119`,
parent `93560502ed96ac92821e16cc7e24458a80fd885a`. This document is the
**scientific outcome and handoff**, not a blanket publication admission or global
claim about microscopic unification.

## 1. Requested science and success conditions

Question: Can original, or minimally first-law-corrected, cold crusts be combined
with a coherent homogeneous RMF core via *genuine sharp Maxwell* or *CSS*
matching; do resulting three-layer (core / inner / outer crust) barotropes support
TOV and nonradial relativistic Cowling oscillations? Vary cubic/quintic `P(mu_B)`
interpolation as **input representation only**; CSS is distinct physics and NOT an
interpolation variant. Microscopic phase minimization/unification is NOT a blocker
for the requested phenomenological full-star model.

Finite matching, full three-layer EOS, 1.4 M_sun TOV, and numerically admitted
f/p1 modes are separate assertions. Original provider tables, new same-functional
homogeneous cores, thermodynamically corrected provider pieces, and secondary
BBP table cannot be conflated. All `P`, `mu`, and `epsilon` in MeV/fm^3 or MeV,
density in fm^-3, radius in km, and frequency in Hz.

## 2. Native physics and numerical algorithm

Metric `(-,+,+,+)`; geometric enclosed mass `m=GM/c^2`, length km,
`omega_g=2*pi*f/c_km_s`, `h=integral dP/(epsilon+P)`,
`s=sqrt(h_start-h)`. Fortran DVODE integrates each phase's TOV `r,m` and
Cowling scaled `w=W/r^(ell+1)`, `v=V/r^ell` independently. No density-jump
interval is ever interpolated. Use separately evaluated one-sided
`(epsilon,P,c_s²)` from the original/declared source.

For `a=omega_g² exp(lambda-2Phi)r`, `b=dPhi/dr`,
`Q=av+bw`, and `H=epsilon+P`, conservation of Lagrangian traction requires
`[H Q]=0`. The native slow conversion map enforces `[w]=0`, the rapid map
`[v]=0` in this Cowling scaled convention; both independently check
`[DeltaP]=0`, with slow `[xi^r]=0` and rapid
`[xi^r-DeltaP/P'_0]=0`. In each case `P'_0=-H Phi'`.
See Tonetto & Lugones, PRD 101, 123029 (2020),
https://doi.org/10.1103/PhysRevD.101.123029 .

**Critical numerical fix:** A normalized surface shooting residual can change
sign through normalization/amplitude singularities without crossing zero.
The implementation now shoots from both center and surface to an interior
matching point and applies normalized determinant with an independent vector
matching residual. Only roots with small determinant, join residual, valid
node identity, surface BC, and one-sided junction residuals are admitted.
Netlib DVODE COMMON storage is non-reentrant: serialization is enforced with
a per-process mutex. Do not share one library state across concurrent threads.

Code:
- `fortran_eos/src/compose_eos_cowling_dvode.f90`
- `fortran_eos/src/compose_eos_cowling_tov_dvode.f90`
- `python_frontend/compose_eos/stellar_sharp_native.py` and
  `stellar_cowling_native.py`.
- `python_frontend/tests/test_cowling_maxwell_native_20261119.py`.

## 3. Same-functional GPPVA RMF surveys and interpolation

Source census: `same_rmf/SURVEY.json`; script
`scripts/survey_same_rmf_original_unified_20261119.py`.
The DDME2 homogeneous branch was recomputed using Fortran native DD-RMF;
TM1e and FSU2 use the existing Python RMF field-evaluation oracle with
published parameters. These are same-functional **generated homogeneous
states**, not additional downloaded original-provider core tables.
For published source crusts, preserve original `n,P,epsilon` with original
and first-law-corrected tiers separate.

| Same functional pair | Source tier | Native Maxwell | Native CSS | Finite Pt Maxwell / CSS [MeV fm^-3] |
|---|---|---|---|---|
| GPPVA(FSU2) + FSU2 | original crust + new core | PASS | PASS | 0.18318479756 / 0.25467999952 |
| GPPVA(FSU2) + FSU2 | first-law corrected crust + new core | PASS | PASS | 0.18319243303 / 0.25641857405 |
| GPPVA(TM1e) + TM1e | first-law corrected crust + new core | PASS | PASS | 0.33073692335 / 0.35869474845 |
| GPPVA(TM1e) + TM1e | original crust + new core | not found/admitted in scanned finite components | not found/admitted | — |
| GPPVA(DDME2) + DDME2 | original or corrected crust + new core | not found/admitted | not found/admitted | — |

The original FSU2 positive density jump is
`n_low=0.04827425147`, `n_high=0.04879235170` fm^-3 at
`mu_t=946.6140298453 MeV`. Native cubic and quintic EOS-input interpolants
both return positive-jump crossings. Their agreement demonstrates limited
interpolation robustness **only in scanned finite connected components**.
Never translate a negative sign screen into global proof of no phase
coexistence. Earlier 2026-11-18 work also covered 635 source-origin pairs,
7 finite double matches (28 cubic/quintic × Maxwell/CSS evaluations);
do not duplicate or relabel those historical receipts as new runs.

## 4. Origin-preserved GPPVA FSU2 three-layer Maxwell AND CSS

Original provider GPPVA(FSU2) crust: first **217 source BSk22 outer rows**,
next **52 source TF inner rows**. Actual input missing density interval
`0.000269021 -> 0.003 fm^-3` is explicitly filled by a validated
causal thermodynamic connector: `c_s²/c² ≈ (6.78e-4,1.02e-3)`.
Original outer `n` and `P` retained; pressure-primary first-law correction
to outer `epsilon/n` ≤0.004142 MeV (not zero). The source inner rows
are unchanged until the native Maxwell/CSS phase transition. Core states
computed from the same FSU2 functional, and are not claimed to be the
provider's unified homogeneous table.

At `n_c=0.30 fm^-3`:
- Maxwell `M=1.177195698 M_sun,R=14.085700712 km`; slow (f,p1)
  (1809.531483,4898.004388) Hz; rapid (1809.514982,4896.331069) Hz.
- CSS `M=1.177193520 M_sun,R=14.085695815 km`; slow
  (1809.573263,4901.321342) Hz; rapid (1809.526256,4896.508911) Hz.
- Both native TOV and independent SciPy Python DOP853 TOV parity are
  checked; all four mode node counts and junction residuals accepted.
- `data/cowling_maxwell_native_20261119/three_layer_original_fsu2/RESULT.json`
  and separate Maxwell/CSS `REPRESENTATION.json` are SSOT.

## 5. Official full GPPVA unified EOS, separated core / inner / outer

Actual official CompOSE ZIPs from
`https://compose.obspm.fr/download/1D/NS/RMF/{FSU2,TM1e,DDME2}/eos.zip`,
pages 215,221,218. An *independent HTTP download* of `eos.thermo`
matches the zipped member bytes for all three, **BUT their separately
published eos.zip_checksum.txt SHA-256 entries differ from observed ZIP
hashes**. Keep `PROVENANCE_ZIP_CHECKSUM_CONFLICT`; provider identity is
not established to the stronger advertised archive hash standard.
Never override original source bytes.

Unified source crust/inner/core source boundary:
- FSU2: 217 / 52 / 157 rows; source inner->core density
  0.054000001 -> 0.057862874 fm^-3, pressures 0.2546799995 ->
  0.2941949183 MeV fm^-3.
- TM1e: 217 / 86 / 152 rows; 0.089000002 -> 0.10700153 fm^-3,
  pressures 0.3455899893 -> 0.4270154463 MeV fm^-3.
- DDME2: original 415 rows, last two core densities DECREASE.
  Preserved all source bytes; created separate unmodified-in-order
  first-413-row projected copy. The applied central density is within
  this monotone source prefix; no high-density global claim.

Raw original split: native Maxwell fails (among others original Q3 Euler
thermodynamic inconsistency, or no overlapping µ region), and raw CSS
fails in these source intervals. First-law reconstruction **preserves
all source n,P knots**, changes `epsilon,mu` by integrating
`d(epsilon/n)/dn=P/n²` anchored to actual source (not arbitrary
energy shift); same phase model no global optimizer assumed.
Corrected split CSS finite native success:
- FSU2 Pt=0.2546799995, n-=0.054000001,n+=0.0548292595;
- TM1e Pt=0.3455899893, n-=0.089000002,n+=0.0890746507;
- DDME2 Pt=0.4067599943, n-=0.071999997,n+=0.0722158979.
- None of these is an admitted source-original **Maxwell** coexistence.
- All three corrected CSS 3-layer stars have actual phasewise Fortran TOV
  and slow/rapid f,p1 native Cowling modes and independent original-source
  corrected smooth Cowling controls.

At `n_c=0.30 fm^-3`, source-corrected smooth reference:
`FSU2 M=1.163885227,R=14.012872679,f=1815.759246,p1=4936.723394`;
`TM1e M=1.034416844,R=13.060202243,f=1986.544132,p1=5652.847591`;
`DDME2` in official `ORIGINAL_CONTROLS.json`.
Sharp CSS correction to f/p1 depends on conversion: for FSU2
slow f=1815.795053,p1=4940.596738 and rapid f=1815.758201,
p1=4936.716714 Hz, at the same central density.
Do NOT infer a Maxwell effect from differences between unlike EOSs.

## 6. BBP+BPS original source and conservative correction

Original historical BPS outer crust was extracted only below neutron drip
from original CompOSE APR attachment; BBP numerical 35-row Table 10.4
from secondary reference Diener (2008), arXiv:0806.0747,
`https://arxiv.org/abs/0806.0747`, PDF displayed page 97 (PDF page 109).
No sealed *primary* BBP original digital file obtained. Raw source values
have ~3 significant figures and non-monotone derived baryon chemical
potential. The unmodified digitization fails native Maxwell and CSS; this
is **not** evidence for global absence of phase equilibrium.

Independent of core, retain only 26 BBP inner rows up to density 0.0747
fm^-3. Preserve each original P knot; derive epsilon/n via the cold first law,
anchor at the published first energy density, and use causal BPS->BBP
connector after adjusting *only the first BBP density* by
`-4.16111482e-7 fm^-3`, within its printed ±5e-7 fm^-3 rounding.
Published original and derived tables are strictly distinct.
The resulting first-law correction reaches **5.3469 MeV/baryon**,
so it is an explicit non-negligible EOS change, not microscopic BBP
verification. Source PDF SHA and conversion details are recorded in
`bbp_bps/RESULT.json`.

For corrected BBP+BPS, native **CSS** finite matches pass with all
FSU2/TM1e/DDME2/GM1/PK1 cores; native Maxwell does not pass in the
scanned actual support. At central density 0.30 fm^-3, all five
complete 3-layer EOS have successful native TOV. Slow+rapid f/p1 node
admission was complete for FSU2/TM1e/PK1; DDME2 and GM1 originally
had slow p1 node-ID ambiguity despite valid TOV and rapid mode solutions.
At target **1.4 M_sun**, the newly rerun 10-case survey resolved
the full f/p1 spectra under a tighter TOV integration tolerance;
do NOT reinterpret that as a proof of global mode-family continuation.

## 7. Explicit fixed 1.4 M_sun branch calculation

`scripts/run_fixed_mass_maxwell_css_20261119.py`, output
`fixed_mass_1p4/RESULT.json`. Fortran native TOV
`rtol=3e-10` (tightened from 3e-9 after a measured ~8e-12 km
tiny negative mass increment at dense output, rejected by unchanged
monotonicity check). Test bracket `n_c in [0.25, 0.45] fm^-3`,
requires locally positive `dM/dn_c`. All **10** candidate EOS models
reach `M=1.4 M_sun` within ~3e-6 and return native slow/rapid
node-identified f/p1, interior-match residual <1e-5 and interface
traction residual <1e-9. This is **local** mass-branch selection only,
not a general relativistic radial-stability theorem.

| EOS | R1.4 [km] | slow f [Hz] | rapid f [Hz] | slow p1 [Hz] | rapid p1 [Hz] |
|---|---:|---:|---:|---:|---:|
| same-FSU2 Maxwell | 13.9403775 | 1916.9739 | 1916.9636 | 5129.2852 | 5128.2010 |
| same-FSU2 CSS | 13.9403735 | 1916.9975 | 1916.9671 | 5131.4545 | 5128.2691 |
| official FSU2 corrected CSS | 13.8705946 | 1928.4415 | 1928.4188 | 5184.7191 | 5182.2348 |
| official TM1e corrected CSS | 13.1641337 | 2084.1801 | 2084.1792 | 6498.4524 | 6498.0556 |
| official DDME2 corrected CSS | 13.2620260 | 2059.2318 | 2059.2265 | 6290.3527 | 6288.3180 |
| corrected BBP+BPS + FSU2 CSS | 13.8321635 | 1944.5813 | 1944.1648 | 5504.7217 | 5460.5861 |
| corrected BBP+BPS + TM1e CSS | 13.0575108 | 2078.6108 | 2078.5867 | 6580.2473 | 6571.8809 |
| corrected BBP+BPS + DDME2 CSS | 13.1470351 | 2064.8359 | 2064.5217 | 6609.5733 | 6496.6135 |
| corrected BBP+BPS + GM1 CSS | 13.7599553 | 1941.5760 | 1941.0632 | 6161.4447 | 6071.2481 |
| corrected BBP+BPS + PK1 CSS | 14.3022788 | 1861.1508 | 1860.7897 | 5443.9870 | 5394.1080 |

## 8. Independent numerical accuracy / gaps

`validated_convergence/RESULT.json`: tests native input grids
2049,4097,8193 samples on 5 cases. All four modes verified at each
resolution. For example, same-FSU2 Maxwell 4097->8193:
`|delta f|≈0.0047 Hz`, `|delta p1|≈0.063 Hz`.
Other phase-interface EOSs show `|delta p1|≈0.0001–0.094 Hz`;
no asymptotic convergence-order proof or error bound.
For same-FSU2 Maxwell original-source full DOP853 Python orbital mode
reference at `n_c=0.30` with native 8193:
slow/rapid `delta f≈+0.00293/+0.00325 Hz`,
`delta p1≈+0.02360/+0.02432 Hz`. Numerical reference is
independent integration but uses the **same physical EOS contract**.
Differences smaller than the observed sampling scale must not be
promoted to resolved astrophysical mode shifts.

Causal/pressure/chemical-potential/one-sided `cs²` checks are local
to admitted cold EOS contracts. Frozen-composition perturbations,
elasticity, superfluidity, surface tension, full-GR metric
perturbations/quasinormal damping, microscopic phase competition and
global radial stability remain outside this claim.

## 9. Status, runnable reproduction and follow-on

Core success: **SCOPED_IMPLEMENTATION_VERIFIED** for phasewise native TOV and
Cowling, Maxwell vs CSS separate physical constructions, tested source inputs,
and stated numerical tolerances. Full source-authentication = **PARTIAL**
due mismatch with official provider ZIP SHA; original archive bytes
preserved, original thermo independently checked. Source-microscopic
unification = **NOT REQUIRED, NOT CLAIMED**. Independent BBP primary
digital table = **UNRESOLVED**, secondary numeric table handled explicitly.

Commands from the project root (all outputs recorded in
`data/cowling_maxwell_native_20261119/`):
```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend
# Fetch and VERIFY exact observed source bytes; retain provider checksum conflict.
python3 scripts/fetch_cowling_m03_sources_20261119.py
python3 scripts/survey_same_rmf_original_unified_20261119.py
python3 scripts/build_original_fsu2_same_rmf_three_layer_20261119.py
python3 scripts/run_official_original_unified_controls_20261119.py
python3 scripts/run_original_unified_css_sharp_20261119.py
python3 scripts/run_bbp_bps_matched_core_20261119.py
python3 scripts/validate_maxwell_css_native_20261119.py
python3 scripts/run_fixed_mass_maxwell_css_20261119.py
pytest -q python_frontend/tests/test_cowling_maxwell_native_20261119.py \
 python_frontend/tests/test_cowling_maxwell_sources_20261119.py
```
These scripts require original downloaded CompOSE ZIP members and the exact
secondary BBP PDF, whose locations and SHA identities are given in the
respective output JSON. External sources are **not** fabricated if missing.

**Next candidate stage:** `COWLING-EOS-SOURCE-CALIBRATION-04`: obtain and
authenticate original BBP digitized data; resolve official ZIP checksum
conflicts with the publisher, compare Q3 Euler first-law repairs and
provided compositions, refine input discretization around sharp interfaces;
extend physical mode diagnostics to independent fluid composition closures
before publication. Do not add new gate inflation or equate BBP's secondary
roundoff-repaired EOS to the author's microscopic model.

### Actual executed regression check

Executed on the connected Ubuntu host (GNU Fortran 13.3, Python 3.12),
with OpenBLAS and OpenMP limited to one thread. The focused, inherited
and source-integrity regression suite using
`test_cowling_maxwell_native_20261119.py`,
`test_cowling_maxwell_sources_20261119.py`, earlier smooth Cowling/TOV and
independent polytrope tests, and the Maxwell source checks finished with
**99 passed, 3 skipped, exit code 0**. The raw log is
`data/cowling_maxwell_native_20261119/verification/TEST_RUN.log`.
Both modified Fortran modules compile with
`-O0 -fcheck=all -Wall -Wextra -ffree-line-length-none`. The unused
argument warnings in the DVODE dummy-Jacobian callback (MF=10/22) do not
affect the numerical results. The source rehydration readback of the 3
provider archived ZIPs and the secondary BBP PDF was also actually run
with pinned observed hashes. Full-source ZIP SHA conflict remains flagged.

**Interpretation:** A root with all four selected mode labels at one
source grid is implementation-verified under the listed numerical controls.
It is not a claim of unlimited frequency precision or independent mode-family
tracking through an avoided crossing. The older 0.30 fm^-3 BBP/DDME2 and
BBP/GM1 slow p1 mode identity ambiguity is retained, despite the new valid
1.4 M_sun snapshots; these are different equilibrium stars and are not
contradictory numerical evidence.
