# COWLING-EOS-SOURCE-CALIBRATION-05: original GPPVA core source reconciliation

Date: 2026-11-19 KST. Implementation branch:
`research/cowling-eos-source-calibration-05-20261119`.
Parent from 04: `f191abe843d81c6a1e5797dd9edc9a1558f095fa`.
Current result code and machine-readable SSOT:
`data/cowling_rmf_reconciliation_20261119/SOURCE_RECONCILIATION_OUTCOME.json`.

This is a bounded, **source-preserving theoretical/numerical audit**
with actual Fortran DVODE TOV/Cowling calculations. Source equivalence,
finite Maxwell matching, CSS effective three-layer stellar construction
and experimental 1.4-M_sun mode identity remain **distinct claims**.
Unpublished original provider RMF table-generation code was not obtained.
No parameters or energy offsets were fit to force agreement with the
provider. No source ZIP, PDF, `eos.nb`, `eos.thermo` or `eos.compo`
bytes were edited.

## 1. Theory, original source contracts and numerical conventions

Use the signature `(-,+,+,+)`, keep `c,\hbar` in the physical
discussion. Original cold, neutral npe(\mu) matter with conserved
baryon density `n_b` [fm^-3], total pressure and energy density
`P,epsilon` [MeV/fm^3], chemical potential `mu` [MeV]. For
barotropic matter at zero temperature and lepton number unconstrained:

```text
 P = n_b^2 d(epsilon/n_b)/dn_b ,
 mu_B = (epsilon + P)/n_b ,
 n_p = n_e + n_mu ,   mu_n = mu_p + mu_e.
```

The original CompOSE Q1 through Q7 definition is not inferred.
See the *CompOSE Reference Manual*:
https://link.springer.com/article/10.1140/epja/s10050-022-00847-y

```text
 Q1=P/n_b;  Q3=mu_B/m_n-1;  Q4=mu_q/m_n;  Q5=mu_l/m_n;
 Q6=f/(n_b*m_n)-1; Q7=e/(n_b*m_n)-1.
 mu_e = mu_le-mu_q = m_n*(Q5-Q4).
```

For the three **original** unified FSU2/TM1e/DDME2 tables,
`Q5==0` in the examined homogeneous phase and `e=f`.
For populated zero-temperature free electrons, the independently
calculable chemical potential is

```text
 n_e = n_b * Y_e,
 mu_e,Fermi = sqrt[(m_e*c^2)^2 + (hbar*c)^2*(3*pi^2*n_e)^(2/3)].
```

A source `mu_Q3` or `mu_Q4` consistent with its original
thermodynamic output, its original lepton fraction and these
assumptions must satisfy both algebraic identities. This requires
no local nucleonic RMF implementation, and can be tested directly
against the **original** Zenodo-matched CompOSE bytes.

Three observed ZIPs were independently byte-verified in stage 04:
https://zenodo.org/records/14808903 (FSU2),
https://zenodo.org/records/14808958 (TM1e),
https://zenodo.org/records/14508967 (DDME2).
The separate publisher `eos.zip_checksum.txt` entries still
differ from the independently byte-identical Zenodo records.
Preserve both facts rather than treating them as contradictory.

## 2. Same-density original core versus local published RMF

No B-spline, root-fit or input-density interpolation was used.
Evaluate the existing local RMF solver at each provider original
homogeneous `Iphase=6` knot with `n_b in [0.08,0.45] fm^-3`.
The original provider and local beta-equilibrium compositions
match to O(10^-6), but this does not demonstrate identical
bulk EOS or identical field configurations.

| Quantity | FSU2 | TM1e | DDME2 |
|:--|--:|--:|--:|
| Examined homogeneous rows | 75 | 63 | 75 |
| max |Y_{p,e,mu}^{local}-Y_{p,e,mu}^{provider}| | 5.92e-6 | 5.72e-6 | 7.28e-6 |
| max relative pressure mismatch | 4.2213% | 15.5182% | 0.005156% |
| median relative pressure mismatch | 1.5799% | 1.2288% | 0.004385% |
| original max |Q3 mu_B - (epsilon+P)/n_B| | 6.2715 MeV | 8.5449 MeV | 3.718e-5 MeV |
| original median relative P-vs-n²d(epsilon/n)/dn mismatch* | 1.035% | 0.675% | 1.82e-6 |

*The pressure derivative is a **sampled cubic-spline diagnostic**
on the interior density knots, not an analytic error theorem.
The corresponding independent local RMF relative residuals
are approx 1.69e-8, 1.53e-8 and 2.59e-8 respectively.
Source finite digit precision, spline choice and boundary
extrapolation are not suppressed or silently reinterpreted.

Published nuclear-matter saturation properties are checked
independently from the local RMF Lagrangian using symmetric
`Yp=0.5` matter and `K=9 dP/dn|n0`:

| Model | Local n0 [fm^-3] | local E/A-m_n [MeV] | local K [MeV] | published K [MeV] |
|:--|--:|--:|--:|--:|
| FSU2 | 0.15040186 | -16.26488 | 237.685 | 238 |
| TM1e | 0.14522294 | -16.26312 | 281.165 | 281 |
| DDME2 | 0.15194724 | -16.14174 | 250.902 | 251 |

Source parameter references: FSU2 (Chen & Piekarewicz 2014,
parameters tabulated by Tolos et al. 2017),
https://diposit.ub.edu/dspace/bitstream/2445/116038/1/669453.pdf .
TM1e (Shen, Ji, Hu, Sumiyoshi 2020),
https://user.numazu-ct.ac.jp/~sumi/eos/table4/guide_EOS4.pdf .
The local RMF theory is not grossly missing its saturation
point. An ad-hoc constant shift of baryon rest energy
`epsilon->epsilon+n_b*Delta_m` leaves the cold pressure
exactly unchanged, so it cannot explain this pressure gap.

## 3. Independent lepton/CompOSE Q4 inconsistency (source-only evidence)

Read **original** phase-6 `Y_e,Q4,Q5` from the three
byte-authenticated GPPVA tables. Verify
`mu_e^source=m_n(Q5-Q4)` against the free-electron
Fermi energy derived from the original source's own `n_e`.
Numerical findings:

| Model | max |Delta mu_e| [MeV] | Q4 rows compatible with printed Ye rounding |
|:--|--:|--:|
| FSU2 | 2.3693 | 0/75 |
| TM1e | 1.9775 | 0/63 |
| DDME2 | 0.00514 | 64/75 |

The source `Y_e` is given to about five decimal places.
The rounding compatibility comparison directly calculates the
electron Fermi-energy interval corresponding to
`Y_e+-5e-6`, with an additional small finite-`Q4`
table precision allowance. Under the explicit free
electron/zero-`Q5` assumptions, the FSU2 and TM1e
discrepancies cannot be attributed to displayed `Y_e`
rounding. This is a **source thermodynamic metadata
consistency issue** independent of the local nucleonic
RMF model; the precise generating pipeline cause is
still unresolved. DDME2 is much closer, but even its
full table is not claimed exactly source-identical.

Evidence: `charge_chemical_potentials/RESULT.json`;
source-level algebraic implementation:
`python_frontend/compose_eos/rmf_source_reconciliation.py`.

## 4. Discriminating independent checks of the pressure mismatch

1. Hold the provider's printed `Y_p` fixed and recompute the
   baryon mean fields plus exactly neutral free leptons at
   its original density. Compare the resulting pressure with
   the local beta-root solution. The *largest* pressure
   perturbation due to the printed-composition difference:
   FSU2 `2.2678e-4`, TM1e `1.2748e-4 MeV/fm^3`.
   The **smallest** provider-versus-local P discrepancy
   among those controls is `0.05143`, `0.06627`
   respectively. The maximum *fraction-roundoff effect
   divided by actual source gap* is <`0.00075`.
   Therefore beta-root rounding does not account for the
   measured original-vs-local pressure discrepancy.
2. Force muons absent (an explicit **hypothesis control**).
   The original source has populated muons in the high
   homogeneous phase. Removing them changes `Y_p`
   at `n~0.45 fm^-3` by roughly `-0.0347` (FSU2)
   and `-0.0252` (TM1e), far above the actual
   few-ppm source fraction mismatch. This is not an
   admissible source-equivalent convention.
3. For TM1e, replace only the printed isoscalar RMF
   couplings/masses by the more precise published
   values in the historical TM1 table, retaining
   the TM1e isovector couplings and mixed quartic
   interaction. This is NOT a fit to the source.
   At n~=0.446 fm^-3 the existing absolute pressure
   error is `0.206568 MeV/fm^3`; the high-precision
   variant produces `0.207332 MeV/fm^3`, i.e.
   rounding precision does not resolve the discrepancy.
4. Numerical derivative of original `epsilon/n_B`
   shows the bulk thermodynamic output itself is
   systematically less self-consistent for FSU2/TM1e
   than for DDME2. First-law pressure-primary correction
   and charge-chemical-potential metadata correction
   are **different operations**; neither is automatically
   equivalent to regenerating the original nuclear RMF
   field solution or correcting a publisher's source code.

These controls disfavor an isoscalar pressure error due
solely to electron fraction rounding, missing muons,
or the extra published decimals in TM1e isoscalar
coefficients. There is not enough evidence to identify
a unique coding term or original generator version.

### 4.1 Pressure-gap partition without modifying or fitting original EOS

The observed same-density total pressure discrepancy has an **exact**
algebraic partition (for a common interpolated density variable)

```text
ΔP = P_local - P_provider
   = n_b² d[(epsilon_local-epsilon_provider)/n_b]/dn_b
     + R_local - R_provider ,
R_i = P_i - n_b² d(epsilon_i/n_b)/dn_b .
```

We evaluated the derivative using the same cubic spline procedure for
both source and local baryon energy-per-particle arrays, excluding two
boundary knots from the RMS comparison. No normalization, fitted EOS
coefficients or changes to any original `eos.thermo` input occur.
The following ratios are **RMS amplitudes divided by RMS(ΔP)**, not
signed contributions that must sum to one:

| Model | RMS ΔP [MeV/fm³] | energy-functional term / RMS ΔP | source identity residual / RMS ΔP | local RMF identity residual / RMS ΔP |
|:--|--:|--:|--:|--:|
| FSU2 | 0.1916540 | 1.28373 | 0.34542 | 4.85×10^-5 |
| TM1e | 0.1223170 | 0.66269 | 0.50481 | 3.23×10^-5 |
| DDME2 | 0.00098665 | 1.00093 | 0.03787 | 0.01042 |

Independent code/tests reproduce the algebraic balance to maximum
absolute residual below `7e-12 MeV/fm³` for all three original-source
rows. This is a **numerical identity validation**, not an inference that
any specific publisher code function was incorrect. In particular,
the source-vs-local energy functions and the source first-law residual
both materially contribute to FSU2/TM1e mismatch, and they partially
cancel in some density ranges. The local RMF energy and pressure are
nearly consistent under this sampled derivative; this says nothing
about whether its chosen couplings are precisely the provider's.

Evidence under `EXACT_DENSITY_AUDIT.json` field
`models.<MODEL>.pressure_difference_decomposition`, including the
source-density-resolved four pressure terms; the non-trivial analytic
polytrope and source-evidence tests are in
`test_cowling_rmf_reconciliation_{20261119,evidence_20261119}.py`.

## 5. Same-crust paired native CSS/Maxwell stellar runs

Construct two three-layer EOSs from identical source:
`surface`, `outer`, `lower_bridge`, and ALL actual
original inner crust first-law corrected `n,P,epsilon`
samples including their endpoints are numerically
**exactly identical**. Replace only upper homogeneous
core and apply the explicitly named native matching
model. We used the original phase-local CSS source
from stage04 and the newly generated same-named
RMF core from stage03. All source identities preserved.

At target `1.4 M_sun` with locally ascending M(n_c):

| Shared crust, CSS | Official source-core R1.4 [km] | Newly computed-core R1.4 [km] | Difference [km] | Change in slow f [Hz] |
|:--|--:|--:|--:|--:|
| FSU2 | 13.84915677 | 13.94037473 | +0.09121796 | -12.616 |
| TM1e | 13.12670682 | 13.20075233 | +0.07404551 | -9.387 |
| DDME2 | 13.23083414 | 13.23093987 | +0.00010573 | -0.0227 |

The FSU2 slow `p1` is `~100.924 Hz` lower
for the generated-core CSS than for original-core CSS
at the same mass. TM1e original-core slow W1 mode
has the independently documented near-interface
V-node classification ambiguity; never fabricate its
strict slow-p1 counterpart. Both source branches'
rapid p1 can be compared numerically but the
family-continuity caveat remains.

All three paired *CSS* stars pass native Fortran TOV,
both slow/rapid f/p1 modes in the generated-core model,
interface traction and node gates. The original
phase-local TM1e source-core case retains its
strict-node W1/V2 ambiguity. Because we held
all crust source points fixed, the differential
result is more directly attributable to the
**different upper core source output** than comparisons
made with different crust families.

Generated FSU2 homogeneous core + same source
inner/outer crust admits a **new finite Maxwell**
at approximately:
`P_t=0.18656666 MeV/fm^3`,
`n_low=0.04856440`,
`n_high=0.04907366 fm^-3`, followed by actual
1.4-M_sun TOV and slow/rapid f,p1.
Generated TM1e and DDME2 cores do **not**
pass finite Maxwell on this phase-local original
inner crust. Most crucially, the *official
provider original homogeneous core* did not
admit a Maxwell crossing under the stage04
source phase-local finite-support tests, and
**must not be equated to the generated core**.
Do not use CSS/Maxwell interchangeably or
call a phenomenological model a microscopic
three-phase unification.

## 6. Unexpected original eos.mr reuse: not independent TOV evidence

Original provider full *unified* ZIPs were compared with
their corresponding **crust-only** GPPVA source packages.
The `eos.thermo` contents and EOS density counts
are different, but the **eos.mr bytes are identical**
in all three source pairings, including identical SHA-256
on the actual ZIP-extracted file:

| Model | Unified original density knots | Crust-only density knots | Same `eos.mr`? |
|:--|--:|--:|:--|
| FSU2 | 426 | 269 | Yes |
| TM1e | 455 | 303 | Yes |
| DDME2 | 415 | 288 | Yes |

As an additional warning, the source full-unified
FSU2 `eos.mr` maximum mass is
`2.06965 M_sun`, while its public *unified*
model metadata reports `2.37 M_sun`:
https://compose.obspm.fr/eos/215 .
We cannot establish whether this is an intentional
legacy crust-core continuation or a packaging error
without the generating scripts. It is *not an
independently recomputed TOV run on the published
unified source* merely because it appears within
that ZIP. Never use this `eos.mr` file as the
independent numerical TOV oracle for this task.
Source identity was **independently** verified
against Zenodo in stage04; byte identity does
not establish identity of the physical provenance
of an optional derived data product.

Evidence: `provider_mr/RESULT.json`;
source schema:
https://link.springer.com/article/10.1140/epja/s10050-022-00847-y .

## 7. Outcomes and stop conditions

**Scoped source and implementation confidence**
- Original FSU2/TM1e/DDME2 archive byte identity:
  VERIFIED IN STAGE04 with Zenodo and source SHA.
- Exact-density beta-equilibrated nuclear fractions:
  NUMERICALLY CHECKED to a few parts in 1e-6.
- Original charge chemical metadata Q4 vs
  independent ideal electrons:
  QUANTIFIED EXPLICIT INCONSISTENCY for FSU2/TM1e.
- Official full source EOS vs local published
  same-name RMF source equivalence:
  **FSU2/TM1e NOT ESTABLISHED / DISFAVOURED**;
  DDME2 HIGH NUMERICAL AGREEMENT, not equality.
- Same-crust controlled 1.4-Msun CSS TOV and
  Cowling: IMPLEMENTATION VERIFIED, with
  matched original crust n,P,epsilon source points.
- Generated FSU2 Maxwell finite crossing:
  NUMERICALLY VERIFIED in the stated model only;
  original source-core Maxwell: UNRESOLVED/NOT ADMITTED.
- Publisher-provided full-unified `eos.mr` as
  an independent source-original TOV reference:
  **BLOCKED**, reused crust-only bytes.
- Original GPPVA homogeneous field equations,
  original output Q1/Q3/Q4 conversion code:
  NOT OBTAINED, exact culprit still OPEN.
- BBP primary digital source: remains separately
  UNRESOLVED, not part of this narrow RMF correction.
- Microscopic unification of all layers:
  NOT REQUIRED and NOT CLAIMED for our explicitly
  constructed three-layer phenomenological
  TOV/Cowling models.

The smallest next constructive science operation,
if authorized, is a direct original-provider field/chemical
potential data or generator source comparison at one
exact `n_B` to resolve the Q1/Q3/Q4 outputs,
not another round of arbitrary EOS fitting or
a repeated terminal audit.

## 8. Reproduction, code and provenance

Previous source parent:
`f191abe843d81c6a1e5797dd9edc9a1558f095fa`.
Precise source raw ZIPs and 35-row BBP secondary PDF
are not republished in the Git commit; use
`scripts/fetch_cowling_m03_sources_20261119.py`
to rematerialize with exact source hashes.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend
python3 scripts/audit_rmf_core_provenance_20261119.py
python3 scripts/audit_compose_charge_thermo_20261119.py
python3 scripts/audit_fixed_fraction_rmf_20261119.py
python3 scripts/run_paired_core_source_stars_20261119.py
python3 scripts/audit_compose_mr_provenance_20261119.py
python3 scripts/build_rmf_source_reconciliation_outcome_20261119.py
pytest -q python_frontend/tests/test_cowling_rmf_reconciliation_20261119.py \
          python_frontend/tests/test_cowling_rmf_reconciliation_evidence_20261119.py
```

Detailed evidence under
`data/cowling_rmf_reconciliation_20261119/`:
`EXACT_DENSITY_AUDIT.json`,
`charge_chemical_potentials/RESULT.json`,
`fixed_composition/RESULT.json`,
`paired_core/RESULT.json`,
`provider_mr/RESULT.json`, and
`SOURCE_RECONCILIATION_OUTCOME.json`.
Each new source reading is explicitly coupled
to the exact original source `eos.thermo` SHA
already verified in stage04.

### Fresh verification receipt (actually executed)

After an initial environmental failure caused by the separately isolated
stage-05 worktree lacking the stage-03 source-locked BBP secondary PDF,
we restored that PDF from the previous worktree with unchanged SHA-256
`2f12c3c6208709ace642e30b484cb56d1b86b7de48bd6540598165a1f8f3803d`.
The failing source-intake test passed again without weakening its assertion.
The **full** inherited and new stage-05 regression suite was then rerun after the bounded pressure-gap decomposition extension:

```text
python3 scripts/verify_rmf_reconciliation05_20261119.py
144 passed, 3 skipped in 43.58s (exit=0)
```

The final fresh regression included two new analytic pressure-split tests,
one source-evidence assertion, and one aggregate-result assertion; their
red/green checks are recorded in the current session. The new pressure
partition does **not** identify a unique original-generator code defect.

The Python compileall pass and GNU Fortran 13 build with
`-O0 -fcheck=all -Wall -Wextra -ffree-line-length-none` both returned
exit status zero, with only unused dummy-Jacobian callback warnings.
All original `eos.thermo` SHA-256 values match the official, independently
verified stage-04 source receipt; all tested code/test input SHA-256
values were stable before and after the run. This is bounded
**SCOPED_IMPLEMENTATION_SOURCE_LOCKED**, not an EOS physics theorem.
The complete machine-readable evidence, command, current base SHA and logs
are in `data/cowling_rmf_reconciliation_20261119/verification/
VERIFICATION_RECEIPT.json` and `TEST_RUN.log`.

**Stage-05 disposition:** `SCOPED_SOURCE_RECONCILIATION_DIAGNOSED`.
The original-generating FSU2/TM1e source functional has **not** been
reconstructed to source identity, and no parameter fitting was used;
source-original Maxwell publication claim remains UNRESOLVED.
The lower layers are verified unchanged in actual generated-core
paired stellar runs and previous raw-source data is preserved.

Do not claim a merged PR or altered original source
if no such operation was requested/performed.
