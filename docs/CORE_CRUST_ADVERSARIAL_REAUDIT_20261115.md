# Core–crust adversarial re-audit, 2026-11-15

Baseline: `d90d3dca033f759c79ed034b3b8ff939acb43145`, PR #15,
`remediation/2026-11-15-m5a-maxwell-eos-contracts`. The earlier handoff at
`6b28fadc...` was superseded before this audit: the baseline already had the
independent report-admission layer and an 80-attempt census.

## Scientific correction

The most consequential defect was a physical energy shift mislabelled as
normalization. The old matcher added `dm=mn_core-mn_crust` to crust F and mu,
and `n_B*dm` to energy density. For GPPVA_FSU2H this changed every retained
free-energy row by +0.56533 MeV per baryon. The new implementation preserves
physical F, epsilon, n, P and mu_Q, with Q columns re-encoded at the output
reference. The existing disclosed Q3 Euler re-derivation remains distinct.
See `MASS_NORMALIZATION_CORRECTION_20261115.md` and the corrected current
`data/remediation/MASS_NORMALIZATION_CONTRACT.json`.

A uniform physical energy shift changes
`(c_s/c)^2=(dP/dn)/(d epsilon/dn)` to
`(dP/dn)/(d epsilon/dn+dm)`; reversibility does not imply physical identity.
Five native matching methods failed a newly added same-physics/different-header
regression before repair and pass afterwards. No source table was modified.

## Other corrected findings

| Object | Counterexample | Correction |
|---|---|---|
| Native Maxwell root | Exactly represented bisection midpoint is discarded | Return `fmid==0` immediately |
| Python Maxwell root | Two causal-source roots at 1049.99 and1050.01 MeV hidden from800-point scan | Union-knot quartic pressure isolation using cubic density stationary partitions |
| Multiple Python roots | Implicit first-root selection could overstate a complete branch | Default requires explicit choice; metadata records all roots and selected-interface scope |
| Maxwell output | Nonzero phase labels could not be written without a sidecar | Supported phase-only eos.compo, no invented abundance masks; complete one-point axes |
| Native cold input | NaN T/Yq or nonzero Q5 with trusted beta flag accepted | Finite metadata and existing beta bound enforced |
| Admission | Corrupt outer rows, reordered/empty bridge, malformed speeds, contradictory extrema could pass | All-row finite/Euler checks, branch secants, exact indices and schema, structured rejection |
| Pair runner | Direct repeated invocation deleted work; unknown method ran CSS | Refuse nonempty output and unknown methods; preserve environment failures |
| Persisted verification | stdout PASS ignored failed/missing detailed gates | Check saved schema, n, overall, failed checks, C1/C2/C3/C4/C5/C6/C10/C11 applicability |

No physical tolerance was relaxed. Admission secants are necessary row checks;
they do not independently prove a higher-order interpolant's continuous
causality or differential first law.

## Executed verification

- Local generator suite: **91 passed**,7.35s (Python3.12, SciPy1.17, NumPy2.3.5).
- Native build: `bash fortran_eos/build.sh build/compose_eos_fortran`, exit0.
- Native frontend+pipeline: **357 passed,9 skipped**,28.69s.
- Focused new native/root/admission/mass regressions: **29 passed**,2.14s.
- Readback regression checks physical arrays, Q1–Q7, phase/status/presence
  masks and C3/C4 applicability. Synthetic Maxwell now passes the actual census
  path; this is not evidence for a real astrophysical Maxwell survivor.
- The initial 5-knot exact-root fixture failed C10/C11 because every derivative
  point was masked (large log spacing plus the phase jump). The identical
  analytic EOS was sampled on129 dyadic knots, retaining the exact midpoint
  root. No verifier bound changed.
- Corrected-code exploratory census:80 configured,60 evaluated,20 incomplete
  source attempts;10 CSS survivors,6 distinct source combinations;0 Maxwell.
  The committed-code census receipt is recorded separately. Old output values
  must not be reused merely because the count happens to agree.
- New CI steps run the configured census and preserve its complete work tree,
  successful and failed pair reports, numerical tables and hashes as an artifact.
  Exact final CI evidence belongs to the live PR/run and final handoff receipt.

Native execution used an isolated `/tmp/NS_Cowling_audit_20261115_d90d3dc`
clone on the authorized device. The existing Dropbox checkout was observed at
`9179edda...` with ongoing uncommitted generator/Fortran/DAG edits and was not
changed. Reconcile those changes deliberately before a local agent resumes.

## Remaining research obligations

1. Sharp Maxwell mass-pair admission remains unchanged:30 real candidate
   attempts stop at the convention gate. This is not a no-coexistence result.
   A provider-backed common physical reference is needed; do not copy the
   removed additive energy shift into that path.
2. Sampled branchwise first-law/causality acceptance in Python remains. Root
   isolation is now interval based; that does not complete continuous
   thermodynamic validation. APR's high-density acausality is not repaired.
3. Two crust directories lack eos.nb/eos.t/eos.yq;20 attempts are not evaluated.
   RG cores fail their beta contract; do not override it.
4. Reconstructed SPG_M1/M2/M4 sources are byte-identical; raw Iphase7 convexity
   and reconstruction provenance remain unresolved. The six distinct survivors
   are implementation-admitted composites, not certified microscopic crust models.
5. Phase-only export is not composition-aware Maxwell v2. Abundance-preserving
   native matching, generated CSS→phase-segmented Cowling adaptation, elastic
   crust and solid-fluid perturbation conditions remain open.
6. Missing raw M0 causes source-dependent skips; full spatial TF/pasta and the
   complete390-case scientific campaign remain outside this bounded repair.

## Primary definitions and independent checks

CompOSE reference manual **v3.01 (2022-04-01)**, despite URL filename v3.00:
https://compose.obspm.fr/download/pdf/manual_v3.00.pdf,
§4.1.2 equations4.5–4.9 and§4.2.2. Q3 is mu_b/mn−1; Q4 is mu_q/mn;
Q5 is mu_l/mn. For Il=0, Euler enthalpy also includes Yq*mu_q. For Il=1,
local neutrality means Yq=Ye+Ymu, not generally Yq=0. The existing Yq0 route
is a scope restriction, not an independent neutrality proof.

Alford–Han–Prakash, PRD88,083013(2013):
https://arxiv.org/html/1302.4732v3. CSS first-law identities were independently
checked symbolically with Wolfram: dP/dmu−n, epsilon+P−mu*n and dP/depsilon−s
all vanish. The Seidov threshold is a small-core onset test, not whole-branch
stability. Exa/Parallel/SciSpace assisted primary-source discovery.

Verdict: **PARTIALLY-COMPLETE** for the scientific project; bounded corrections
are implemented and executed. Keep PR15 draft/open/unmerged.

## Committed-code census receipt

`data/receipts/core_crust_audit_20261115/CENSUS_b11a530.json` records a clean
run of code commit `b11a5305c461872f43a0a17611fd8cea796ff1d9`, all80 pair
report hashes, source identities and10 survivor output hashes. Counts remain
80/60/20/10/6 as above. Numerical physics changes are material: for
SPG_M1+BBB_BHF2, the CSS latent energy falls from4.853739572 to
1.775592014 MeV fm^-3 and s=(c_s/c)^2 changes from0.0873574303 to
0.0855269965. For SPG_M1+GM_GM1 latent energy changes from0.434594362 to
0.807502668 MeV fm^-3. Equal survivor counts do not mean equivalent EOS.

This receipt uses native Python3.12.3, not the GitHub Actions Python3.13
environment. Final CI repeats the full census and uploads its complete evidence.
