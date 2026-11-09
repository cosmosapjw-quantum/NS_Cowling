# Mass-dependent crust-core matching response: research plan

Date: 2026-11-09  
Track: historical smooth crust-core matching cohort; physically separate from native sharp-Maxwell branch  
Status: implementation staged on a separate research branch; exact-head pilot CI is the next acceptance gate

## 1. Scientific objective

Determine how thermodynamically admitted crust-core matching prescriptions alter neutron-star structure and the relativistic Cowling f- and p1-mode spectrum as a function of stellar mass along the initial stable TOV branch.

The central question is:

> How does thermodynamically admissible crust-core matching uncertainty propagate through neutron-star structure and non-radial oscillation spectra across the stellar mass sequence, and how does that response depend on core and crust families?

The intended contribution is not another single-mass comparison. It is a controlled, paired, crossed-cohort map of mass-dependent response and its core/crust dependence, with explicit numerical and mode-identity gates.

## 2. Frozen input authority and scope

The new execution branch is research/mass-dependent-cowling-v1, forked from live native-Maxwell feature head a09ec30d74b363652d2d1aefe990557c468679dd. PR #3 remains open, draft and unmerged; this track does not modify it.

The frozen smooth-cohort authority is:
- oscillation_release_v1/data/spectrum_summary_234.csv (234 rows; target mass 1.3300106263 solar masses);
- oscillation_release_v1/data/eos_projected/ (projected tables);
- oscillation_release_v1/provenance/release_summary.json;
- docs/V5_PAIRED_MATCHING_SENSITIVITY_20261109.md;
- docs/EXISTING_SPECTRUM_REUSE_VERIFICATION.json.

The crossed cohort is 13 crusts × 3 admitted core families × 6 smooth matching prescriptions: 39 fixed core/crust pairs and 234 method rows. The three cores are BBB_BHF2, GM_GM1 and RG_SLY4. The six methods are css, cubic_F, cubic_P, maxwell, quintic_F and quintic_P. These are not 234 independent microscopic EOS draws.

The historical method token maxwell is a smooth near-pressure-plateau bridge. It has no density discontinuity or latent heat and is not a first-order Maxwell phase transition. The native sharp-Maxwell solver is a separate physical track with distinct interface transfer conditions and conversion assumptions.

The prior spectrum-reuse audit reports 234 rows and recomputed_modes=0. All old mode frequencies are fixed at approximately 1.3300106 solar masses. They are a frozen reference anchor, not mass-dependent evidence.

## 3. Hypotheses and falsifiers

H1: low-mass p1 amplification. The paired six-method p1 range may increase toward low masses and contract at higher masses. Canullán-Pascual et al. (2025) already report low-mass sensitivity, so reproducing that qualitative behaviour is an external-consistency check, not novelty by itself.

H2: core-dependent response shape. The mass where the six-method range peaks, the width of the high-sensitivity interval, or the sign pattern of method deviations may depend on the core family. This is a candidate differentiator, not an assumed result. It is falsified if core-conditioned curves collapse within numerical uncertainty.

H3: crust-conditioned response. At fixed core, the method response may depend on crust family and may be concentrated in particular core/crust pairs rather than follow a universal mass-only curve. It is falsified if the response is consistent with a common mass-only curve.

H4: mode contrast. The p1-mode may remain more matching-sensitive than the f-mode across a substantial part of the stable branch. It is falsified if the contrast disappears after mode-identity and convergence gates or is confined to a small subset.

H5: physical mechanism. A later analysis may connect response to crustal acoustic travel time and mode amplitude near the matching region. No sensitivity-kernel result is assumed. A kernel claim requires a variational derivation for the relativistic Cowling system, including boundary and matching terms.

## 4. Confirmatory mass design and stable branch

The predeclared target masses are
\[
M/M_\odot = 0.8,\ 1.0,\ 1.2,\ 1.33,\ 1.4,\ 1.6,\ 1.8,\ 2.0.
\]

For every EOS:
1. Scan a logarithmic central-pressure grid below the projected EOS pressure ceiling, anchored to the frozen 1.33-solar-mass central pressure.
2. Integrate TOV at every sample. Preserve failures explicitly.
3. Extract only the initial connected rising-M(Pc) branch. The current implementation resolves a first maximum after two consecutive fractional mass decreases larger than 1e-4. A failed sample is a hard barrier; later points are not joined across it.
4. If the scan terminates before a maximum is resolved, report a branch endpoint, not Mmax.
5. For each mass bracketed inside the branch, pass the adjacent sampled central pressures to the production Fortran target-mass solve.
6. Do not extrapolate a target beyond the resolved branch or across a failed-sample barrier.

The turning-point criterion is a barotropic branch/stability proxy. Unqualified radial stability claims require a separate radial-pulsation calculation, especially for genuine phase transitions or disconnected branches.

## 5. Spectrum solve and mode identity

The frozen f and p1 frequencies seed search windows using only the mean-density scaling
\[
\nu_{\rm dyn}\propto\sqrt{M/R^3}.
\]
The predeclared window is 0.72–1.28 times the scaled reference frequency. This is a root-bracketing device, not a fitted physical prediction.

The `cell_solve` EOS-cell integrator now accepts a configurable `r_limit_m` in its namelist; its historical default remains 1e6 m. For this mass scan, the explicit exploratory ceiling is 1e8 m and is recorded in every model manifest. The wrapper adds `solve_cowling_tov_star(..., task=star)` so sequence scans use the same pressure-cell integration and central-start contract as the mode calculation. The generic `case_driver real_table_star` API is retained separately. The Python wrapper accepts optional explicit pc_bracket_pa=(Pc_lo,Pc_hi). When omitted, its former 0.65–1.45 multiplicative bracket remains unchanged. When supplied, the bracket is validated and passed into native cell_solve. This constrains the target-mass search to the sampled rising branch.

A spectrum row is admitted only if:
- the target-mass residual closes;
- both shooting residuals satisfy the declared tolerance;
- the smooth l=2 Cowling W/V node signature is established in both profile segments:
  - f: zero W and V internal nodes in inner and outer segments;
  - p1: one W and V internal node in the inner segment, zero in the outer segment.

The identity contract is smooth_l2_WV_nodes_v1. It is not transferable to the sharp-Maxwell interface branch.

## 6. Paired estimands

For y in {f,p1}, at fixed core c, crust s, mass M and method m, define
\[
\Delta_m^{(y)}(M;c,s)
=100\left[\frac{y_m(M;c,s)}{\bar y_6(M;c,s)}-1\right],
\qquad
\bar y_6=\frac{1}{6}\sum_{m=1}^{6}y_m.
\]

The six-method paired range is
\[
{\cal R}_6^{(y)}(M;c,s)
=100\,\frac{\max_m y_m-\min_m y_m}{\bar y_6}.
\]

An estimand enters a paired summary only if all six methods for that fixed core/crust/mass cell have valid, identified roots. Missing methods are not imputed. At every mass report the complete-pair denominator, median, 90th percentile and maximum range, plus core/crust-resolved maps and method deviations.

The 39 fixed core/crust pairs are the physical pairing units. Do not treat the 234 method rows as independent EOS samples, posterior samples or population-level draws. The existing 1.33-solar-mass variance decomposition remains descriptive for its selected finite cohort, not causal or posterior.

## 7. Numerical validation gates

N1: frozen baseline reproduction. Recompute TOV M,R at every reference central pressure using the same `cell_solve task=star` EOS-cell integrator as the spectrum calculation, with the recorded v5 settings `r0_m=10 m`, `rtol=1e-9`, `step_cap=0.05` and the identical first-row finite surface-pressure cutoff. The separate `case_driver real_table_star` path is not an acceptable direct reproduction comparator because its ODE and central-start contracts differ. Current pilot acceptance requires relative M,R differences no larger than 2e-5.

N2: target-mass solve. Require relative mass residual no larger than five times the configured mass_rtol; default mass_rtol is 1e-10.

N3: shooting and identity. Require both matching residuals no larger than 1e-6, plus the exact W/V node signature above. A small determinant alone is insufficient.

N4: branch resolution. Preserve all pressure samples and failures; state whether the first turning point was resolved. A truncated endpoint is never relabelled Mmax.

N5: refinement versus effect size. For representative low-mass tail, 1.33-solar-mass anchor and high-mass cases, vary TOV/Cowling tolerance, maximum step and root brackets independently. Frequency movement must be at most one tenth of the smallest paired method effect reported. Otherwise weaken or withhold the corresponding physical claim.

N6: independent benchmark. Compare representative EOS/mass points with published relativistic Cowling/full-GR results or a separately maintained solver. Cowling/full-GR systematic offsets are distinct from within-Cowling matching shifts and must be tabulated separately.

N7: recognized matching prescriptions. The six historical smooth bridges are not a sufficient substitute for direct literature baselines. Add (i) simple EOS concatenation as a deliberately non-thermodynamic diagnostic control only, (ii) linear P(epsilon) matching, (iii) the thermodynamically consistent P(epsilon) interpolation used by Canullán-Pascual et al. (2025), reconstructed from the published prescription, and (iv) a CUTER-style or otherwise documented thermodynamically consistent crust reconstruction. For every method preserve endpoint data, units, source checksums, chemical-potential/first-law closure, causality and interpolation receipts. Do not label simple concatenation physically admissible, and do not call an approximate CUTER-like implementation CUTER itself.

## 8. Structural and mechanism extensions

This increment produces TOV M(Pc), branch endpoints/turning-point candidates, R(M), and f(M),p1(M). It does not yet produce crust mass, crust thickness, tidal Love number/Lambda, full-GR complex mode frequencies or damping times.

Next mechanism work:
- extract the transition radius and enclosed mass from pressure-event output or validated profile interpolation;
- compute crust acoustic travel time with the correct relativistic Cowling radial measure and EOS sound speed;
- compare normalized eigenfunctions and overlap near the matching region;
- derive a variational/sensitivity kernel only after all boundary and matching terms are included and checked.

Next astrophysical work:
- implement Love-number/tidal deformability and source-backed observational likelihoods;
- combine NICER mass-radius posteriors, maximum-mass and tidal constraints;
- admit CompOSE/source-backed EOSs only with source checksums, units, phase labels and thermodynamic admission receipts.

## 9. Publication decision

JKPS / EPJA: viable only after the full crossed cohort has mass-dependent paired spectra, explicit stable-branch coverage, refinement below one tenth of the smallest reported effect, and independent benchmarks. The distinguishing result must be a reproducible core/crust-conditioned mass-response structure, not simply “p1 is more sensitive than f”.

A&A / MNRAS: requires realistic EOS provenance, structural observables, tidal deformability and observational inference.

PRD: requires a new relativistic perturbation/sensitivity derivation or systematic Cowling/full-GR response result.

## 10. Implementation status

Stable-branch extraction, adjacent-Pc bracketing, six-method complete-pair estimands, explicit Pc and radius-limit controls, the shared cell_solve TOV-only backend, the all-cohort runner and the pilot/full workflow are staged on the new branch. The first integration run (Actions 42598463512) exposed a non-comparable baseline: the runner used `case_driver real_table_star` (`r0=1 m`, `rtol=1e-8`) while the v5 reference metadata records `cell_solve` (`r0=10 m`, `rtol=1e-9`, `step_cap=0.05`). The same run also reached the old hard-coded 1e6-m direct-cell radius ceiling. The wrapper now uses the same `cell_solve task=star` backend and exposes the ceiling without changing its default. The validator syntax defect is corrected. These changes await a new exact-head CI run; the failed run yielded zero valid mass-dependent spectra and no Mmax result.

References:
1. Canullán-Pascual et al. (2025), Astronomische Nachrichten, DOI 10.1002/asna.20240150, arXiv:2502.02373.
2. F. Nola (2026), A&A 713, A34, DOI 10.1051/0004-6361/202659810.
3. Kunjipurayil et al. (2022), Phys. Rev. D 106, 063005, arXiv:2205.02081.
4. Pradhan and Chatterjee (2022), Phys. Rev. C 106, 015805, DOI 10.1103/PhysRevC.106.015805.
