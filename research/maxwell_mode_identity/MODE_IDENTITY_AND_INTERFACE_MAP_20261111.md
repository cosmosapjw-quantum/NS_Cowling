# Sharp-Maxwell Cowling spectra: interface maps and mode-identity protocol

**Date:** 2026-11-11  
**Status:** DERIVED_INTERFACE_MAPS + SYNTHETIC_HELPER_TESTED  
**Physical labels:** UNRESOLVED  
**EOS scope:** the existing analytic two-phase fixture only; no source-backed CompOSE pair is admitted by this work.

## 1. Research question and scope

The C7X crust-composition D4 blocker prevents source-backed EOS pair admission. It does not prevent a separate methodological extension of the already implemented analytic two-phase Cowling boundary-value problem. The open issue on that fixture is no longer simply finding roots: the existing native/reference comparison found two positive-frequency roots per conversion prescription, but the root sets were sorted and compared as sets. Their f-, p1- or discontinuity-g-mode identities have not been established.

The source implementation already stores phase-resolved W,V profiles and one-sided interface values. The Python reference already counts sampled W and V sign changes separately in each phase. This increment therefore adds profile normalization, interface diagnostics, phase-aligned overlaps and conservative branch assignment rather than duplicating the node counter or re-integrating the ODEs.

Inspected provenance:
- source-authority branch: research/maxwell-native-fortran-v1
- source-authority head: a09ec30d74b363652d2d1aefe990557c468679dd
- Python phase-aware solver blob: 8959f90a37aba63191227bb363e97302988864e9
- native Cowling blob: ce4fac829d3a5c211139d5174541a7053793afc6
- native spectrum module blob: 726a03b219827e77412175c68e3d7c28ef32ef84

The existing analytic-fixture roots recorded at that source head are:
| Conversion | Root 1 | Root 2 |
|---|---:|---:|
| Slow | 1764.8989380477 Hz | 4937.6449322607 Hz |
| Rapid | 4937.0162338857 Hz | 19493.6326018255 Hz |

Existing convergence/reference checks apply to those roots on that fixture only. They do not establish family identity.

## 2. Interface convention and exact maps

Use the static spherical metric

\[
ds^2=-e^{2\nu(r)}c^2dt^2+e^{2\lambda(r)}dr^2+r^2d\Omega^2
\]

and the solver's polar displacement variables

\[
\xi^r=e^{-\lambda}\frac{W}{r^2}Y_{\ell m}e^{i\omega t},\qquad
\xi^\theta=-\frac{V}{r^2}\partial_\theta Y_{\ell m}e^{i\omega t},
\]

with the corresponding azimuthal component. Define

\[
A=\left(\frac{\omega r_t}{c}\right)^2e^{\lambda-2\nu},\qquad
Q=AV+\nu'W,\qquad
\Delta P=-h e^{-\lambda}\frac{Q}{r_t^2},\qquad h=\epsilon+P.
\]

Here A is dimensionless, nu' has units m^-1, W has units m^3 and V has units m^2. Inner denotes the high-density side; outer denotes the low-density side. Pressure and metric potentials are continuous at the interface, while h and density may jump.

Continuity of Lagrangian pressure requires

\[
h_oQ_o=h_iQ_i.
\]

Set \(\alpha=h_i/h_o\), so \(Q_o=\alpha Q_i\). At cold Maxwell coexistence, \(h=\epsilon+P=\mu_Bn_B\); equality of \(\mu_B\) on the two sides therefore gives \(\alpha=n_{B,i}/n_{B,o}\).

### 2.1 Slow/material conversion

The material interface has continuous normal displacement, \(W_o=W_i\). Solving pressure continuity for V_o gives

\[
\begin{pmatrix}W_o\\V_o\end{pmatrix}
=
M_s\begin{pmatrix}W_i\\V_i\end{pmatrix},\qquad
M_s=
\begin{pmatrix}
1&0\\
(\alpha-1)\nu'/A&\alpha
\end{pmatrix}.
\]

### 2.2 Rapid/equilibrated conversion

For rapid conversion, continuity of the isobaric displacement gives

\[
W_o-\frac{Q_o}{\nu'}=W_i-\frac{Q_i}{\nu'}.
\]

Together with pressure continuity this implies \(V_o=V_i\), and

\[
\begin{pmatrix}W_o\\V_o\end{pmatrix}
=
M_r\begin{pmatrix}W_i\\V_i\end{pmatrix},\qquad
M_r=
\begin{pmatrix}
\alpha&(\alpha-1)A/\nu'\\
0&1
\end{pmatrix}.
\]

Both maps satisfy det(M)=alpha, pressure-jump residual zero, and reduce to the identity at alpha=1. They are not interchangeable at finite density jump. Wolfram Language symbolic evaluation verified both pressure-continuity residuals, the rapid isobaric-displacement residual, both determinants and the zero-jump limit. The Python helpers implement these exact matrices and dimensionless residual checks.

These maps are derived under the stated sharp, cold, one-conserved-charge Maxwell assumptions. They do not supply a finite-rate phase-conversion law.

## 3. Why frequency ordering is not mode identity

A sharp density discontinuity can support a discontinuity gravity mode. Miniutti, Pons and Berti studied this mode in relativistic neutron stars and its association with density jumps (2003, MNRAS; [SciSpace record](https://scispace.com/papers/non-radial-oscillation-modes-as-a-probe-of-density-5gugy5xabr)).

Tonetto and Lugones derived different junction conditions for slow and rapid conversion. In their ideal sharp-interface treatment, the slow/material limit supports a finite-frequency discontinuity g-mode, whereas the rapid-equilibrium discontinuity-g frequency tends to zero because the phase boundary adjusts and the associated buoyancy restoring force is removed (2020, [arXiv:2003.01259](https://arxiv.org/abs/2003.01259)). Thus a positive rapid-interface-g root must not be invented by analogy with the slow branch.

A 2025 relativistic-Cowling study classifies discontinuity g-modes across strong transitions and twin-star categories ([SciSpace record](https://scispace.com/papers/classification-of-g-modes-for-neutron-stars-with-a-strong-pi816jxwynzb)). A 2026 full-GR study of sharp interfaces reports conversion-dependent f/p1 spectra and associates p1 sensitivity with tangential displacement behavior near the interface ([arXiv:2609.33002](https://arxiv.org/abs/2609.33002)). That full-GR result is motivation for recording tangential polarization, not a numerical prediction for the present Cowling model.

Avoided crossings create a second identity problem. Rodriguez et al. compared node/frequency nomenclature, generalized Cowling nomenclature and classification based on modal properties for protoneutron-star oscillations; modal-property tracking can retain physical behavior when simple ordered labels become inconsistent (2023, [arXiv:2304.00033](https://arxiv.org/abs/2304.00033)). This is methodological support from an adjacent problem, not direct evidence that the same crossings occur in this cold fixture.

The policy is therefore to store **spectral branch identity** separately from **instantaneous physical character**.

## 4. Profile diagnostics

### 4.1 Positive displacement norm

The displacement convention and redshifted kinetic-energy weighting give the positive profile norm

\[
N=\sum_{s\in\{i,o\}}\int_s h e^{\lambda-\nu}
\left(\frac{|W|^2}{r^2}+\ell(\ell+1)|V|^2\right)dr.
\]

The core and envelope fractions are \(F_i=N_i/N\) and \(F_o=N_o/N\). A shell integral around r_t reports interface localization. This is a normalization and continuation measure; it is not asserted to be a canonical orthogonality product between different stellar backgrounds.

### 4.2 Existing node count and one-sided interface values

The project Python reference already reports sampled W and V nodes separately in each phase. Preserve those outputs and test their mesh/tolerance convergence. Do not join the inner and outer W samples across r_t: a finite sign change between W_i(r_t) and W_o(r_t) is an interface jump, not an interior radial node. The rapid map can generate W_o != W_i by construction. V-node counts remain secondary diagnostics.

Each root record should retain (W_i,V_i,DeltaP_i) and (W_o,V_o,DeltaP_o), plus matching residuals.

### 4.3 Interface polarization

The physical tangential-to-radial displacement ratio on either side is

\[
\Pi_\pm=\sqrt{\ell(\ell+1)}\,r_t\frac{|V_\pm|}{|W_\pm|}.
\]

The factor r_t is required by the dimensions and displacement convention. If W is near its numerical floor, the ratio is flagged unreliable and raw one-sided amplitudes remain available. Pi alone cannot label a mode.

### 4.4 Phase-aligned profile similarity

To compare stars with different interface radii, align the two phases separately:

\[
x_i=r/r_t,\qquad x_o=(r-r_t)/(R-r_t),\qquad
\rho_s=h_se^{\lambda_s-\nu_s}\frac{dr_s}{dx},
\]

\[
\mathbf q_s=(W_s/r_s,\sqrt{\ell(\ell+1)}V_s).
\]

For adjacent backgrounds a,b define

\[
\mathcal O_{ab}=
\frac{\left|\sum_s\int_s dx\,\sqrt{\rho_{a,s}\rho_{b,s}}\,
\mathbf q_{a,s}^{*}\cdot\mathbf q_{b,s}\right|}
{\sqrt{N_aN_b}}.
\]

Report core and envelope overlaps as well as the global value. The global absolute value is taken only after adding the complex core and envelope overlap contributions. This removes one arbitrary global eigenfunction phase, but does not discard a relative sign/phase flip across the sharp interface. The overlap is symmetric under model exchange and invariant under global amplitude rescaling. It remains a chosen continuation similarity, not an observable or a unique inner product.

### 4.5 Branch assignment

Across adjacent models, a global Hungarian assignment minimizes a cost dominated by profile-overlap loss, with a weak logarithmic frequency penalty. If the chosen pair has low overlap, is not the local best overlap, or has a small row/column margin, the output is AMBIGUOUS_DO_NOT_FORCE_ID. The current helper tracks only continuous spectral branches. It keeps physical_character unresolved, because character can exchange at an avoided crossing even when a branch is continuously tracked. Complex profiles receive no binary real sign-reversal label.

## 5. Seed and continuation protocol

1. Set the density jump to zero and verify that the split-interface calculation converges to the smooth-EOS roots and profiles. Both matrices must tend to identity, but matrix convergence alone is insufficient.
2. Seed smooth-limit f and p1 solutions using converged regular eigenfunctions, phase-local W-node conventions and profile morphology jointly. Do not use frequency order alone.
3. Track slow and rapid conversion separately. A slow discontinuity-g branch may emerge from zero frequency as the jump increases; do not assume a finite positive rapid interface-g branch in the ideal rapid-equilibrium limit.
4. At each parameter step solve a sufficiently broad root window, evaluate profiles, and assign globally by overlap. Frequency is a tie-breaker only.
5. On low overlap or small assignment margin reduce the parameter step or track a multi-mode subspace; do not force labels.
6. Repeat forward and backward continuation. Away from flagged degeneracies, the branch mapping should agree.
7. Report character/localization diagnostics separately from branch IDs. Intermediate conversion rates require a finite-rate constitutive model and cannot be represented by an arbitrary interpolation between slow and rapid maps.

## 6. Validation state and acceptance gates

**Implemented and tested:** exact interface maps; pressure and kinematic residual diagnostics; positive profile norm; core/envelope trapping and interface-shell fractions; interface polarization; phase-aligned global/profile overlap; global branch assignment with ambiguity refusal. The existing phase-local node counter is reused rather than duplicated.

**Local validation:** 22 deterministic synthetic tests passed on 2026-11-11. Wolfram symbolic evaluation returned the exact identities listed above. Python source syntax compilation passed.

**Not done:** no production Fortran code changed or rebuilt; no actual four fixture-root profiles were fed through the helper; no zero-jump physical sequence or nontrivial jump continuation was executed; no f/p1/interface-g label was promoted; no realistic CompOSE EOS or source-backed TOV/Cowling run was used.

Acceptance before any physical label:
- zero-jump root/profile convergence against smooth-star f/p1;
- slow/rapid spectra continued independently over a nonzero jump sequence;
- profile and root convergence under half parameter step and doubled profile resolution;
- forward/reverse branch agreement away from ambiguity;
- native/Python phase-resolved profile comparison including one-sided interface values;
- explicit treatment of the slow interface-g branch and the zero-frequency rapid limit;
- physical character reported separately from spectral branch identity.

The outer profile endpoint in the helper is the current finite-pressure-cutoff background endpoint, not necessarily the exact P=0 surface. Comparisons must use a consistent cutoff policy.

## 7. Reproducibility and references

Run from repository root:

\`\`\`bash
PYTHONPATH=research python -m pytest -q research/maxwell_mode_identity/tests
\`\`\`

Dependencies are listed in requirements.txt. Machine-readable profile and status contracts are in profile_schema.json and MODE_IDENTITY_STATUS.json.

References:
1. G. Miniutti, J. A. Pons and E. Berti (2003), [SciSpace record](https://scispace.com/papers/non-radial-oscillation-modes-as-a-probe-of-density-5gugy5xabr).
2. L. Tonetto and G. Lugones (2020), [arXiv:2003.01259](https://arxiv.org/abs/2003.01259).
3. M. C. Rodriguez et al. (2023), [arXiv:2304.00033](https://arxiv.org/abs/2304.00033).
4. M. C. Rodriguez, J. C. Jiménez and I. F. Ranea-Sandoval (2025), [SciSpace record](https://scispace.com/papers/classification-of-g-modes-for-neutron-stars-with-a-strong-pi816jxwynzb).
5. Z. Z. Dong, H.-J. Kuan and N. Yunes (2026), [arXiv:2609.33002](https://arxiv.org/abs/2609.33002).

**Conclusion:** the method layer is ready for fixture-level continuation testing. The physical identities of the four existing roots remain unresolved.
