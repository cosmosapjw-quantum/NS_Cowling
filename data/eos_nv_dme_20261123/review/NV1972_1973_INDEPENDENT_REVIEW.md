# Independent NV1972/NV1973 source and implementation review

## Scope and identities

This review independently inspected the rendered original NV1972 Table I (journal p1487), Eq. (3.10)–(3.11) (p1478), and NV1973 pp307–311 and p317. NV1972 PDF SHA-256: `fd981fe5c9d7d914329c533fb57a2ac9375ed3d641f76ecc01b307bf82cd500b`; DOI `10.1103/PhysRevC.5.1472`. NV1973 DOI: `10.1016/0375-9474(73)90349-7`.

`COEFFICIENT_VISUAL_REVIEW.json` records an exact-decimal comparison against an independent transcription from the rendered original. It checks all 64 B coefficients, 56 C coefficients, and 28 independent symmetric-D coefficients. There are zero mismatches. The comparison script is `independent_coefficient_check.py`. This is an input identity check, not proof of the validity of the DME physical approximation.

## Functional variation

Set x=rho_n, y=rho_p, T_q=(3/5)(3*pi^2)^(2/3)*rho_q^(5/3), B_n=B(x,y), B_p=B(y,x), and analogously C_n,C_p. The nuclear functional, excluding rest masses and Coulomb, is

```
H = t0*(tau_n+tau_p)+V(x,y)
  + B_n*(tau_n-T_n)+B_p*(tau_p-T_p)
  + C_n*|grad x|^2+C_p*|grad y|^2+D*grad x.grad y
  + (W/2)*grad(x+2y).J_p.
```

B and C are not symmetric. D is symmetric and its cross-gradient term occurs once. NV1973 changes V to its own Table 1 while retaining the original NV1972 B,C,D. The common bare kinetic coefficient convention and rest-mass constants must be separately recorded when modern constants are used.

Holding tau and density gradients fixed, the local neutron derivative is

```
V_x + B_n,x*(tau_n-T_n)+B_p,x*(tau_p-T_p)-B_n*T'_n
+ C_n,x*|grad x|^2+C_p,x*|grad y|^2+D_x*grad x.grad y.
```

The full neutron gradient contribution, after subtracting the divergence of its gradient derivative, is

```
- C_n,x*|grad x|^2+(C_p,x-D_y)*|grad y|^2
- 2*C_n,y*grad x.grad y-2*C_n*laplacian(x)-D*laplacian(y).
```

The proton expression follows by the corresponding explicit species exchange. The spin-orbit contributions to the local potentials are `-W/2 div(J_p)` and `-W div(J_p)` respectively. Omitting the neutron LS splitting does not remove the neutron density derivative of the proton-only spin-orbit energy.

The kinetic operator has coefficient `a_q=t0+B_q`. In uniform matter, even when `tau_q=T_q`, the local potential contains `-B_q*T'_q=-B_q*k_Fq^2`. It cancels the B contribution to the Fermi energy, recovering the derivative of `t0*T+V`. This is a useful independent test of rearrangement terms.

The code `compose_eos_nv_dme.f90` as first reviewed has the correct coefficient shapes and transposition, analytically differentiated Horner expressions, species derivative order, A derivatives, and gradient/spin-orbit variational derivatives. A second density derivative can diverge on a density axis where B remains nonzero; the code explicitly marks the A Hessian undefined there. First derivatives remain finite.

## Radial kinetic operator and boundary conditions

For the continuum operator `-div(a(r)*grad)` and reduced radial function u=r*R,

```
h_l u = -(a*u')' + [a'/r+a*l*(l+1)/r^2+U]u.
```

Omitting `a'/r` is incorrect. A conservative radial R formulation with spherical cell-volume weights includes it implicitly and can be preferable to explicitly discretizing the reduced-u form.

The paper specifies even parity R(Rcell)=0 and odd parity R'(Rcell)=0. For u these become u(Rcell)=0 and u'(Rcell)=u(Rcell)/Rcell. The latter is not a homogeneous Neumann condition on u. For a free constant-mass spectrum the reference roots are `j_l(k*Rcell)=0` for even l and `j_l'(k*Rcell)=0` for odd l.

In the proton radial Hamiltonian, the LS contribution is

```
(W/2)*(rho_n'+2*rho_p')/r * [j*(j+1)-l*(l+1)-3/4].
```

The neutron channels instead carry degeneracy `2*(2*l+1)` with no LS split.

An independent Wolfram Language evaluation returned exact zero for (i) the difference between the symbolic spherical Euler–Lagrange derivative and the neutron gradient expression above, (ii) the radial reduction identity, and (iii) the Robin transformation. The tool emitted an undefined-global-symbol warning for the assigned intermediate symbol `L`, but returned the complete correct expression and all three exact zeros; this warning is retained in the tool record rather than hidden.

## Table 3 conventions

N is the neutron number, Z the proton number, and A=N+Z. The first cell is N=140,Z=40,A=180. This follows both the table header and the filling prescription on p310. In Fig. 2, the unusual isotope mass labels denote the total baryon count in the entire cell, including free neutrons.

Table 3 orbital chemical potentials exclude rest masses. The electron kinetic formula on p307 subtracts electron rest mass. The beta variation on p310 restores `m_e+m_p-m_n`. Therefore, if the implementation uses the total relativistic electron energy and total chemical potential,

```
E/A-m_n = (E_nuclear+E_Coulomb+E_e,total+Z*(m_p-m_n))/A
beta residual = mu_p,orb-mu_n,orb+mu_e,total+(m_p-m_n).
```

The electron Coulomb contribution must use the same neutral-cell potential convention as the proton potential. It is the volume-average electron potential with the opposite charge. Net charge neutrality makes the combined beta condition invariant under a potential gauge shift.

The 1973 Table 3 row at n_b=5.77e-3 fm^-3 prints E/A-m_n=1.996 MeV. This edition-specific value must not be silently changed to the later 1.966 value. Table 4 is a fit to these discrete ground-state calculations, not an independent self-consistent solver.

## Boundary smoothing and scientific admission

NV1973 p310 replaces rho and tau within a width 1/k_F of the cell edge by their averages over the rest of the uniform gas region on every iteration. The precise numerical definition of that gas region is not supplied. This is a documented author numerical prescription, but it is not identical to solving the unmodified finite-cell variational problem. An implementation should record the chosen averaging region, width, normalization action, and whether smoothing is enabled. It should compare interior profiles and energy with unsmoothed or enlarged-cell calculations before claiming the boundary correction is negligible.

A fixed published (N,Z,n_b) HF solution establishes a reconstructed cell, subject to convergence and grid tests. It does not establish a newly minimized ground-state EOS. A finite catalog search can establish only the lowest accepted energy among the explicitly searched cells and initial guesses. The paper itself emphasizes multiple dissimilar initial guesses because multiple local minima exist.

An EOS built from these calculations needs a separately identified energy-to-pressure representation. Its pressure must obey `P=n_b^2 d(E/A)/dn_b` on the declared branch/interpolant. Changes of minimizing composition and possible loss of convexity need explicit treatment. Failure of a bounded search, failure to converge, failure of a discrete first-law check, or lack of resolution must remain different statuses. A matched table that passes matching constraints does not repair an unresolved HF source.

## Additional source-integrability finding

The successful transcription audit does **not** make the printed coefficient set internally consistent. An independent visual reading of Eqs. (3.8), (3.11), and the construction prescription on p1488 gives, with `s=(x+y)/2`, `L=V_D^L`, and `U=V_D^U`,

```
F(x,y)=x*L(x)+y*U(s)-B(x,y)/4
C-B_x/4=-2L-xL'-(y/2)U'(s)
D-[B_y(x,y)+B_y(y,x)]/4=-4U(s)-sU'(s).
```

Writing `U(s)=sum u_k*s^k`, for k=i+j and j>=1,

```
C_ij-(i+1)*B_(i+1,j)/4 = -k*binomial(k-1,i)*u_k/2^k,
D_ij-[(j+1)*B_(i,j+1)+(i+1)*B_(j,i+1)]/4
    = -(k+4)*binomial(k,i)*u_k/2^k.
```

The exact-rational last-printed-digit interval audit in `TABLE_INTEGRABILITY_INTERVALS.json` finds compatible intersections at k=0,1,2, but incompatible ones at k=3,4,5,6. Every tested coefficient relation above total degree 6 contains zero within its rounding interval. The inconsistencies are therefore localized and far larger than printed rounding. The discussion of higher-order corrections in Sec. VI places those corrections into V_NM; it does not resolve these B/C/D polynomial identities.

These identities do not uniquely identify the author-intended corrections. `EXACT_DEGREE3_AMBIGUITY.json` contains two exact-rational counterexamples with all other degree-3 entries inside their original printed rounding intervals: changing only D03 to -662699.325, or changing only B04 to +143292.335, makes all five degree-3 identities exact. This proves that an apparent D03 exponent correction is conditional on retaining B, even though that correction is typographically plausible.

An exploratory floating-point linear-program completion, with exact source-decimal intervals as inputs, found 2,3,8,24 smallest-support repair families at degrees 3,4,5,6 respectively, with 1,2,3,5 changed printed entries. These are numerical ambiguity witnesses, not an authenticated erratum. The complete enumerated supports are in `SPARSE_REPAIR_AMBIGUITY.json`. Keeping every B entry and repairing six C and five independent D entries is one additional assumption and must receive a distinct model identity.

For a high-frequency uniform-matter variation, the C/D matrix alone is insufficient to infer instability. Normalized occupied plane waves multiplied by `1+eta_q*cos(Q*x)`, with Q above all occupied momentum differences, provide an explicitly orthonormal variational family. Its leading Q^2 contribution is controlled by

```
M = [[C_n+(t0+B_n)/(4*rho_n), D/2],
     [D/2, C_p+(t0+B_p)/(4*rho_p)]].
```

A negative eigenvalue gives a constructive high-frequency energy-lowering direction of the truncated functional. Such a statement concerns that functional; it does not establish a physical no-go theorem outside the slow-variation regime of the DME.

## Initial review status (superseded by final review below)

- Original coefficient transcription: **implementation-verified within the 148 printed entries**.
- Continuum functional variation and radial transformation: **derived and symbolically checked**.
- New DME local module first review: **implements the printed table and its variational derivatives consistently; the printed B/C/D table has independently established internal coefficient-identity conflicts**.
- HF end-to-end implementation, nonlinear convergence, original boundary prescription, finite composition search, and EOS thermodynamic admission: **pending actual implementation/results review**.


## Final read-only scientific and code review, 2026-11-23

The final reviewed HF source SHA-256 is
`69c3ecd5ab5e045c4fa0157759701e68b4897a20bb8fdb903764faaef8f72946`;
the DME source is
`27477ef266ee3a9a67e716daebf1e0404510099ad67c3eeca25e44fd2c8ef9f0`.
The structured final receipt is `FINAL_INDEPENDENT_HF_REVIEW.json`.
Within this review scope, no remaining code blocker was found for publishing
these explicitly bounded research results. The original source and both
potential EOS admissions remain blocked for the scientific reasons below.

The initial implementation findings were checked again in the actual final
source. Native success now requires both `status=ok` and convergence; later
field failures clear the convergence flag. Density, kinetic density and proton
spin-current residuals are measured, and a final raw-field solve verifies the
fixed point. This final solve refreshes the spectral flag, occupied/addition
frontier levels and beta interval. The frontend requires finite scientific
success fields, complete profile length, matching configuration identities,
unchanged executable/source bytes and zero process exit. It always records
`source_eos_admitted=false`. Empty plan IDs and unknown seed configurations
are rejected. `HF_INITIAL_CODE_REVIEW.json` remains a historical finding
record; its pending statements are superseded by this final receipt.

### Verification actually read

The independent reviewer read the final native logs with 43 free-spectrum,
normalization, adjoint and electron checks, plus 24 orbital variation and
pressure-dilation checks: 67 in total. The free-spectrum error ratio for
240 to 480 radial cells is 0.25000343881. The maximum displayed scaled orbital
variation error is 3.994084744634e-9. These checks ran with array/runtime checks
and invalid, zero and overflow traps. The focused frontend/legacy regression
log contains 140 passing tests in 31.49 s; the 11 new protocol cases are part of
that scope and must not be counted a second time.

The enlarged-basis test uses N=130, Z=50, n_b=2.79e-4 fm^-3 and mesh 160 in the
D03x100 diagnostic. Increasing (l_max, radial states) from (16,11) to (18,13)
changes E/A by 4.2355e-12 MeV and pressure by 2.0761e-9 MeV fm^-3. This is useful
evidence for that cell. The `spectral_complete` flag is still a cutoff sentinel,
not a proof that all omitted channels are harmless for every density.

### Actual nonlinear campaigns and numerical limits

The 28-job D03x100 campaign contains 22 published fixed-composition cell/seed
combinations and six A=180 even-charge neighbors. It has 23 `ok` results and
five iteration failures. All eleven clustered published-cell starts converged;
six uniform starts converged and five did not. The failed starts are preserved.
Across successful outputs the maximum raw replay residual is 1.32274e-6,
particle-count error 4.54747e-13, and absolute dilation-pressure discrepancy
2.05389e-10 MeV fm^-3. All 23 declared successful profile row counts were checked.
These are fixed-cell stationary points, not ground-state EOS rows.

The literal printed-table campaign contains 22 attempts: one converged uniform
first cell, twelve source-domain step failures, eight iteration failures, and
one out-of-domain raw-density failure. The converged uniform result proves that
the cluster failures alone cannot justify an absence-of-solutions claim.

For N=140,Z=40,n_b=2.79e-4 fm^-3 in the diagnostic model, meshes 160,320,480 give
E/A-m_n = -1.5059714251, -1.4575748282, -1.4487528817 MeV, respectively.
Pressures are 3.7025864923e-4, 3.6106337269e-4, 3.5941293230e-4 MeV fm^-3.
The material mesh bias cannot be hidden behind a tight nonlinear residual.
A diffuse start agrees with the clustered energy at mesh 320. The refinement
trend is consistent with second order at this one cell, but no EOS-wide
continuum error bound has been established.

At fixed A=180 and n_b=2.79e-4 fm^-3, mesh 320, the tested charges are
{40,42,44,46,48,49,50,51,52}. Z=50 has the lowest accepted energy in that set.
Relaxed neighboring Z=49,50,51 energies are -1.5616116789, -1.5634759227 and
-1.5352944240 MeV per baryon. The Z=50 first-order beta interval includes zero.
This verifies a restricted charge comparison; A, cell size, other shapes,
other densities and unsearched starts have not been minimized globally.

### Final executable replay and source provenance

The 28-job campaign was compiled from HF SHA `8b8739d5...` and DME SHA
`0d267566...`, with binary SHA
`a0b3ae94260aee443d4e34abb194be2deef399ce2a6b3f394c20f7a963949700`.
Its seven original build files were retained. Independent file comparison
confirmed that the subsequent DME change only added recursive declarations;
the HF change refreshed the post-replay status/frontier reports. The energy,
pressure and iteration equations did not change. A launch-time source hash is
not a compilation manifest, and these two identities are kept distinct in
`../hf/BUILD_PROVENANCE.json`.

An initial final-binary frontend smoke run used mixing .03 and at most 600
iterations. It correctly returned nonconvergence: its spin residual
2.201360979e-6 exceeds the required 2e-6, even though its energy differs from the
converged value by only 2.56e-13 MeV per baryon. This failure remains evidence,
not a result to relabel or erase. Repeating the original campaign configuration
(.08 mixing, maximum 1200 iterations) with the final executable converged in
225 iterations, exit zero, unchanged input identities. Its E and P decimal
strings match the historical campaign. The raw replay residual is 8.24458e-7.
The refreshed beta interval is [-7.318522981,-3.745206954] MeV and still excludes
zero. The final executable hash is
`8a2cf49d3eaed986403ce1402ec6206940d6f7a36bda2f752adda54b9c0f7dbf`.
The apparent discrepancy was an iteration-budget/configuration difference;
no tolerance was relaxed to obtain a success label.

### Source and EOS admission decision

The exact transcription audit and source-identity contradiction are compatible
findings: the code transcribes the printed table faithfully, while the table
fails its own generating identities. D03x100 is a conditional diagnostic and
does not repair the remaining C/D identity conflicts or authenticate an author
erratum. Numerical convergence of that modified functional cannot establish
source identity. Raw unsmoothed HF also differs from the paper's partly
unspecified boundary averaging prescription.

A 32-source-identity by 3-crust-contract by 7-method census records 672
`BLOCKED_SOURCE` cells, with zero native matcher calls and zero new survivors.
The 32 input identities comprise 19 unique core model names. The three crust
categories contribute 224 blocked cells each: missing Appendix D, diagnostic
fixed cells not admitted as a ground-state EOS, and original-table internal
inconsistency. This is an admission census, not 672 physical matching failures.

The release can preserve these source failures, derived ambiguity witnesses,
Fortran implementations and bounded numerical successes. Admission of a new
original or diagnostic ground-state EOS still requires an authenticated source
contract, sufficient resolution/basis control, composition and metastability
search, and an energy/pressure representation with its first-law and convexity
checks. Discrete coordinate-dilation agreement does not by itself satisfy those
EOS-wide requirements. No new matched EOS is certified by this review.

Public supporting evidence is under `../hf/`, `../validation/`,
`../source_contract/` and `../remaining_census/`. Full raw execution receipts are
retained separately using logical archived-runtime paths `hf_runtime/...` and
`root_runtime/...`; no user device or backup location is exposed here.
