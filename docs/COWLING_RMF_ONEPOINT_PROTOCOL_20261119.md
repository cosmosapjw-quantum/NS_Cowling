# COWLING-RMF-ONEPOINT-06 | Frozen research protocol

Prepared before the new one-density numerical comparison (2026-11-19 KST).
Parent: `6050cf11cd7f090fcdeb82f59039366cdff62678`.
Existing Stage-05 source EOS identity: original GPPVA(FSU2/TM1e/DDME2)
CompOSE tables, copied to an isolated working tree with exact inherited SHA.
No original files, physics couplings, or starting EOS archives are modified.

## Research question and falsifiable hypotheses

**Question:** Can an original provider's exact single-density homogeneous RMF
mean-field quantity, component breakdown, or generator source be obtained and
compared with local FSU2/TM1e fields? If the GPPVA original generator is
not public, how far can source-backed independent reference data take us
without falsely promoting same-named RMF to identical provider physics?

- **H01 (source-provenance, confirmatory)**: The official Pais GPPVA(2023)
  CompOSE source packages expose P/epsilon/Y and no individual sigma,
  omega, rho or generator. *Falsifier:* inspect official `eos.micro` or
  repository archive and find the raw fields, term breakdown or generator.
- **H02 (independent original-theory author data, confirmatory)**:
  The independent Shen-Ji-Hu-Sumiyoshi 2020 EOS4 original TM1e T=0
  fixed-composition uniform-matter table can be compared at its **exact**
  published grid point (log10 rho_B = 14.5 [g cm^-3], Y_p=0.10), with no
  density/proton interpolation. Evaluate the identical T=0, fixed Yp
  baryons-only RMF state with the code's published TM1e parameters.
  Compare published M*, pressure, F (per baryon rel M=938 MeV),
  mu_n and mu_p. *Falsifier:* disagreement beyond table-roundoff and
  numerical tolerance; preserve failed comparisons.
  **Critical:** Shen is the original TM1e functional's independent
  author, NOT the Pais GPPVA CompOSE table-generation code.
- **H03 (exact source GPPVA, confirmatory):** At one closest-to-0.2
  original phase-6 density knot separately for FSU2 and TM1e (and DDME2
  control), original provided P/epsilon/Yp/Ye/Ymu can be compared with
  native beta-state values, while the local meson fields, contributions
  and EOM residuals are recorded as LOCAL ONLY. *Falsifier:* nonfinite
  fields, violated Euler/beta/charge/meson residual, or source identity drift.
- **H04 (local decomposition, confirmatory):** The local epsilon and
  P split exactly into neutron+proton Fermi gases, electron+muon gases,
  sigma meson, omega and rho fields, vector quartic and mixed omega-rho
  contributions. Check EOM with fixed model couplings and recomposition.
  *Falsifier:* any composition sum fails relative 1e-9, meson field EOM
  >1e-6 fm^-3, or energy/pressure difference above 1e-8 MeV fm^-3.

## Methods and explicitly chosen sample

- **E01, original-author research:** Search official Pais/CompOSE GPPVA
  datasheets and source bundles, original Grill/Pais/Providencia papers,
  and publicly available author code repositories. Record URL, checksum,
  status, source authority and limitations. Do not infer absence of
  private data from public web search.
- **E02, primary Shen TM1e:** Download immutable official 2020
  `eos4.t00.gz` from
  `https://user.numazu-ct.ac.jp/~sumi/eos/table4/eos4.t00.gz`.
  Existing local download SHA-256:
  `346b2701f7372817c516c3c945fdee0b3da5af9f08e2fe83f743951d1b5ebaef`.
  Parse 16 fixed columns per Table VI of
  `https://user.numazu-ct.ac.jp/~sumi/eos/table4/guide_EOS4.pdf`,
  selecting log10 rho=14.5 and Yp=.10 with exact source file values,
  physical uniform core checked by heavy nucleus fraction = 0.
  Source provides baryons ONLY. Do not add leptons or compare raw Shen
  baryon-only pressure to GPPVA total pressure.
- **E03, exact GPPVA one-point:** Select *nearest original tabulated*
  `n_B` to 0.2 fm^-3, in the original homogeneous phase 6 and within
  [0.08,0.45], no source interpolation or source repair. Evaluate
  `gppva_rmf.beta_state()` in parallel as an independent calculated
  state, no fit. DDME2 control separately from native Fortran routine;
  lack of direct DDME2 meson field values remains explicit.
- **E04, component audit:** For FSU2 and TM1e, compute field potentials
  using the independent local stationary `baryons()` and
  `free_fermion()` values, both canonical and on-shell expression:
  `epsilon_vector = -W + g_w omega n_B + g_r rho n3/2`,
  `P_vector = W`, with every potential divided by (hbar c)^3.
  Decompose into sigma quadratic/cubic/quartic, omega quadratic/quartic,
  rho quadratic, omega-rho mixing; verify Euler and meson field
  identities, local composition, zero T. Do not call this original
  provider's pressure component output.
- **E05, output and verification:** Machine-readable stage-specific
  source index, summary and comparisons; unit tests with known polynomial
  identities and real one-point source provenance. Run scoped native
  regressions and keep true PASS/FAIL states. Publish only derived outputs,
  source URLs/checksums and code. Never publish or overwrite original
  full EOS source table or original Shen archive.

## Stop conditions / claim gate

Stop once E01-E05 are run and independently verified, or on a concrete
network/source access blocker. Original Pais micro-field identity requires
an actual original author field row or code; an independently verified
Shen uniform TM1e reference supports only *theoretical functional*
agreement. Missing Pais source remains UNRESOLVED and may prompt a
concrete minimal author data request, not an inferred identical origin.

Explicit exclusions: no parameter fitting, no new microscopic three-layer
unification claim, no misleading Maxwell source-coexistence claim, no new
TOV/Cowling reruns unless exact source change warrants it.
