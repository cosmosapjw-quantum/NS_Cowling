# Core–crust re-audit continuation, 2026-11-15

Baseline: `b0a3f10e44488d3fc73313beb03f8e570288a2af` (PR #15, draft, unmerged).
Each item below followed counterexample → minimal correction → execution →
receipt. No physical tolerance was relaxed and no raw source was modified.

## 0. Dropbox checkout reconciliation

`/home/cosmosapjw/Dropbox/NS_Cowling` (local HEAD `9179edd`) has 23 dirty
tracked files. Every one equals an older committed version (`877229b`,
`f707572`, `bee7fd3`) or is an empty file; none contains unique uncommitted
content. It is a stale/rolled-back working copy, not separate work. It was left
unchanged. Receipt: `data/receipts/core_crust_audit_20261115/DROPBOX_RECONCILIATION.json`.

## A. Interval extrema instead of samples (implemented)

| Object | Counterexample (former code admitted it) | Correction |
|---|---|---|
| Python Maxwell branch (`_branchwise_first_law_report`) | 9-knot branch: sampled max c_s² = 1+5.0e-7, declared PCHIP max = 1+6.14e-6 (gate 1+1e-6); `_phase_arrays` admitted it | exact per-interval extrema: dn/dμ quadratic; c_s² ≤ 1+tol ⇔ g=(1+tol)μ dn/dμ − n ≥ 0 with g cubic (closed-form stationary points); reported extrema add quintic stationarity roots |
| same | derivative residual compared `ns(q)` with itself (vacuous) | removed; n=dP/dμ reported as holding by construction; knot P/ε reconstruction gates unchanged |
| TOV `find_causal_domain` | sampled max 1+5e-9 (gate 1+1e-8) vs interpolant 1+9.2e-5 | certified ordered branch-and-bound with exact bounds exp(min w)/max z′ ≤ c_s² ≤ exp(max w)/min z′ |
| same | certified ln P(ln ε) while the TOV RHS integrates the different PCHIP ln ε(ln P): 0.9995 vs 1+4.6e-4 | certification, cut point, ε and c_s² at the limit all use the RHS interpolant |

Tests cover narrow pockets, luminal tangency, constant speed on both sides of
the fixed gate, interior compressibility minima, a pocket in the interval
ending at the one-sided jump row, refusal to span a jump, and dense-grid
bounds. The native Maxwell kernel already decided exactly (quadratic density
and margin per Hermite interval); its 33-point loop only fills a summary.
Secant row checks in `scripts/match_report_admission.py` remain necessary row
checks, not continuous certificates.

## New defect: reconstructed crusts carry n_B in the pressure slot

`python_frontend/compose_eos/reconstruct_cld.py::emit_compose_set` writes
thermo rows `iT inb iYq n_B T Q3 Q4 Q5 Q6 Q7`; CompOSE expects
`iT inb iYq Q1 Q2 Q3 …`. The reader therefore uses Q1 = n_B and
**P = n_B²** for all 8 crusts in `download/reconstruction/` (e.g. 0.021025 =
0.145² MeV fm⁻³ at every SPG endpoint), 10²–10⁵ below n²dF/dn of the same
table. The CLD solver also leaves `P_cell = 0.0` ("filled by emitter"), and the
emitter never fills it. Q3 is encoded against 939.56536 while the header
says 939, and μ_p − M_P is written into the μ_q slot.

Consequence: the provider-Euler residual |μ_src − (F + P/n)|/μ is 1.5e-2–3.7e-2
for SPG and O(10) for GPPVA_DD2. The CSS matcher re-derives output μ from
F + P/n, so the defect vanished from serialized rows and native verification;
derivative checks masked the crust rows. **The 10 CSS survivors of
`CENSUS_b11a530.json` (and of every earlier census) were built on this
fabricated crust pressure and are withdrawn.** Earlier GPPVA_DD2 "physical
no-window" rejections were made on the same defective input.

Fail-closed correction (`5a4d0af`): `match_crust_core` rejects any source row
with provider-Euler residual > 5e-4 *before* the window search, reusing the
tolerance the sharp-Maxwell kernel already declares. Real providers pass with
large margin (APR/BBB ≤ 4e-16, GM1 ≤ 3e-11, RG ≤ 1.2e-4, GPPVA_FSU2H ≤ 1.6e-6).
The admission layer checks the new `source_euler_relative_residual` too.
Evidence: `data/receipts/core_crust_audit_20261115/CLD_COLUMN_DEFECT_EVIDENCE.json`.

Census at `5a4d0af` (`CENSUS_5a4d0af.json`, manifest
`91dc9087…77ba37`): 80 configured, 60 evaluated, 20 source-incomplete,
**0 survivors**. Failures: 30 Maxwell implementation_convention, 18 CSS
input_contract (source Euler: SPG×{APR,BBB,GM1}, GPPVA_DD2×{APR,BBB,GM1}),
12 CSS input_contract (RG beta flag), 20 data_provenance.

Not fixed: the writer/solver themselves and regeneration of the crusts. A
correct regeneration needs the solver's cell pressure (thermodynamically
consistent with E_cell), gas–cluster μ_n equality, the CompOSE column layout,
one normalization reference per table, and μ_q = μ_p − μ_n.

## B. Maxwell energy reference (adjudicated from sources)

All 30 Maxwell attempts pair a reconstructed crust (header 939/939) with a
provider core (939.56533/938.272; RG_SLY4 939.5654). Source evidence:

| Crust model | Hamiltonian nucleon masses (source) | Repo solver mass | Header written |
|---|---|---|---|
| GPPVA_DD2 | physical m_n, m_p (Typel et al. 2010) | physical | 939/939 |
| GPPVA_DDME2, GPPVA_TW | common 939 (Lalazissis 2005; Typel–Wolter 1999) | 939 | 939/939 |
| SPG_M1 | common 938.9 (Gögelein et al. PRC 77 025802, 2008) | physical (parameter set records none) | 939/939 |
| SPG_M2–M5 | not stated in repo; M2/M4 are placeholders equal to M1 | physical | 939/939 |

Verdict: the 939/939 header is a writer constant, not a provider or
Hamiltonian quantity, and the crust tables are defective anyway. These 30
attempts cannot be adjudicated as normalization-resolvable and are **not**
evidence about coexistence. Distinguish two future cases after regeneration:
(1) equal Hamiltonian masses (e.g. DD2 vs physical-mass cores): header
normalization only, physical ε, μ, P unchanged, gate passes legitimately;
(2) different constituent masses (DD-ME2/TW 939, SPG 938.9 vs physical
cores): a Hamiltonian-level difference that normalization cannot remove;
it must stay blocked or be studied as an explicit, labelled model choice. The
removed additive energy shift must not be reintroduced on either path.

## C. Source items (status)

- GPPVA_DDME2/TW: still missing eos.nb/eos.t/eos.yq; not manufactured.
- RG_SLY4/RG_KDE0v1: read back `beta_equilibrated=F`; contract not overridden.
- SPG_M1/M2/M4 duplicates: root cause found — M2 and M4 parameter sets reuse
  Gögelein 2008 Table I couplings because Scurto 2024 publishes them only as
  curves. They are not distinct models.
- Raw SPG Iphase-7 convexity and APR high-density acausality (c_s² ≈ 1.47):
  unchanged.

## D. Not started

Composition-aware Maxwell v2 and the generated CSS → phase-segmented
TOV/Cowling adapter. Phase-only output is not abundance preservation or
elastic-crust validation. With the crusts withdrawn, no generated matched EOS
currently exists to feed such an adapter.
