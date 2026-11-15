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

## Regenerated crusts and census on corrected inputs

**Writer/solver fix (`489d87f`, `6a2417e`).** `emit_compose_set` now builds
eps = E_cell + m_n n_n + m_p n_p with the model's own masses (E_cell is
measured from them), derives mu_B = m_n + R'(n_B) from one not-a-knot cubic
spline of R = eps − m_n n_B and P = n_B mu_B − eps, writes the CompOSE layout
with the model masses in the header, and refuses (writing nothing) unconverged
cells or tables whose eps, mu_B (convexity) or P are not strictly increasing,
P ≤ 0, or row secants exceed 1. Rows whose CLD minimum lies on the optimizer
box (f ∈ {1e-4, 0.95}, R_d ∈ {2, 40} fm) are not stationary, so the envelope
derivative is not thermodynamic there: contiguous end runs are trimmed and
recorded, interior ones refused. The v2 census exposed this — all 6 v2
survivors attached at a box-bound last row — so v2 is kept only as evidence.

**Regeneration** (`scripts/regenerate_cld_crusts_v2.py`, σ₀ = 1, b = 1, which
reproduces the stored SPG_M1 row exactly; deterministic, byte-identical reruns):

| Model | v3 result |
|---|---|
| SPG_M1 | 13 rows, n_B 0.0048–0.1116 (dropped 0–5 on R_d = 2 fm, 19 on f = 0.95) |
| SPG_M3 | REJECTED: box-bound minima inside the table (rows 1–2) |
| SPG_M5 | 15 rows, n_B 0.001–0.039 (dropped 15–19) |
| SPG_M2, SPG_M4 | NOT_GENERATED: placeholders of SPG_M1 |
| GPPVA_DD2 | 19 rows (last row dropped) |
| GPPVA_DDME2, GPPVA_TW | 14 rows each (last row dropped); eos.nb/eos.t/eos.yq now exist |

Diagnostic, not a gate: the solver's cluster μ_n differs from the derived μ_B
by up to 3.1 MeV (DD2) because the CLD prototype does not impose gas–cluster
neutron equilibrium. SPG reconstructions use the solver's physical nucleon
masses although Gögelein 2008 specifies M = 938.9 MeV; the SPG crusts are
therefore not the published Hamiltonian (source-provenance defect, open).

**Census at `e434bee`** (`CENSUS_e434bee.json`, manifest `dac5e45a…04c6`,
crust root `download/reconstruction_v3`): 80 configured, 50 evaluated,
30 not evaluated (20 SPG_M2/M4 placeholders = data_provenance, 10 SPG_M3 =
source_reconstruction_rejected). **10 CSS survivors, 10 distinct**, all
independently ADMITTED: {SPG_M1, SPG_M5, GPPVA_DD2, GPPVA_DDME2, GPPVA_TW} ×
{BBB_BHF2, GM_GM1}. Rejections: 25 Maxwell mass convention, 10 CSS RG beta
contract, 5 CSS APR C10 (max c_s² 1.4709, APR core). Seidov ratios exceed 1
for SPG_M5 (1.52, 3.79), DD2+GM1 (1.57), DDME2+GM1 (1.80), TW+GM1 (1.73):
small-core onset instability per Seidov is reported, not gated.

Scope of these survivors: thermodynamically admissible composite cold
barotropes built from a CLD prototype crust (no gas–cluster μ_n equality,
optimizer box, SPG mass deviation) plus one CSS segment. They are not
certified microscopic crust models and carry no composition.

## B re-adjudicated on v3 crusts

Headers now equal the Hamiltonian masses used to compute each crust, so the
sharp-Maxwell gate (1e-12 relative) compares constituent masses:

| Crust (Hamiltonian m_n / m_p, MeV) | vs APR/BBB/GM1 (939.56533 / 938.272) | vs RG_SLY4 (939.5654 / 938.272) |
|---|---|---|
| GPPVA_DD2, SPG_M1, SPG_M5 (939.56536 / 938.27203) | Δ = 3e-5 MeV each; mass-constant version difference | Δm_n = 4e-5 |
| GPPVA_DDME2, GPPVA_TW (939 / 939) | Δm_n = −0.565, Δm_p = +0.728 | same order |

All 25 Maxwell attempts remain blocked by the gate, now for a source-backed
reason: different constituent masses in the two Hamiltonians. For DDME2/TW
this is a sub-MeV Hamiltonian difference that normalization cannot remove.
For DD2/SPG it is a 3e-5 MeV difference of tabulated constants; admitting it
would require an explicit, labelled decision on mass-constant equivalence,
not a tolerance change. No coexistence result is claimed or excluded.
