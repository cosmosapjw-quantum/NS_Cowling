# Core–crust matching census — 2026-11-15 (code HEAD 843a399)

Receipt for the full configured census after the independent admission layer
(`scripts/match_report_admission.py`) was made a survivor requirement.

## Provenance

| item | value |
|---|---|
| code HEAD (clean tracked tree) | `843a399487b49fa90e177278609d9b663a5ed2fe` |
| native executable SHA-256 | `650057fe2a332633edc166639ef0151e003b5abb885dbc1237ee0b2b55822400` |
| run root (local, untracked, 28 MB) | `download/matched_runs/20261115_head843a399_census/` |
| `matching_manifest.json` SHA-256 | `1a7ed300d299b97893bcb2f73d2db96b573597ce47587156da324f320fd83a2e` |
| compact receipt | `data/core_crust_census_20261115_head843a399.summary.json` (per-pair outcome, category, `pair_report.json` SHA-256, survivor output SHA-256) |
| command | `NS_COWLING_MATCH_OUT_ROOT=download/matched_runs/20261115_head843a399_census python scripts/match_core_crust.py` (explicit empty root; no prior census touched) |

Two earlier exploratory runs from uncommitted code (`20261115T130356Z`,
`20261115T130510Z`) gave identical counts; they are kept locally and are not
evidence for this receipt.

## Counts

* Configured: 8 crusts × 5 cores × 2 methods = **80**. All 13 model
  directories exist, so missing directories do not reduce the count.
* **20 NOT_EVALUATED** (data provenance): `GPPVA_DDME2` and `GPPVA_TW` contain
  only `eos.compo` and `eos.thermo`; `eos.nb`, `eos.t`, `eos.yq` are absent.
  These are not failed physics cases.
* **60 evaluated** → **10 survivors, all CSS, all ADMITTED** by the
  independent layer; 0 Maxwell survivors.
* Native-PASS and independently-ADMITTED sets coincide (no admission
  rejections among native survivors).
* **Only 6 distinct survivors**: `SPG_M1`, `SPG_M2` and `SPG_M4` reconstructed
  crusts are byte-identical (`eos.nb`, `eos.t`, `eos.yq`, `eos.thermo`), so
  their survivors duplicate one another. All five SPG reconstructions share
  one `eos.nb` and have 21 rows.

| failure category | CSS | Maxwell | evidence |
|---|---|---|---|
| data_provenance | 10 | 10 | missing `eos.nb/eos.t/eos.yq` (GPPVA_DDME2, GPPVA_TW) |
| implementation_convention | 0 | 30 | `MASS_CONVENTION_FAIL`: the sharp-Maxwell kernel refuses crust/core tables with different nucleon rest masses; the CSS path applies a representation-only gauge shift instead. No Maxwell physics was tested. |
| input_contract | 12 | 0 | `RG_SLY4`, `RG_KDE0v1` cores read back with `beta_equilibrated=F`; native cold matching refuses them. Not overridden. |
| native_consistency_gate | 5 | 0 | all `*+APR_APR`: C10 causality, max c_s² ≈ 1.47, 24 acausal points; two independent derivative estimators agree (1.4709 / 1.4707). APR core source is superluminal at high density. |
| physical | 3 | 0 | `GPPVA_DD2` + {APR, BBB_BHF2, GM_GM1}: no causal single-speed CSS window with a positive crust density jump. |

## Outcome grid

| crust \ core | APR_APR | BBB_BHF2 | GM_GM1 | RG_SLY4 | RG_KDE0v1 |
|---|---|---|---|---|---|
| SPG_M1 | css: native_consistency_gate<br>mx: implementation_convention | css: ADMITTED<br>mx: implementation_convention | css: ADMITTED<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention |
| SPG_M2 | css: native_consistency_gate<br>mx: implementation_convention | css: ADMITTED<br>mx: implementation_convention | css: ADMITTED<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention |
| SPG_M3 | css: native_consistency_gate<br>mx: implementation_convention | css: ADMITTED<br>mx: implementation_convention | css: ADMITTED<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention |
| SPG_M4 | css: native_consistency_gate<br>mx: implementation_convention | css: ADMITTED<br>mx: implementation_convention | css: ADMITTED<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention |
| SPG_M5 | css: native_consistency_gate<br>mx: implementation_convention | css: ADMITTED<br>mx: implementation_convention | css: ADMITTED<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention |
| GPPVA_DD2 | css: physical<br>mx: implementation_convention | css: physical<br>mx: implementation_convention | css: physical<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention | css: input_contract<br>mx: implementation_convention |
| GPPVA_DDME2 | css: data_provenance<br>mx: data_provenance | css: data_provenance<br>mx: data_provenance | css: data_provenance<br>mx: data_provenance | css: data_provenance<br>mx: data_provenance | css: data_provenance<br>mx: data_provenance |
| GPPVA_TW | css: data_provenance<br>mx: data_provenance | css: data_provenance<br>mx: data_provenance | css: data_provenance<br>mx: data_provenance | css: data_provenance<br>mx: data_provenance | css: data_provenance<br>mx: data_provenance |

## Survivors (CSS, single sound speed, crust-side density jump)

Admission residuals are recomputed from the serialized readback rows. Gates
are native declared values: latent identity 1e-10, CSS rows 5e-9, CSS→core
endpoint 2e-10·max(1,|x|). Observed residuals are ≥4 orders below gate.

| pair | c_s² | Δn_B (fm⁻³) | Δε (MeV fm⁻³) | Seidov ratio | latent-id rel. resid. | CSS-row max rel. resid. |
|---|---|---|---|---|---|---|
| SPG_M1+BBB_BHF2 | 0.0874 | 5.097e-03 | 4.854e+00 | 0.070 | 2.6e-15 | 3.2e-14 |
| SPG_M1+GM_GM1 | 0.1928 | 4.564e-04 | 4.346e-01 | 0.006 | 6.5e-15 | 1.4e-14 |
| SPG_M2+BBB_BHF2 | 0.0874 | 5.097e-03 | 4.854e+00 | 0.070 | 2.6e-15 | 3.2e-14 |
| SPG_M2+GM_GM1 | 0.1928 | 4.564e-04 | 4.346e-01 | 0.006 | 6.5e-15 | 1.4e-14 |
| SPG_M3+BBB_BHF2 | 0.1296 | 5.419e-02 | 5.199e+01 | 0.747 | 6.8e-16 | 3.7e-14 |
| SPG_M3+GM_GM1 | 0.5128 | 2.736e-02 | 2.625e+01 | 0.377 | 4.1e-16 | 4.3e-14 |
| SPG_M4+BBB_BHF2 | 0.0874 | 5.097e-03 | 4.854e+00 | 0.070 | 2.6e-15 | 3.2e-14 |
| SPG_M4+GM_GM1 | 0.1928 | 4.564e-04 | 4.346e-01 | 0.006 | 6.5e-15 | 1.4e-14 |
| SPG_M5+BBB_BHF2 | 0.1216 | 4.753e-02 | 4.556e+01 | 0.655 | 4.7e-16 | 3.8e-14 |
| SPG_M5+GM_GM1 | 0.6379 | 2.519e-02 | 2.414e+01 | 0.347 | 1.0e-15 | 2.8e-14 |

## Not claimed

* Composition is dropped in both methods; no composition-preservation claim.
* The CSS survivors are smooth only inside the bridge; the crust→CSS interface
  is a first-order jump by construction.
* SPG_M1 survivors use the reconstructed crust in `download/reconstruction/`;
  the raw-source Iphase-7 convexity issue (99/147 intervals) is unresolved
  and is not addressed by this census.
* No Maxwell survivor exists, so this census provides no evidence for or
  against sharp-Maxwell crust–core matching; the blocker is the rest-mass
  convention contract, not a physics result.
* APR+CSS thermodynamic-potential SSOT, composition-aware Maxwell v2,
  crust-CSS→phase-segmented TOV/Cowling adapter and elastic-crust physics
  remain open.
