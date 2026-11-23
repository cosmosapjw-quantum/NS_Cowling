# Primary-source crust generation checkpoint — 2026-11-23 KST

Baseline: 7ad3040d93215db67b23041987e17ea19b6c6715 (PR20).
Branch: research/eos-primary-crust-generation-20261123.
Production physics remains Fortran. Python handles frontend, source registration, parallel dispatch, receipts and independent validation.

## Completed work

- Four explicit models: BPS1971_sequence, BBP1971_table5, BBP1971_pressure_integral, NV1974_author_fit. No old source identity was silently replaced.
- Native catalog expanded 18 -> 22; public dispatch/readback support; source diagnostics sealed.
- All 32 existing core source identities x 4 new representations x 7 methods executed: 896/896, with 672 previously blocked target cells and 224 separate derived-sensitivity cells.
- Native: 547 positives, 349 no admitted match, no runtime/protocol errors.
- Source adjudication: 410 representation/interface candidates, 224 raw-BBP first-law source rejections, 262 no admitted match. All 137 native raw-BBP positives overturned by a retained-prefix certificate.
- BPS 132, BBP pressure-integral 140, NV author fit 138. Inner-source 278 require source-accuracy and neutron-gas interpretation qualifications.
- 393 public candidate archives packaged; 17 private candidates retained locally. Every qualified archive contains SOURCE_ADMISSION.json.
- Full regression: 1985 passed, 35 skipped, 4 subtests passed. Logs in evidence/.
- Independent source/code/archive/publication-boundary review completed in review/. No remaining actionable finding within declared scope.
- Original source/catalog/journal and matching tolerances preserved; matcher executable identical to baseline.

## Canonical records

- docs/EOS_PRIMARY_CRUST_GENERATION_20261123.md: complete technical report, equations, qualifications, sources and historical execution commands.
- public/qualified_campaign/: publicly distributable adjudicated ledger and CSV catalog.
- public/native_campaign/: unadjudicated results, for comparison only.
- delivery/: two lossless shards and extraction index for 393 public candidates.
- research_bps/, research_bbp/, research_nv/: primary facts, derivations and Fortran validation evidence.
- Full mixed SOURCE_CATALOG.json, campaign_v1/, adjudicated/, adjudicated_v2/ are deliberately Git-ignored. They are preserved in the private local backup, not suitable for public publication.

## Remaining exact-model inputs

- BBP original full CLDM: DOI 10.1016/0375-9474(71)90281-8, or original unrounded thermodynamic table with conventions.
- NV1973 original HF paper and comparison against 1974 author fit: DOI 10.1016/0375-9474(73)90349-7.
- Full original BPS mass-candidate catalog: Myers & Swiatecki UCRL-11980, 27 May 1965 (verified report DOI unavailable).
- No full-star/TOV/Cowling validation is claimed for these new source-qualified candidates.

## Shipping and recovery

Publication is authorized to this research branch and a Draft PR based on research/eos-core-crust-all-methods-20261123; no merge or release.
Backup directory: /home/cosmosapjw/Dropbox/NS_Cowling_EOS_Primary_Generation_Backups/20261123/.
Final shipping receipt records actual commit, PR, backup hashes and readback results separately. An incremental Git bundle requires a repository already containing the baseline commit.
Original dirty /home/cosmosapjw/Dropbox/NS_Cowling is not used for edits.
