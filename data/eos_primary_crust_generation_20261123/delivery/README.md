# Source-qualified primary-crust candidates

This delivery contains 393 public representation/interface candidates in two lossless solid tar.xz shards. Another 17 candidate archives depend on a private core source and are retained in the private backup only. The index preserves metadata for all 410 positives.

These are source-qualified computational candidates. BPS is restricted to the original Table 2 nuclear sequence. BBP pressure-integral is a separately derived EOS, not the original full CLDM or a rounding-preserving repair. NV is Negele's 1974 author fit, not a rerun of the 1973 HF calculation. Every candidate carries SOURCE_ADMISSION.json and retained-source qualifications. No raw BBP Table 5 candidate is included.

From the repository root, extract one candidate using the checked extractor:

    python3 scripts/extract_eos_reaudit_candidate_20261123.py --index data/eos_primary_crust_generation_20261123/delivery/DELIVERY_INDEX.json --cell cell_000056 --out /tmp/ns_cowling_primary_cell_000056

The output directory must be new. The extractor validates shard and member SHA256 identities. Read SOURCE_ADMISSION.json before using the EOS. Cell IDs identify method/source trials, not unique physical equations of state.

See docs/EOS_PRIMARY_CRUST_GENERATION_20261123.md and public/qualified_campaign/SURVIVOR_CATALOG.csv for scope and counts. The mixed private/full journals are intentionally absent from Git.
