# Lossless EOS delivery packaging — 2026-11-23

## Completed result

The completed 15,190-cell campaign contains 8,712 surviving candidates. This delivery publishes all **8,318 public survivors**, preserving each complete original candidate tar byte-for-byte. The **394 private survivors** and all private failure traces remain outside the public shards. The index retains only whitelisted identity/status metadata for all **700 private cells**.

Public output directory:

`data/eos_matching_reaudit_20261123/delivery_v1`

| Measurement | Bytes |
|---|---:|
| Original public per-cell gzip archives | 717,997,877 |
| Original decompressed candidate tar content | 5,083,412,480 |
| 33 published solid-XZ shards | 66,247,512 |
| DELIVERY_INDEX.json | 6,455,803 |
| Complete published delivery | 72,703,315 |
| Largest single shard | 3,321,184 |

The shards alone are **10.8381 times smaller** than the existing individual gzip archives, a **90.7733%** reduction. Including the index, publication is **9.8757 times smaller**, a **89.8742%** reduction. This is lossless compression of complete retained candidate data, not bridge-only reconstruction or numerical approximation.

Packaging used 256 cells per shard, XZ preset 6, and four concurrent shard writers. The recorded packaging interval was 199.5896 seconds; the full CLI process exited 0 in 201.67 seconds. An observed resource sample showed 369% CPU and approximately 578 MiB RSS. This was an observation, not a peak-memory bound.

## Identity and privacy contract

Each shard entry is named `cell_XXXXXX.tar` and equals `gzip.decompress(original_cell_archive_bytes)`. The index maps each public cell to its shard/member and records the original gzip SHA256 and byte count, raw tar SHA256 and byte count, and shard SHA256 and byte count.

Original gzip files are retained unchanged for full backup; the delivery does not claim to regenerate their exact compressed bytes. The original tar bytes and every archived scientific/provenance file are preserved exactly. No candidate source coverage, physics gate, admissibility classification or numerical content is recomputed by these frontends.

Publication requires a COMPLETED campaign with unchanged input/binary identities, complete cell accounting, matching status counts, matching SURVIVORS and RESULTS records, and identical SETTINGS. SUMMARY, SURVIVORS, RESULTS, SETTINGS and INPUT_IDENTITIES hashes are verified again before publishing the final fresh output directory. Candidate source gzip hashes and internal PROVENANCE identities/private flags are checked before writing any public member. Existing outputs are not reused or overwritten.

Private index entries contain only cell/source/method/status/coverage/view identifiers and a metadata-only publication label. They contain no private native report, artifact path, numerical EOS, or failure trace. Public shard membership is disjoint from all private cell IDs.

## Independent verification actually executed

`VERIFICATION.json` records an independent full-delivery verification, exit 0 in 6.74 seconds (verification interval 6.4149 seconds):

- Rechecked all 8,318 public original gzip SHA256 values against both the campaign and delivery index.
- Rechecked all 33 compressed shard SHA256 values.
- Streamed every shard and rechecked all 8,318 original raw tar hashes and byte counts, covering 5,083,412,480 uncompressed bytes.
- Checked exact public membership and private exclusion against the completed campaign.
- Verified all campaign JSON/journal hashes remained unchanged.
- Extracted Maxwell `cell_001561` (14 files), CSS `cell_000223` (15 files), and first-law bridge `cell_000218` (15 files); every restored file name, byte count and SHA256 matched the original gzip archive.
- Attempted to extract private survivor `cell_000428`; the frontend rejected it without creating an output directory.

Final delivery-index SHA256:

`e2fbcc6c42433039fde00684d45ba5b3b8268e099507b2882c58b19656c5ea31`

## Frontend usage

From the repository root, restore one full public candidate into a fresh directory:

```sh
python3 scripts/extract_eos_reaudit_candidate_20261123.py \
  --index data/eos_matching_reaudit_20261123/delivery_v1/DELIVERY_INDEX.json \
  --cell cell_001561 \
  --out /tmp/restored_eos_cell_001561
```

The extractor validates the complete shard hash and selected raw-tar hash before restoring files. It rejects absolute paths, parent traversal, Windows drive/backslash paths, symlinks, hardlinks, special files, sparse members, duplicate paths and file/parent collisions. It uses a fresh staging directory and atomically publishes only to a newly reserved output directory. It never runs archived content.

Reproduce packaging from a complete campaign, using a new output path:

```sh
python3 scripts/package_eos_reaudit_20261123.py package \
  --campaign data/eos_matching_reaudit_20261123/campaign_v1 \
  --out /tmp/reproduced_eos_delivery \
  --shard-cells 256 --workers 4 --xz-preset 6
```

A separate `benchmark` mode permits a bounded completed-public-cell snapshot only under /tmp and emits BENCHMARK_INDEX.json with status BENCHMARK_ONLY. The extractor refuses benchmark indexes as public completed deliveries. The final publication above used `package`, not benchmark mode.

## Tests and implementation provenance

New files:

- scripts/package_eos_reaudit_20261123.py
- scripts/extract_eos_reaudit_candidate_20261123.py
- python_frontend/tests/test_eos_reaudit_delivery_20261123.py

Actual tests: new-feature RED observed first, then **18 passed in 0.21 seconds** in `GREEN_FINAL_PY310.log`. Tests exercise complete opaque-byte roundtrip, private exclusion, incomplete/changed campaign rejection, original gzip and published hash corruption, traversal/link/special-file/duplicate-path rejection, and existing-output preservation. A malformed duplicate test fixture was corrected to contain two actual file records; its initial result remains in GREEN.log.

The project declares Python>=3.10. A dedicated RED simulation showed that the first SHA helper depended on Python 3.11's hashlib.file_digest; it was replaced with the equivalent streaming hashlib.sha256 interface and all 18 tests passed. The already running packaging process retained its loaded SHA helper; the code snapshot and SHA lists at packaging start and after the compatibility correction are separately preserved as EXTRACTOR_AT_PACKAGE_START.py, CODE_SHA256.txt and CODE_SHA256_FINAL.txt. Both paths compute SHA256; the independent final verification ran with the final extractor. No EOS calculation was added to Python.

The running campaign program, campaign JSON/journal and source tables were not edited. No commit was made by this packaging worker. Original public/private gzip archives remain available for the parent's complete backup.
