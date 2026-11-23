# Immutable campaign_v2 delivery preservation

Status: all archive identities and representative extractions verified. This report
certifies byte preservation of campaign_v2. It does not certify exhaustive
matching search; the parent is separately triaging 145 post-bridge verification
failures. No new physical admission is performed by packaging or extraction.

## Contents and compression

| Quantity | Value |
|---|---:|
| Source census | 105 = 32 core + 73 crust |
| Campaign cells | 16,352 |
| Survived cells | 9,377 |
| Public survivors fully archived | 8,972 |
| Private survivors metadata only | 405 |
| All private cells metadata only | 728 |
| Solid tar.xz shards | 36 |
| Original public gzip bytes | 753,066,466 |
| Preserved raw tar bytes | 5,327,052,800 |
| Published shard bytes | 69,840,608 |
| Delivery index bytes | 6,957,074 |
| Complete delivery bytes | 76,797,682 |
| Largest shard bytes | 3,239,552 |
| Gzip / shard size ratio | 10.78264476 |
| Shard size reduction | 90.725837% |
| Packaging seconds | 188.292196 |
| Independent verification seconds | 4.969727 |

DELIVERY_INDEX.json SHA256:
d9695d8b08bc4ea82c0ba6e815d49bf8d4b53bc972ddbeca5e27ab79a86a06e2

## Verified contracts

All 8,972 original public gzip hashes, all 36 compressed shard hashes, and every
raw tar member hash were independently re-read. Each member is the complete
original tar from gzip decompression. Scientific tables, native reports,
readback certificates and provenance remain in that original tar.

Maxwell cell_001617 (14 files), CSS cell_000230 (15 files), and first-law bridge
cell_000225 (15 files) were extracted by the shipped frontend. Every extracted
file hash and byte length equals its corresponding original gzip member.
Private survivor cell_000442 was rejected without creating an output directory.

All five frozen campaign input files remained unchanged. The original v1
delivery index and all 33 v1 shards were independently hash-verified unchanged.
All original candidate gzip files remain present for backup. No private numeric
table or private trace appears in published shards.

The production packager/extractor/delivery tests remain byte-identical to the
previous verified versions (18 pytest cases passed previously). No production
code was changed for v2. See CODE_SHA256_FINAL.txt.

## Evidence

- PACKAGING.log: completed packager output and wall time.
- VERIFICATION.json: full shard-level independent audit and extraction receipts.
- VERIFICATION.log: compact result.
- verify_delivery.py: independent validator used for this package.
- extracted/: three representative exact file trees in scratch, no private data.

This package is retained while final scientific triage is resolved separately.
