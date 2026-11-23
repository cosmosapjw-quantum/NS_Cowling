# Final scoped review: public projection and delivery package

Date: 2026-11-23.
Scope: publication privacy, census integrity, source-qualification portability,
and exact packaging bytes. No new physics audit or test-suite run was performed.
This does not extend the earlier representation/interface results into full
original-microphysics or complete-star admission.

## Disposition

No actionable discrepancy was found in the actual public projection and package.
The earlier private-numerical-metadata publication gap is closed in the inspected
public records. The portable source-qualification gap remains resolved inside
every inspected public candidate in the delivery shards.

## Reviewed identities

- Projection script:
  scripts/publish_eos_primary_crust_projection_20261123.py
- Script SHA256:
  b77940171d1df0c62618d07b90a053daed95141e99dc255b15a75d7cf6cc83b5
- Public source catalog SHA256:
  29be5f470b3cd637121d93b82698432fd7eb29f16e0da47c94185a0137fca817
- Public qualified journal SHA256:
  9bbe02601f73333341134a26e418bb4b80c21bf38ce2c1102e7c4e3d3ec20bef
- Delivery index SHA256:
  b25c5d95d73e46ae014c077ce547857d57dfdb04a4346d035021df8f9f1bab7f

## Public projection actually checked

Both the native and qualified public ledgers contain896 unique cell entries:
868 complete public records and28 strictly whitelisted private metadata records.

For each private record, the review checked allowed keys and scalar value types,
rejected nested dictionaries/lists/floating-point values, and confirmed the absence
of native_report, verification, artifact, failure_trace, source_view_bounds and
endpoint/table fields. All retained metadata values match their private originals.

The qualified survivor list includes393 complete public entries and17 private
metadata-only entries. The native survivor list has524 public entries and23 private
metadata-only entries. The source-adjudication change therefore remains visible
without publishing private numerical reports.

For each public record, the review compared the full record with its original,
allowing only the documented public artifact path and raw-receipt scope projection.
Statuses, IDs and source/method identity are unchanged.

The public catalog contains35 complete public sources, semantically identical as parsed
records to their originals, plus one private source with the explicit whitelist:
source_id, role, model, private, available, status, publication,
execution_source_available, and empty matching_views. It contains no private
table, density array or endpoint data.

The public catalog's original-catalog SHA and the delivery index's campaign file
SHA inventory match the actual files. Copy-through sidecars contain public-source
certificates, count summaries and path/hash provenance rather than private EOS tables.

## Every delivery payload actually checked

Both .tar.xz shards were opened and all393 nested candidate .tar members were read.
Each member equals gzip.decompress(the corresponding qualified candidate archive)
byte-for-byte. Every raw tar size/hash, original qualified gzip size/hash, index
locator and shard size/hash was checked.

Every candidate contains exactly one SOURCE_ADMISSION.json. Its cell ID, source
IDs, private=False, representation-specific final status, source accuracy
qualification, admission scope and qualification SHA match the qualified ledger.
It explicitly marks full original microphysics and complete-star validation false.

No private candidate archive and no rejected raw-BBP candidate occurs in either
shard. The public package consists of:
- 134 BBP1971_pressure_integral candidates;
- 127 BPS1971_sequence candidates;
- 132 NV1974_author_fit candidates.

Shard evidence:
- shard_00000.tar.xz:256 members,2,089,412 bytes;
  SHA256795a1027effd06db5182f8f04faa7863ec2fe07c8c4442f0c80187bce57d0bed.
- shard_00001.tar.xz:137 members,1,037,164 bytes;
  SHA2560b0b329a3a42c85cc22432b6954b83fe1ee685024fe0c912d4ee0957a0df0580.

Total published shard bytes:3,126,576.
The delivery index's28 private entries contain metadata only and no archive,
native-report or endpoint fields.

## Original evidence and repository boundary

The native journal remains:
4e87c5ef7df2dd38cc25572e6334cfd148b070fe660ea235a75a35cc5ac2365b.

The complete qualified v2 journal remains:
e18cd3da16b3444efc3956bd81c5cebe219e0afa4b5617731b2aa43ea6ff53cf.

The task .gitignore excludes the complete mixed SOURCE_CATALOG.json,
campaign_v1/, adjudicated/, and adjudicated_v2/. Both git check-ignore and a
read-only git ls-files check were performed; none of those mixed paths was
present in the Git index at the inspection time.

Backup creation, commit and remote publication are root's subsequent operations
and are not claimed completed by this review.

Machine evidence: public_delivery_review.json.
