# Bounded follow-up review: primary-source adjudication

Date: 2026-11-23. This extends INDEPENDENT_PRIMARY_REVIEW.md only for the
source-adjudication and delivery representation. The original generator/kernel
review verdict is unchanged; no broad physics audit or test-suite rerun was made.

## Reviewed v1 script identity

scripts/adjudicate_eos_primary_crust_20261123.py
SHA256 ec044d40fd240af7dedc2a59e1e0cfebcfb1e547639daf85fb0105b0b5e529de

Actual output reviewed:
data/eos_primary_crust_generation_20261123/adjudicated/

## Confirmed logic and actual census

- The primary source file identity is pinned to the actual recovered BBP table.
- The Fortran certificate is bound to the generation receipt's recorded SHA.
- The certificate establishes an incompatible first interval, with conservative
  simultaneous printing boxes for n, P, and epsilon.
- Every one of the137 native-positive raw BBP candidates was checked against
  the two retained source knots before its status was changed.
- Original native records remain in campaign_v1; each adjudicated record keeps
  its native_status alongside the separate source-level status.
- Independent recomputation of the expected crust_source_id/core_source_id/method
  tuple set agrees exactly with all896 unique actual records.
- Actual statuses:224 REJECTED_SOURCE_FIRST_LAW,410 SURVIVED,262 NO_ADMITTED_MATCH.
- Remaining positives:132 BPS1971_sequence,140 BBP1971_pressure_integral,
  138 NV1974_author_fit.
- Public/private positives:393/17.
- All278 inner-representation positives occur exactly once in the explicit
  INNER_SOURCE_QUALIFICATIONS candidate ID list.
- All410 positive archive SHA256 and byte counts were independently checked,
  including the17 private archives, which remain outside the public output directory.
- Public artifact filenames equal exactly the393 public positive cell IDs;
  no raw BBP archive or private positive archive occurs in that public directory.
- Current native campaign files match the preserved pre-adjudication SHA manifest.
- Current catalog SHA matches the original native campaign settings.
- Frozen source-prefix CSV/NML hashes match their catalog record.
- No census, source-prefix application, archive-byte or privacy-path mismatch
  was found in the actual v1 outputs.

Exact identities:
- Native journal SHA256:
  4e87c5ef7df2dd38cc25572e6334cfd148b070fe660ea235a75a35cc5ac2365b
- Adjudicated v1 journal SHA256:
  5f3ea5e89953303f9705d78d95aa7f93bc6c2891113f3478d1df10a26bad2233
- Catalog SHA256:
  e4e8d25b3f94b4342950b2ae192b44d76447c1b56da5de1d1656ffe865d10e4d

Machine evidence: adjudication_identity_review.json in this review directory.

## Actionable delivery finding in v1

Severity: material provenance/delivery issue; numerical EOS and kernel bytes
themselves are correct and unchanged.

Location: scripts/adjudicate_eos_primary_crust_20261123.py:149-155.

The public outputs hard-link the original native candidate archives unchanged.
Source qualifications and admission scope were added only to outer catalog
records. A user moving one candidate archive separately would retain its original
PROVENANCE.json with numeric status SURVIVED, without the newly imposed source
accuracy qualification.

Reproduction: inspect adjudicated/artifacts/cell_000673.tar.gz.
Its15 original members include PROVENANCE.json but no SOURCE_ADMISSION.json.
The outer record identifies NV author-fit limitations and near-drip threshold
failure, while the embedded original provenance lacks those fields.

Required correction: preserve the original native gzip and all original member
bytes; deliver a separate qualified archive or wrapper carrying a machine-readable
source-admission sidecar bound to the original gzip SHA and exact cell ID.
The sidecar must state representation/interface scope and original-microphysics/
stellar non-claims. This issue was communicated to root and the bounded v2 patch
is in progress; v2 had not yet been inspected at this snapshot.

## Public-package privacy qualification

The393 public archive paths are correctly separated from17 private archives.
However, full v1 RESULTS.jsonl and SURVIVORS.json contain the private cells'
native_report endpoint numbers as well as their paths. Those all-cell receipts
must not be copied unfiltered into a public package if its privacy contract
excludes derived numeric content from the user's private core sources.
A public-only record subset, together with a metadata-only private count, preserves
the census without disclosing those values. This is a delivery observation;
the review does not reinterpret existing user publication authorization.

## V2 disposition: portable qualification finding resolved

The bounded patch and actual `adjudicated_v2` archives were subsequently
inspected. No actionable defect was found in that patch or its actual outputs.
The generator/kernel physics verdict is unchanged.

Reviewed script SHA256:
11fe35f70de3d3aa9b51b2b9ac3493c2d316d3c873a6ce84ca23e32881535e10

Adjudicated v2 journal SHA256:
e18cd3da16b3444efc3956bd81c5cebe219e0afa4b5617731b2aa43ea6ff53cf

Native journal SHA256, still unchanged:
4e87c5ef7df2dd38cc25572e6334cfd148b070fe660ea235a75a35cc5ac2365b

The new `qualify_archive` function (lines81-118) preserves the original native
gzip and creates a separate deterministic archive. It retains all original
member contents and adds `SOURCE_ADMISSION.json`. That sidecar binds the cell,
both source IDs, privacy flag, original native status, explicitly scoped final
status, source-accuracy limitation, original gzip SHA256, and original member hashes.
It explicitly declares original microphysics and a complete star unvalidated.

All410 qualified archives were opened alongside their original native archives.
For every candidate, the review compared every original member's bytes, both
gzip SHA records, all sidecar identity/scope fields, and the sidecar hash recorded
in the outer catalog. No mismatch was found.

Actual member counts:
- 408 candidates:15 original members ->16 qualified members.
- 2 exact-Maxwell candidates:14 original members ->15 qualified members.
  Their original archives have no PHASE_TRANSPORT.json; this is expected and
  does not indicate a missing qualification or changed data.

All393 public qualified archives occur exactly in the public artifact directory.
All17 private qualified archives remain under the separate private directory
outside the repository. The896-cell census and224/410/262 status split are
unchanged. All278 inner-representation positives retain their qualification IDs.

Machine evidence: `qualified_archive_review.json`.

The earlier v1 delivery finding is resolved for the actual v2 candidate archives.
The full mixed all-cell receipts still require the public projection described
above before public delivery. Root confirmed it is preparing a projection with
868 full public records plus28 private metadata-only outcomes, and35 public
sources plus one private metadata-only source. That separate projection/package
was not part of this bounded patch check and is not claimed verified here.
