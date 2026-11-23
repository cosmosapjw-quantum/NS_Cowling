# BBP raw fixed-prefix rejection certificate — 2026-11-23

## Scope and result

This certificate resolves exactly the 145 native post-bridge trials in 135 cells of the frozen v2 campaign that used `secondary_literature:BBP_Diener2008_raw:crust`. It does not change the v2 journal, source catalog, source numbers, derivative operator, admission thresholds, or survivor lists. It does not decide Maxwell trials or any unlisted future input.

The complete set of later endpoint candidates in each affected prescribed phase view cannot pass the **existing native C11 contract**, because every such candidate retains the same first three crust rows and the first-point derivative is necessarily checked and negative. This is a constructive finite-grid exclusion under this representation and verifier, not a no-go result for the BBP physical model.

## Fixed rows and operator

The selected source has 26 rows. Native non-Maxwell endpoint search requires at least three retained crust rows: Fortran j ranges from 26 down to 3. For every j, the output starts with source rows 1, 2, and 3. Matching code reencodes their physical free energy and pressure using the target core mass header; it does not change these physical source quantities. The output chemical potential is the Euler value, with the Il=0 charge contribution treated according to the actual native encoding.

Let h1=n2−n1, h2=n3−n2, s12=(mu2−mu1)/h1, and s23=(mu3−mu2)/h2. The frozen verifier uses the quadratic one-sided endpoint derivative

```
D1 = ((2*h1+h2)*s12 - h1*s23)/(h1+h2).
```

The first three chemical-potential values themselves increase:
938.8768044781031, 938.9861021938987, 940.1564152841778 MeV (last displayed digits depend on lossless representation arithmetic).
The secants are approximately 2376.0373 and 14272.1109 MeV fm^3.
The one-sided quadratic derivative is nevertheless approximately **−1899.11413469 MeV fm^3**.
Thus the C11 failure must not be described as a decrease between the first two tabulated chemical-potential values.

The inspector applies the actual matching reencoding to every one of the **25 distinct affected core views**. The resulting D1 lies in
[−1899.1141346935183, −1899.1141346897930] MeV fm^3, below the frozen C11 lower threshold −1e−10 for every view. These four floating-point values reproduce the four values in the 145 original failure receipts. All affected views currently have Il=1 and Yq=0; the inspector nevertheless reproduces the native Il-dependent formula. This does not claim support for untested Il/Yq configurations.

## The first point cannot be masked

Only an interface on source interval 1 or 2 could mask the first-point three-row stencil. Neither interval qualifies:

- Delta log n = 0.1589545321399033 and 0.2326825743609691, both below the fixed `max_log_step=0.25`.
- All three source Iphase values are 0; no source phase boundary is present.
- Relative pressure differences are 0.05508982035928164 and 0.08241758241758232, far above the native equal-pressure jump trigger 1e−10; these intervals cannot be classified as a genuine Maxwell interface.
- A new bridge begins after at least source row 3. Its boundary can mask its own adjacent rows, but cannot alter the first-point stencil.
- For all 25 target-view reencodings, first-point dP/dn, dF/dn, d epsilon/dn, and d mu/dn are finite, and d epsilon/dn is positive (approximately 937.1792227907893 MeV). Therefore the invalid-derivative masking branch does not remove the first point.

The explicit source-prefix calls additionally report C11 failure for all 24 eligible prefix lengths. This check is corroborative only: the proof uses the unchanged **first** endpoint stencil, not a shortened prefix's final stencil. A whole-prefix consistency result must not be reused as a general prefilter, because its final stencil would differ from the assembled EOS.

Consequently every retained endpoint alternative in each affected view necessarily fails C11, independently of the new bridge and core tail. The 135 cell outcomes remain NO_ADMITTED_MATCH; the 145 initially incomplete-search trials can be adjudicated as `FIXED_RETAINED_PREFIX_C11_REJECTION` by this separate certificate.

## Separate cold first-law check

This section is a direct derivation plus native numerical check, not a newly introduced production gate.

For a stable, cold, one-charge branch, dP=n dmu and n increases along the branch. If mu2>mu1, then

```
n1*(mu2-mu1) <= P2-P1 <= n2*(mu2-mu1).
```

The exact immutable Euler source values violate even this necessary bound:

| Source interval | Delta P, MeV fm^-3 | n1 Delta mu, MeV fm^-3 |
|---|---:|---:|
| 1→2 | 2.8710941742519621e−5 | 2.9182490117426708e−5 |
| 2→3 | 4.6811318058455668e−5 | 3.6630799725736991e−4 |

Causality gives the stronger lower bound
n1*(mu2^2−mu1^2)/(2*mu1), which is also violated. Thus merely replacing the endpoint derivative by a monotonic interpolation would not establish exact cold first-law consistency of these immutable source knots. The table is a secondary-literature transcription at approximately three significant digits. These calculations do not identify the author's original rounding error or prove the underlying BBP model inconsistent. A higher-precision original or an explicitly new reconstructed input identity is required to change that source claim.

## Implementation and provenance

`inspect_bbp_prefix.f90` performs all physical calculations. Python is used only to map frozen journal records to source paths and create hashes/receipts. The inspector links the frozen type, IO, and verifier objects without editing production sources.

Build:

```sh
gfortran -O2 -ffree-line-length-none -fcheck=all \
  -I/tmp/ns_eos_reaudit_20261123/build/mod \
  inspect_bbp_prefix.f90 \
  /tmp/ns_eos_reaudit_20261123/build/obj/compose_eos_types.o \
  /tmp/ns_eos_reaudit_20261123/build/obj/compose_eos_io.o \
  /tmp/ns_eos_reaudit_20261123/build/obj/compose_eos_verify.o \
  -o inspect_bbp_prefix
```

Run with the frozen BBP internal directory and the 25-view path manifest:
`inspect_bbp_prefix BBP_INTERNAL_DIR AFFECTED_CORE_VIEW_PATHS`.
The path manifest remains in private intake because it may contain private source locations; the public log contains only source IDs, mass headers, and generic scalar diagnostics. It contains no private numerical table.

Evidence: `NATIVE_PREFIX_CHECK.txt`, `CERTIFICATE.json`, and matching-agent `V2_POST_BRIDGE_145_TRIAGE.json`. Exact source, executable, object, catalog, journal, and certificate artifact hashes are recorded in the JSON receipt.
