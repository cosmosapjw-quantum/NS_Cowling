# Independent integration review — 2026-11-16

Scope: one bounded review of the new public generation diagnostics, native build
integration, receipt hashing, and core-domain/crust-attachment API risk points.
This review did not change production code, invoke another reviewer, repeat the
model-family research, or rerun the parent's full suites.

## Findings and disposition

### 1. Sealed optional diagnostic deletion — fixed and independently checked

The initial frontend diff hashed optional diagnostic files only when they existed.
The final comparison iterated only the files currently found. Consequently,
deleting a diagnostic already recorded in `RESULT.json.output_sha256` removed
that key from the comparison and verification still returned `rows=3`.

This was actually reproduced in a fresh Python process with
`_read_native_table` mocked to isolate receipt control flow. The mock supplied
only the fixed three-node grid, five native-file hashes and empty consistency
diagnostics; it did not supply or override optional-file handling.

The parent added an explicit check for any sealed optional key missing from the
current hash set, and a parameterized deletion regression covering
`sound_speed.csv` and `attachment_diagnostics.json`. The same independent
control was rerun after reading the updated source:

| Case | Observed result |
| --- | --- |
| Delete sealed sound_speed.csv | ValueError: sealed native diagnostic missing |
| Delete sealed attachment_diagnostics.json | ValueError: sealed native diagnostic missing |
| Historical receipt with no optional diagnostic key/file | Accepted |

Disposition: **CLOSED**. This is evidence-integrity validation; physical values
such as c_s²/c²>1 and an unsupported attachment interval remain diagnostics rather
than new generator acceptance gates.

### 2. New direct sound-speed APIs accepted out-of-domain points — fixed and independently checked

The CLI writer guards finite-difference stencils using declared model bounds.
However, the new public native routines initially called `solve_core_state`
without their own declared-support check. That existing function obtains model
bounds but does not itself enforce them, and several model kernels evaluate a
larger mathematical interval.

An actual standalone native probe, built against the reviewed source, returned:

| Model | Declared upper n (fm⁻³) | Queried n | FD status, c_s²/c² | Analytic status, c_s²/c² |
| --- | ---: | ---: | --- | --- |
| RG_SLY4 | 1.4 | 1.5 | ok, 1.10113886516250 | ok, 1.10113886517076 |
| GM_GM1 | 1.0 | 1.1 | ok, 0.77368947344736 | Unsupported analytic model, correctly non-ok |
| APR_APR | 1.34 | 1.5 | ok, 1.56184481484275 | ok, 1.56184481484016 |

Probe source/executable: `/tmp/ns_integration_api_probe_20261116/probe.f90`
and `probe`. Compilation used only types, core-native, core-domain and this
small program, with a separate module directory and one process.

The applied minimal correction queries `core_model_bounds` in the new routines;
it rejects an analytic evaluation outside [lo,hi] and any FD stencil with
n−2h<lo or n+2h>hi. Preserve existing physical values for valid queries.
The new `core_domain_bounds.f90` includes outside central points, inside central
points with stencils crossing both endpoints, inclusive analytic endpoints,
nonfinite inputs, and valid interior/crossing controls (20 cases).

Impact: the current CLI request and FD-wrapper checks protect the parent public
smoke outputs. This finding does **not** invalidate those numerical tables or
require repeating the full model campaign. It concerns the new direct API's
support/status contract.

The fixed source was read at SHA256
`683ef46c86fe72a9b535dc22d0511e15c93905b3879311ef8abc3e9658930c50`.
The exact original standalone probe was rebuilt against it and rerun. All six
FD/analytic out-of-domain calls now return
`density outside declared native core model bounds`. The bounds-test source was
inspected; the core lane separately reports its 20-case execution with exit 0.
This reviewer did not rerun that full test after independently confirming the
original failure cases.

Disposition: **CLOSED**. No unresolved integration blocker remains in this review's
bounded scope; final build-signature and end-to-end execution ownership stays
with the parent.

## Inspected integration contracts

- `fortran_eos/build.sh` compiles core-domain after core-native and
  crust-attachment after lower-crust-native, then includes both in the native
  executable. Source-signature scanning already includes new native source
  files. No dependency-order or missing-object defect was found.
- `write_core_domain_diagnostics` writes unavailable results as empty CSV
  fields and preserves explicit status strings. It does not emit a fake zero
  sound speed on failure. At BBB energy knots, phase crossings and muon onset
  crossings, the native FD routine returns a non-ok status. APR interior mixed
  states use the constrained Gibbs response.
- The single-phase analytic formula was checked against the independent
  density/composition variables. With total H(n,y), eliminating equilibrium
  composition gives H_nn−H_ny²/H_yy. The implemented lepton additions
  y²/chi, mu+n*y/chi, n²/chi have the correct derivatives and units.
  The output n*(H_nn−H_ny²/H_yy)/mu_B is dimensionless c_s²/c².
- The Gibbs transformation from (n,y) to (n,q=n*y), the inverse phase
  susceptibilities, positive lepton charge susceptibility, and coexistence
  tangent (Delta q,−Delta n) have consistent signs. The result
  n*(Delta q)²/[mu_B*(v.K.v)] matches the stated global-neutrality contract.
  This was a formula/code inspection; the lane's numerical comparisons and
  parent suites supply the separate execution evidence.
- `eos_endpoint` reconstructs epsilon, P and mu_B from each segment's own
  neutron-mass header. It does not mistake a header normalization change for a
  new physical energy shift. The source conventions remain those of the
  underlying models; the attachment code does not fit an arbitrary offset.
- Attachment diagnostics require finite positive endpoints and bracketing of
  the nominal source boundary. Pressure/chemical-potential differences are
  retained rather than forced to zero. An accepted source prescription never
  marks equilibrium or the unsampled interval as established.
- `source_attachment_supports` certifies only the two recorded endpoints,
  matching its explicit cross-seam scope. Neither a gap nor an overlap becomes
  an automatically interpolable EOS.
- The BSk reference-integral wrapper fixes its anchor from the source energy,
  keeps the same rest-energy convention, and distinguishes a conditional
  pressure-envelope bound at equal baryon density from empirical fit errors
  at equal mass-energy density. The factor r/(1−r) is consistent with the
  documented relative-to-true-pressure inequality. Quadrature estimates stay
  separate and are not advertised as rigorous enclosures.

No additional blocker was found in these inspected paths. In particular, the
new attachment `c4dd6f0762fad107997f31c0ae3b1ea39c1cb10ed6f2b80813c13ba38b8955c2`
and unchanged pre-existing lower module
`56056c76bfcd66a6b36ae5bdd5e016494547d955f296030237fd9705540ba6b5`
are different files; there was no baseline mismatch.

## Publication hygiene performed on the BBB lane

Full arXiv TeX/archive and full external RNS code/data were moved into ignored
`build/blocker_bbb_sources`. The repository retains source URLs, hashes,
original numeric-table extraction and measured results. A GitHub metadata/root
listing check established no explicit RNS redistribution license
(license metadata null; no license/copying file in the pinned root), so the
complete RNS files are not included in publication.

`SOURCE_RELOCATION.json` records old-to-cache paths and byte identities.
The historical scientific `result.json` remains byte-identical. The BBB runner
now uses the ignored source cache, can download missing originals with required
SHA256 checks, and supports `--source-cache`. A replay after relocation exited
0 and reproduced four retained scientific output files byte-for-byte.
The refreshed BBB checksum manifest was verified.

## Review boundary

The parent owns the public five-core plus DD2/lower smoke campaign, full suites
and final build-signature check. Their results are not relabelled as independent
tests performed by this reviewer. Only the two targeted receipt controls,
direct API probe, source-cache replay and inspected formulas/contracts described
above are this review's execution evidence.
