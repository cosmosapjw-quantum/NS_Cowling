# NS_Cowling EOS390: Cloud Execution Developer Contract

**Version:** 1.3.0  
**Date:** 2026-11-11  
**Repository:** \`cosmosapjw-quantum/NS_Cowling\`  
**Execution branch:** \`cloud/eos390-fortran-unification-v1\`  
**Parent feature branch:** \`research/eos390-fortran-unification-v1\`  
**Pinned parent commit:** \`b12ba451b6755bf8f13d22c882075a23eae800a7\`  
**Machine authority:** \`data/eos390_native_fortran_unification_dag.v1.json\`  
**Model-paper registry:** \`data/eos390_source_bibliography.v1.json\`  
**Existing sharp-Maxwell 65-pair DAG:** separate, immutable, not replaced by this plan.

## 0. Non-negotiable task definition

The end goal is to regenerate exactly the original **390 crust × core × method** outputs with one scientifically authoritative Fortran backend and a thin Python frontend. The historical 390 outputs and the raw EOS source tables are immutable reference material. The objective is not to force every model through every method: an input EOS that fails causality, composition closure, or source-data validation must remain a failed or diagnostic-only result with its original evidence.

The execution target for this handoff is an **ephemeral cloud sandbox** (GitHub Actions or an authorized online coding-agent workspace). Do not use the user's desktop, WSL, local Dropbox mount, or a guessed local filesystem path. A tool that merely reads GitHub metadata is not an execution sandbox; use the cloud runner for compiler, tests, numerical runs, and artifacts.

The single source of truth is:

1. Current checked-out branch HEAD for code.
2. \`data/eos390_native_fortran_unification_dag.v1.json\` for machine-readable authority, dependencies, statuses, gates, and next actions.
3. This document for the human-readable contract and copyable agent prompt.
4. \`data/eos390_source_bibliography.v1.json\` for EOS-specific literature provenance.
6. \`data/eos390_reconstructed_source_inventory.v1.json\` for the surviving replacement source cohort.
5. The surviving v3.2 source ZIP, its delivery manifest, the checked-in 18-model/158-file inventory, and a freshly generated per-file SHA-256 manifest.
6. The v3.2 Dropbox release and old \`/ns_eos/matched_eos_all_combinations\` as immutable historical baselines, not as new-run evidence.

Do not create another EOS390 plan, another competing schema, or another “improved” output authority. Extend the existing DAG in place.

## 1. Frozen 390 task space and historical truth

The exact matrix is:

| Axis | Count | IDs |
|---|---:|---|
| Crust models | 13 | \`SPG_M1\`, \`SPG_M2\`, \`SPG_M3\`, \`SPG_M4\`, \`SPG_M5\`, \`GPPVA_DD2\`, \`GPPVA_DDME2\`, \`GPPVA_FSU2\`, \`GPPVA_FSU2H\`, \`GPPVA_FSU2R\`, \`GPPVA_NL3wrL55\`, \`GPPVA_TM1e\`, \`GPPVA_TW\` |
| Core models | 5 | \`APR_APR\`, \`BBB_BHF2\`, \`GM_GM1\`, \`RG_SLY4\`, \`RG_KDE0v1\` |
| Smooth EOS390 method IDs | 6 | `first_law_constrained_hermite_pressure_bridge`, `css`, `cubic_F`, `cubic_P`, `quintic_F`, `quintic_P` |
| Separate sharp-transition method ID | 1 | `maxwell` (existing 65-pair sharp-Maxwell workflow; not part of the 390 smooth matrix) |
| **Total** | **390** | \`13 × 5 × 6\` |

The recorded v3.2 baseline contains 390 constructed/read-back combinations, with 234 whole-table hard-gate PASS and 156 diagnostic-only whole-table C10 causality failures. The 156 failures are associated with the high-density source tails of \`APR_APR\` and \`RG_KDE0v1\` (78 pair-method combinations each). The historic causal-domain subrange count does not convert a full-table C10 failure into PASS. These are recorded historical values, not a new independent execution from this session.

Never overwrite or “clean up”:

- \`/ns_eos/matched_eos_all_combinations\`
- \`/ns_eos/Crust_Core_Matching_390_Summary.csv\`
- any original/raw EOS input directory
- the v3.2 release archive, checksums, or per-model papers
- prior attempt outputs or receipts.

New generated data must be written under a new immutable attempt path, for example:

\`\`\`text
<cloud-workspace>/artifacts/eos390/<UTC_RUN_ID>/
\`\`\`

The path is relative to the cloud workspace. Do not use the previous local Dropbox absolute path in a cloud job.

## 2. Branch and cloud bootstrap contract

The active development branch is \`cloud/eos390-fortran-unification-v1\`. Its parent is the reviewed E390 feature head \`b12ba451b6755bf8f13d22c882075a23eae800a7\`, not \`main\`. The branch is intentionally separate from both \`main\` and the existing feature branch. Do not merge or push to \`main\`.

A cloud agent must begin from a clean shallow checkout:

\`\`\`bash
set -euo pipefail
git clone --depth 1 --branch cloud/eos390-fortran-unification-v1 \\
  https://github.com/cosmosapjw-quantum/NS_Cowling.git NS_Cowling
cd NS_Cowling
test "$(git branch --show-current)" = "cloud/eos390-fortran-unification-v1"
git status --porcelain
git rev-parse HEAD
git rev-parse HEAD^{tree}
git log -1 --format='%H%n%T%n%P%n%s'
\`\`\`

The agent must record the exact branch, commit SHA, tree SHA, parent SHA, clean/dirty state, compiler version, Python version, OS image, and command line before editing. If the branch or parent differs from the contract, stop and reconcile before continuing.

Expected cloud build and tests:

\`\`\`bash
set -euo pipefail
gfortran --version
python --version
python -m pip install --upgrade pip
python -m pip install -e '.[test]'
bash fortran_eos/build.sh build/compose_eos_fortran
PYTHONPATH=python_frontend \\
COMPOSE_EOS_FORTRAN_ROOT="$PWD" \\
COMPOSE_EOS_FORTRAN_EXE="$PWD/build/compose_eos_fortran/compose_eos_fortran" \\
python -m pytest -q python_frontend/tests pipeline/tests
\`\`\`

Also compile and run the native analytic unit tests exactly as the repository workflow does:

\`\`\`bash
gfortran -std=f2018 -O2 -ffree-line-length-none -fimplicit-none -fcheck=all \\
  -Ibuild/compose_eos_fortran/mod \\
  fortran_eos/tests/test_polynomial_extrema.f90 \\
  build/compose_eos_fortran/obj/compose_eos_types.o \\
  build/compose_eos_fortran/obj/compose_eos_reconstruct.o \\
  -o build/compose_eos_fortran/test_polynomial_extrema
build/compose_eos_fortran/test_polynomial_extrema

gfortran -std=f2018 -O2 -ffree-line-length-none -fimplicit-none -fcheck=all \\
  -Ibuild/compose_eos_fortran/mod \\
  fortran_eos/tests/test_css_segment.f90 \\
  build/compose_eos_fortran/obj/compose_eos_types.o \\
  build/compose_eos_fortran/obj/compose_eos_css.o \\
  -o build/compose_eos_fortran/test_css_segment
build/compose_eos_fortran/test_css_segment

sha256sum -c SHA256SUMS_v4.txt
\`\`\`

Never report a command as passed unless its actual return code and logs were inspected. A skipped, unavailable, or unexecuted step is \`NOT_RUN\`, not PASS.

## 3. Surviving source archive and reconstruction contract

Dropbox is not automatically mounted inside a clean cloud sandbox. Do not assume a cloud agent can call the user's Dropbox connector or access the user's local `/home/.../Dropbox` path. The GitHub `eos_generator/eos_data` tree is retained for source development, but normalized GitHub tables are not byte-identical provider inputs.

The previously required archive `NS_Cowling_EOS390_raw_inputs_v1.tar.zst` is **lost**. Do not ask the researcher to upload it again, and do not attempt to recreate its byte identity from extracted text, screenshots, normalized GitHub files, or inferred values. Its exact bytes and original checksum cannot be recovered from the evidence currently available.

A surviving replacement package has been identified in Dropbox:

- Archive: `EoS_complete_core_crust_release_v3_2.zip`
- Dropbox path: `/EoS_core_crust_release_v3_2_20261106/EoS_complete_core_crust_release_v3_2.zip`
- Size: **19,398,255 bytes**
- Dropbox revision: `65fba15b552db0af03d47`
- Dropbox content hash: `a019e13e21de29726aead8c7d75caddf6013432d59db730ec58848825f87ed26`
- SHA-256 stated by `DELIVERY_MANIFEST_v3_2.json`: `d671310eb630ffd13f1b278425d945741ad7e08ce14c0d9aff42309c707e6a3b`
- ZIP member count stated by the delivery manifest: **2,516**
- ZIP CRC test: **PASS as reported by the v3.2 delivery manifest**

The SHA-256, size, member count and CRC above are historical Dropbox/manifest evidence. The ZIP bytes were not staged into this cloud workspace during the inventory operation, so the SHA-256 has **not** been independently recomputed in this turn. The cloud agent must recompute it before admitting the source cohort.

The checked-in inventory is `data/eos390_reconstructed_source_inventory.v1.json`. It records the 18 EOS models and all **158 files** in the surviving `code_package/eos_data/` subtree (3,883,275 bytes in total), including each relative filename and listed size. A second Dropbox listing of the v3.2 release's extracted `code_package/eos_data/` directory has exactly the same 158 relative paths and file sizes: no missing paths, extra paths, or size mismatches. This is a path-and-size comparison only, not a byte-identity proof. A sample text extraction of GPPVA_DD2 `eos.thermo` and `eos.compo` matched between the two Dropbox roots, but text extraction is explicitly not used as binary evidence.

Provider absences are intentional and must remain absent: APR_APR has no `eos.compo`; BBB_BHF2 has no `eos.compo` or `eos.mr`; GM_GM1 has no `eos.mr`. Never synthesize these files or encode missing composition as zero.

### 3.1 Cloud staging and reconstruction

The replacement ZIP, not the lost tar.zst, is the source input for the next cloud run. Stage it outside the Git worktree:

```bash
set -euo pipefail
mkdir -p "$WORKSPACE/inputs/replacement"
# Stage the exact Dropbox file:
# /EoS_core_crust_release_v3_2_20261106/EoS_complete_core_crust_release_v3_2.zip
ZIP="$WORKSPACE/inputs/replacement/EoS_complete_core_crust_release_v3_2.zip"

test "$(stat -c %s "$ZIP")" = "19398255"
echo "d671310eb630ffd13f1b278425d945741ad7e08ce14c0d9aff42309c707e6a3b  $ZIP" | sha256sum -c -
python - <<'PY' "$ZIP"
import sys, zipfile
p=sys.argv[1]
with zipfile.ZipFile(p) as z:
    assert len(z.infolist()) == 2516, len(z.infolist())
    bad=z.testzip()
    assert bad is None, bad
print("ZIP member count and CRC: PASS")
PY
python scripts/prepare_eos390_reconstructed_inputs.py \
  --archive "$ZIP" \
  --inventory data/eos390_reconstructed_source_inventory.v1.json \
  --output-root "$WORKSPACE/inputs/raw_compose" \
  --output-archive "$WORKSPACE/inputs/NS_Cowling_EOS390_reconstructed_inputs_v2.tar.zst"
```

The helper selects only `code_package/eos_data/<model>/<file>`, checks the exact inventory, rejects duplicates/missing/extra/unsafe/nested/symlink/size-mismatched members, preserves extracted file bytes, writes fresh `SHA256SUMS.txt` and `source_manifest.json`, and creates a newly versioned `NS_Cowling_EOS390_reconstructed_inputs_v2.tar.zst`. Its unit tests are in `python_frontend/tests/test_reconstructed_source_prep.py`.

After reconstruction, independently verify the new tar.zst listing, its SHA-256, every file hash, the 18 model directories, and the exact 158-file count. Stage the resulting archive or extracted tree into the cloud workspace's `inputs/raw_compose/` and attach the fresh source manifest to the M0 receipt. Do not overwrite any existing input, output, or historical release.

**Identity limitation:** the new tar.zst is a newly versioned reconstruction from the surviving v3.2 ZIP. It is not byte-identical to the lost `NS_Cowling_EOS390_raw_inputs_v1.tar.zst`; that equivalence is unprovable. Public GitHub should contain the inventory, provenance and reconstruction procedure, not raw EOS data bytes, unless redistribution rights are independently established.

**Current M0 state:** replacement archive identified and inventoried, but binary staging and independent rehash are pending. M0 remains BLOCKED until all cloud-side checks pass.

## 4. Single backend and typed data contract

Fortran owns all production scientific calculations. Python is restricted to:

- command argument and path validation
- input file staging and SHA-256 calculation
- serialization/deserialization
- calling the native Fortran executable
- decoding native JSON/CSV receipts
- presentation and summary generation without new numerical inference.

Python must not be the production implementation for EOS interpolation, EOS reconstruction, composition solving, root finding, branch selection, crust-core matching, causality, TOV, tidal, Cowling, or numerical batch orchestration. No silent Python/SciPy fallback is allowed. Historical Python modules remain read-only characterization oracles until a named native feature gate is passed.

### 4.1 CompOSE key and row rules

The unique row identity is the full integer triple:

\[
(i_T,i_{n_B},i_{Y_q}).
\]

Each individual CompOSE table must contain at most one record for any key. Duplicate keys must fail closed **before writing into the destination arrays**. A row outside the selected grid or a mismatched T/Yq slice must fail before using its indices.

A present numerical zero is not equivalent to missing data. File presence, row presence, and individual field presence are three distinct facts. Do not encode missing fields as zeros and later use a global Boolean to claim all species are available.

The limited current EOS data view is not the general CompOSE ontology. The generic reader must retain row-local fields before projecting them into a cold one-dimensional EOS.

### 4.2 Row-level composition and applicability

The internal raw representation must carry:

\`\`\`text
compose_cell:
  key:
    iT: integer
    inb: integer
    iYq: integer
  row_present: boolean
  Iphase:
    present: boolean
    value: integer
  particles:
    - particle_id: integer
      value: float64
      present: boolean
  nuclear_groups:
    - group_id: integer
      Abar: float64
      Zbar: float64
      YI: float64
      present: boolean
  thermo_extensions:
    - extension_id: provider-defined
      value: typed scalar or raw token
      present: boolean
\`\`\`

Do not turn an absent species into a measured zero. For each row evaluate C3 and C4 independently:

- C3 applicability requires the complete baryon abundance inventory needed for \(Y_n+Y_p+Y_h=1\).
- C4 applicability requires the proton, nuclear charge contribution, electron and muon terms needed for \(Y_p+(\bar Z/\bar A)Y_h=Y_e+Y_\mu\).
- Store the applicability and result per row, not in one table-global Boolean.

Per-row result states are exactly:

- \`CHECKED\`: all required values exist and the check was evaluated.
- \`PARTIAL\`: some required data exist but the check cannot be completed.
- \`NOT_AVAILABLE\`: the provider does not supply the required field and no approved derivation exists.
- \`FAIL\`: all required data exist and the numerical/physical closure fails.

Unknown particle IDs and multiple nuclear-group records must be retained as typed/raw records in the generic table model. If the current reduced EOS view cannot represent them, mark the projection unsupported with a precise reason. Do not silently discard or remap them.

### 4.3 Provider semantics are a separate claim gate

CompOSE v3.01 defines the standard field syntax, particle number fractions, and nuclear-group number abundance. For a group:

\[
Y_I=\frac{n_I}{n_B},\qquad Y_h=\bar A\,Y_I.
\]

For the supported one-group composition closure:

\[
Y_n+Y_p+Y_h=1,\qquad
Y_p+\frac{\bar Z}{\bar A}Y_h=Y_e+Y_\mu.
\]

This standard is not proof that a particular provider export actually encodes a disputed field with the expected meaning. The sampled SPG/GPPVA rows do not satisfy one consistent interpretation of their last group field across outer and inner crust. That provider-specific semantic problem is currently \`BLOCKED\`.

Do not choose between \(Y_I\), \(Y_h\), \(X_h\), a charge fraction, or a provider-specific quantity by making whichever formula pass. Preserve the source token and metadata and require a written model/provider definition plus independent closure. Until then, do not use the disputed field to construct a physical matched composition. This block does not prevent byte-preserving intake or metadata-only row-key tests.

### 4.4 Nadd and categorical extensions

Nadd counts additional \`eos.thermo\` columns. Their meanings are provider-specific. GPPVA_FSU2H uses an observed extension value that behaves as a categorical pasta-geometry label. It must not be numerically interpolated or confused with \`Iphase\`. Preserve the raw extension ID/value and provider declaration. If the extension cannot be represented across a combined output, stop or emit a named loss receipt; never silently set it to zero or treat it as a real-valued thermodynamic variable.

### 4.5 General axes and dimensions

The existing reader is deliberately narrow: one cold T=0, one-dimensional fixed-\((i_T,i_{Y_q})\) slice, a small particle-ID set, at most one nuclear group, and a limited Nadd layout. This is not universal CompOSE support. The M10 generic generator must retain the full T, \(n_B\), and \(Y_q\) axes and their exact bounds, support valid 1D/2D/3D layouts, and make slice selection explicit. No implicit flattening, axis averaging, nearest-index substitution, or temperature/charge-fraction guessing is permitted.

## 5. Physical conventions and numerical claims

Project EOS units:

- \(n_B\): \(\mathrm{fm}^{-3}\)
- \(P,\epsilon\): \(\mathrm{MeV\,fm}^{-3}\)
- \(\mu_B\): \(\mathrm{MeV}\)
- \(c_s^2/c^2\): dimensionless.

Spacetime calculations use metric signature \((- ,+,+,+)\). Do not silently change units or signature between modules.

For cold matter, keep the project master-function conventions:

\[
x=\ln n_B,\qquad F=\frac{\epsilon}{n_B},
\]
\[
Q_1=\frac{dF}{dx}=\frac{P}{n_B},\qquad
P=n_B Q_1,\qquad
\mu_B=F+Q_1=\frac{d\epsilon}{dn_B},
\]
\[
\frac{c_s^2}{c^2}=
\frac{Q_1+dQ_1/dx}{F+Q_1}.
\]

At \(T=0\), test the declared Q6=Q7 relation only where the source contract supports it. For \(I_l=0\), preserve the established charge-potential convention:

\[
\mu_{B,\mathrm{Euler}}=(Q_3+1)m_n+Y_q Q_4m_n.
\]

Do not use a barotropic \(\mu_B\) identity to erase the charge-potential correction when checking the source thermodynamic columns.

### 5.1 First-order sharp Maxwell interface

For two locally charge-neutral cold branches with one globally conserved baryon charge:

\[
P_L(\mu_t)=P_H(\mu_t)=P_t,
\quad
\mu_{B,L}(\mu_t)=\mu_{B,H}(\mu_t)=\mu_t,
\]
\[
n_L=\frac{dP_L}{d\mu_B},\quad
n_H=\frac{dP_H}{d\mu_B},\quad n_H>n_L>0,
\]
\[
\epsilon_i=\mu_t n_i-P_t,\quad
\Delta\epsilon=\mu_t\Delta n_B>0.
\]

Pressure and chemical potential are continuous; density and energy density have a positive finite jump. The two one-sided rows are adjacent. No derivative stencil, spline, interpolation interval, or root search may cross the interface. This is a locally neutral Maxwell construction, not global-neutrality Gibbs matter.

The historical v1 EOS390 token method="maxwell" named the smooth NPPHB pressure bridge. Active code now identifies that kernel as first_law_constrained_hermite_pressure_bridge. Existing v1 manifests, task receipts, output folders and hashes containing __maxwell are immutable historical provenance and must not be rewritten. The public method ID maxwell now means only a true sharp local-neutrality first-order coexistence construction. The common native matcher remains fail-closed for maxwell until composition-aware M5 integration is complete.

### 5.2 CSS high-density branch

A constant-sound-speed high-density branch is distinct from the interpolation name and from a phase-order classification. For a sharp interface with \(s=c_s^2/c^2\), high-side state \((P_t,\epsilon_+,n_+)\), and \(0<s\le1\):

\[
\epsilon(P)=\epsilon_+ + \frac{P-P_t}{s},
\]
\[
n(P)=n_+\left[1+
\frac{(1+s)(P-P_t)}{s(\epsilon_++P_t)}
\right]^{1/(1+s)}.
\]

With \(\mu_t=(\epsilon_++P_t)/n_+\), the branch obeys \(dP/d\epsilon=s\), \(d\epsilon/dn_B=\mu_B\), and \(dP/dn_B=s\mu_B\). The analytic form was symbolically checked with Wolfram Language: \(d\epsilon/dP=1/s\), \(dP/dn_B\;/\;\mu_B=s\), \(d\epsilon/dn_B\;/\;\mu_B=1\), and the stated \(\mu_B(P)\) expression has unit residual. This symbolic check does not substitute for a native test of the compiled kernel.

A zero jump does not establish a second-order transition. Use \`CONTINUOUS_UNCLASSIFIED\` unless a thermodynamic potential, derivative continuity, resolved nonanalytic higher derivative/susceptibility, and refinement robustness establish a stronger classification.

### 5.3 TOV, tidal and Cowling interface conditions

The native stellar solver must be used for production. A true sharp transition requires phase-segmented background integration and the appropriate slow/rapid conversion perturbation interface conditions. Tidal matching must use the finite-pressure jump term from Takátsy & Kovács (2020), arXiv:2007.01139, with its assumptions documented. Do not infer stellar stability from a local Seidov flag alone. The Cowling equations, surface boundary condition, central regularity, matching/interface conditions, and mode-shooting normalization must all be explicit and independently tested before M6 passes.

## 6. Single DAG and dependency order

The JSON DAG is the only machine-readable execution authority. The intended dependency structure is:

\`\`\`text
M0 source lock
  -> M1 legacy semantics and model bibliography
  -> M2 strict native CompOSE I/O and phase metadata
      -> M2A row/field presence masks and independent C3/C4 applicability
      -> M10 generic 1D/2D/3D CompOSE table model
M2 + M2A + M10 + provider semantic resolution
  -> M3 literature/provider composition solvers
  -> M4 five smooth matching methods
  -> M5 true sharp Maxwell matcher
M2 + M3 + M4 + M5 + M5A
  -> M6 native TOV / tidal / Cowling and interface validation
M6
  -> M7 exact manifest-driven 390 task runner
  -> M8 independent physics, differential, convergence and readback validation
  -> M9 immutable output release, checksum and backup
\`\`\`

The implementation task added in this cloud branch is the bounded parser-integrity slice E390-02B: reject duplicate \((i_T,i_{n_B},i_{Y_q})\) keys in \`eos.thermo\` and \`eos.compo\` before any destination array is overwritten. It does not resolve missingness, provider semantics, multidimensional tables, or the 390 output goal.

Current blockers remain explicit:

- E390-02A: per-row/per-species/per-group presence masks not yet implemented end-to-end.
- E390-10: general CompOSE axes and dimensions not yet implemented.
- C6B: SPG/GPPVA disputed provider abundance meaning not resolved by the manual alone.
- M3–M6: native composition solvers, complete matchers, and production stellar integration still need acceptance evidence.
- M7–M9: no new 390 regeneration is authorized yet.

## 7. Test and numerical verification contract

Every code change follows test-first development:

1. Read the exact source module, its callers, serialization types, and relevant tests.
2. Add the smallest regression test that demonstrates the defect.
3. Run it against the unmodified code and record the actual failure.
4. Implement the smallest scientifically neutral correction.
5. Run the targeted test, full \`python_frontend/tests\` and \`pipeline/tests\`, native analytic tests, and checksum verification.
6. Inspect logs and return codes. Record each result as \`PASS\`, \`FAIL\`, \`BLOCKED\`, or \`NOT_RUN\`.
7. A passing compile does not establish scientific validity; a passing regression does not establish that 18 sources or 390 outputs are correct.

For CompOSE round-trip and source-cohort validation, test:

- duplicate, missing, out-of-range, and shuffled row keys
- independent T/Yq selection and grid lower/upper bounds
- zero-valued present field versus absent field or absent row
- empty composition row, phase-only row, partial particle pairs, missing lepton species
- multiple particle IDs, multiple nuclear groups, unknown groups, non-default Iphase
- Nadd=0, Nadd=1 continuous extension, and categorical extension values
- exact Q1..Q7 and raw metadata preservation
- readback of every generated file, not merely existence or row count.

For every physics kernel, include analytic fixtures, dimension checks, monotonicity, limiting cases, boundary conditions, and independent residual calculations. Use scale-aware tolerances fixed before execution. Never tune a tolerance after observing a failure unless a separate mathematical/numerical justification is recorded.

## 8. Mutation gates

| Gate | Required evidence | Abort condition |
|---|---|---|
| G0 Branch lock | exact repo, branch, HEAD, tree SHA, clean worktree | wrong branch, moved/unexpected parent |
| G1 Source lock | file-level input/old-output SHA-256 manifest | hash mismatch, missing source, substituted normalized data |
| G2 Read-before-edit | full module, callers, tests, relevant DAG read | invented symbols or unreviewed authority |
| G3 Physics gate | equations, units, hypotheses, boundary conditions, cited/derived status | unsupported physical interpretation |
| G4 Build/static | compiler and Python versions, clean compile, exact tests, manifest | build/test failure or unexecuted test reported as pass |
| G5 Independent validation | analytic and independent-oracle results | residual mismatch, unexplained divergence |
| G6 Isolated output | immutable attempt directory and run ID | any write to source or historic output |
| G7 Readback | native readback plus per-task JSON and SHA-256 | missing row/file/check or malformed output |
| G8 Git mutation lease | expected branch SHA, commit parent, reviewable diff | force push, main mutation, unexpected branch head |
| G9 Publish | 390 manifest, statuses, readbacks, failures, hashes | incomplete task matrix or unresolved blockers |
| G10 Closeout | DAG status, provenance, anti-drift audit, next action | unsupported PASS or omitted failed task |

GitHub mutations must be serialized per file. Re-fetch the current blob SHA before updating a file. Re-fetch the branch head before commit/ref changes. Never force push. Do not delete or rewrite the history of the base branch.

## 9. Anti-drift / hallucination self-audit

The agent must answer each item with a file path, line/range, receipt, or hash. A bare assertion is not evidence.

1. What exact repository, branch, commit SHA, tree SHA and parent SHA did I inspect?
2. Was the tree clean before the change and clean after tests?
3. Did I read the full target implementation and all affected callers/tests before editing?
4. Is every physics equation cited to a primary source or labelled as a direct derivation?
5. Are units, metric signature, approximations, boundary conditions and interpolation domain explicit?
6. Did I conflate Maxwell, CSS, interpolation law, and transition order?
7. Did I infer physical phase meaning from an integer \`Iphase\` or a categorical \`q_extra\` value?
8. Did I turn missing composition into zero, discard an unknown species/group, or collapse row-level presence into a global Boolean?
9. Did I treat the normalized GitHub table as byte-identical to the Dropbox raw provider file?
10. Did any interpolation/root/derivative stencil cross an unresolved source boundary or sharp interface?
11. Did Python do any production numerical work that belongs in Fortran?
12. Did I modify a source EOS, historic v3.2 output, or previous attempt directory?
13. Does every 390 pair-method key appear exactly once in the task manifest, with a receipt even when it failed?
14. Does each result include input and output hashes, compiler and code version, parameters, physical gate results, and readback state?
15. Did I silently skip, retry-overwrite, smooth, clip, or relabel a failed EOS?
16. Which gates remain \`BLOCKED\` or \`NOT_RUN\`, and what exact evidence would unblock them?

Any unsupported answer is a blocker. Never fabricate citations, paper contents, source hashes, CI states, numerical results, or successful outputs.

## 10. Literature and source provenance

The bibliography JSON contains per-provider references. Important shared and model-specific anchors include:

- CompOSE Reference Manual v3.01, S. Typel et al. (2022), arXiv:2203.03209. Use it for file grammar and standard field meanings; it does not settle undocumented provider-specific export semantics.
- Pearson et al., MNRAS 481, 2994 (2018), DOI 10.1093/mnras/sty2413, for the BSk22 outer crust.
- Grill et al., PRC 90, 045803 (2014), arXiv:1404.2753, for the GPPVA Thomas-Fermi inner crust.
- Gögelein et al., PRC 77, 025802 (2008), DOI 10.1103/PhysRevC.77.025802, for the SPG density-dependent RMF functional.
- Scurto, Pais & Gulminelli, PRD 109, 103015 (2024), arXiv:2402.15548. The local Dropbox data sheets carry an earlier “in preparation” annotation; current CompOSE pages link this later publication to SPG variants. Verify exact variant-to-paper correspondence before using it to infer undocumented columns.
- APR: Akmal, Pandharipande & Ravenhall, PRC 58, 1804 (1998), arXiv:nucl-th/9804027.
- BHF2: Baldo, Bombaci & Burgio, A&A 328, 274 (1997), arXiv:astro-ph/9707277.
- GM1: Glendenning & Moszkowski, PRL 67, 2414 (1991), DOI 10.1103/PhysRevLett.67.2414.
- SLy4: Chabanat et al., NPA 635, 231 (1998), DOI 10.1016/S0375-9474(98)00180-8; Gulminelli & Raduta (2015), arXiv:1504.04493.
- KDE0v1: Agrawal, Shlomo & Au, PRC 72, 014310 (2005), arXiv:nucl-th/0505071; also Gulminelli & Raduta (2015).
- DD2: Typel et al., PRC 81, 015803 (2010), arXiv:0908.2344.
- Other GPPVA parameterizations (DDME2, FSU2, FSU2R/FSU2H, NL3wrL55, TM1e, TW) are in \`data/eos390_source_bibliography.v1.json\`; retain exact per-model mapping rather than treating the common GPPVA reference as a complete parameter citation.
- Tidal phase-interface jump: Takátsy & Kovács, PRD 102, 028501 (2020), arXiv:2007.01139.

Use public arXiv full texts when available and record the exact version. Do not commit publisher PDFs, paywalled full text, or third-party copyrighted course materials to this public repository. Cite source and DOI/arXiv in code comments and bibliography. If a non-open original is essential and not already accessible in the user's Dropbox, ask the user for that PDF rather than scraping or copying a restricted version.

## 11. Required run receipts and release contents

Each task receipt must include:

\`\`\`json
{
  "schema": "ns-cowling.eos390.task-receipt",
  "schema_version": "1.0.0",
  "run_id": "UTC_RUN_ID",
  "task_id": "CRUST__CORE__METHOD",
  "input_files": [{"path": "relative/path", "sha256": "..."}],
  "code": {"branch": "...", "commit": "...", "tree": "..."},
  "environment": {
    "os": "...",
    "compiler": "gfortran ...",
    "python": "...",
    "build_flags": ["..."]
  },
  "selection": {"iT": 1, "iYq": 1, "method": "..."},
  "status": "PASS | FAIL | DIAGNOSTIC_ONLY | BLOCKED",
  "checks": {"C1": "PASS", "C2": "NOT_APPLICABLE"},
  "readback": {"performed": false, "status": "NOT_RUN"},
  "output_files": [{"path": "...", "sha256": "..."}],
  "failure_reasons": [],
  "provenance": {"crust_model": "...", "core_model": "...", "references": []}
}
\`\`\`

The pair-method manifest must have exactly 390 unique IDs and no duplicate task keys. Every task, including failed and blocked ones, has a receipt. The run summary must count statuses by model, method, source family and physical check. Every table generated must be opened by the native reader and compared with the exact output manifest. Publish only the output from one immutable run ID, never a mixture of attempts.

## 12. Static audit and lean cloud branch

A source-level audit on the cloud branch found that the repository is not yet completely free of Python numerical calculations. Treat the following as explicit migration/containment issues rather than silently describing them as a pure Fortran production path:

| File / symbol | Current numerical behavior | Required status |
|---|---|---|
| `python_frontend/compose_eos/api.py::find_density_gap` | NumPy adjacent-density ratio and `argmax`; called by `fill_density_gap` before the native command | Move gap selection and phase-crossing check into Fortran; Python only sends threshold/options and decodes receipt |
| `python_frontend/compose_cowling/solver.py::count_nodes` | NumPy sign-change count on native solver output | Diagnostic post-processing only until native node counter exists; do not put in a production numerical claim without explicit review |
| `pipeline/compose_cowling/adapter.py` | Python PAV projection, Hermite interpolation, pressure integration and sampled causality checks | Diagnostic/reference-only; prohibit import from production EOS390 orchestration |
| `pipeline/compose_cowling/direct_barotrope.py` | SciPy PCHIP, Gauss-Legendre integration, inverse evaluation/root search | Diagnostic/reference-only; port approved interpolation to Fortran before production use |
| `pipeline/compose_cowling/intake.py` | NumPy first-law, monotonicity and sampled sound-speed checks | Diagnostic intake only; authoritative physical checks remain native Fortran |
| `python_frontend/compose_eos/eos_data.py` | NumPy arrays for typed data transport, serialization and copy operations | Allowed frontend representation; no physical solving/interpolation |

The Python adapter/direct-barotrope code is useful as an independent diagnostic or future differential oracle. Its existence does not mean it is the production EOS390 engine. Add import/call-graph tests that prevent these modules from entering the production Fortran pipeline. The remaining numerical selector in `find_density_gap` is small but real: do not waive it merely because the subsequent fill operation is native.

### Lean cloud tree

The cloud branch removes history-heavy payloads from the checkout while leaving the parent feature branch intact. After pruning, the recursive tree contained **1,996 blobs / 109,588,218 bytes (109.59 MB decimal)** versus **4,674 blobs / 360,093,413 bytes (360.09 MB)** before pruning: approximately **250.5 MB / 69.6% fewer tracked bytes**.

Removed from this branch's current tree:

- the old `oscillation_release_v1/` output/release snapshot
- `data/eos_cowling_adapter_v2/` and duplicate `DATA_ADAPTERS_V2.zip`
- the archived v3.2 EOS code-package ZIP (its separate checksum and provenance records remain where tracked)
- old `fortran_eos/differential_corpus/`, `_v2/`, `_v3/`, and per-model v4 differential payload directories.

Retained because current tests or the active source tree depend on them:

- `data/eos_original_normalized/` and `data/eos_cowling_adapter_v1/`
- `eos_generator/eos_data/` and `eos_generator/eos_data_normalized/`
- `fortran_eos/differential_corpus_v4/fortran_compose_differential.json`
- the active `oscillation_source/` Fortran solver and method source.

This is a current-tree trim only. There was **no Git history rewrite**. Removed content remains available from the parent feature branch history. Use a shallow clone for the online agent so old blobs are not fetched into its working checkout.

## Latest verified cloud baseline (2026-11-11)

The last verified parser/code baseline before the source-reconstruction helper was added is branch HEAD `cf89b3bb7ad5d1ee82e2f92b9736cadcdc582656`.

- GitHub Actions run [#43028463767](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/43028463767): **SUCCESS**; PR run [#43028475287](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/43028475287): **SUCCESS**.
- Environment: GNU Fortran 15.2; Python 3.13.
- Native Fortran build: PASS.
- Frontend and pipeline tests: **80 passed in 44.31 s**.
- Native quartic polynomial-extrema test: PASS.
- Native CSS high-density segment analytic test: PASS.
- Handoff byte diagnostic: PASS.
- Runner-measured SHA-256 manifest: PASS.
- Source snapshot packaged and uploaded as artifact ID `13860950517`, 55,105,617 bytes, SHA-256 `fa85568129cbc2b959735014982ee2b9a2c1a0b7b606b214f164358c0059126a`, expires `2026-12-11T13:23:24Z`.

This is a verified parser/code baseline, not an EOS390 production run. The replacement v3.2 ZIP has been identified and its 18-model/158-file metadata inventory is committed, but its binary has not been staged in the cloud workspace and its manifest SHA-256 has not been independently rehashed there. The next hard blocker is staging `EoS_complete_core_crust_release_v3_2.zip` from the Dropbox path documented in section 3, then producing `NS_Cowling_EOS390_reconstructed_inputs_v2.tar.zst` with fresh hashes. No new 390 generation, native readback of regenerated outputs, physical validation, or Dropbox EOS-output backup is claimed.

Duplicate CompOSE row-key guards are CI-verified. Per-row composition masks, generic multidimensional CompOSE support, provider-specific SPG/GPPVA abundance adjudication, and the complete Fortran production path remain incomplete and are tracked in the single machine DAG.

## 13. Current bounded implementation status

**CI-verified parser slice:** duplicate row-key guards in \`fortran_eos/src/compose_eos_io.f90\` reject repeated \((i_T,i_{n_B},i_{Y_q})\) keys for both \`eos.thermo\` and \`eos.compo\` before assignment. The phase index is now range-checked before indexing the destination phase array. Two regression tests were added. The first CI attempt compiled Fortran and reported 78 passed, 2 failed; both new tests failed because the test helper initially wrote literal backslash-n characters instead of line breaks. That helper was corrected in commit \`2265580a4bbdb204a4a4ea0b699d0f9ff5411838\`. A new CI run is required. Do not mark this slice PASS until the new run completes successfully.

**Established blocker:** the prior 18-source native verification receipt showed broad reader/semantic failures for 15/18 input tables, independent C3/C4 coverage problems for GM_GM1, and source causality failure in APR_APR. The source-set review also identified multidimensional table axes for RG_SLY4 and RG_KDE0v1. A prior 70-test/build receipt exists in Dropbox but is not a new run from this cloud branch.

**Not yet done:** row-level masks, generic 1D/2D/3D parser, provider semantics adjudication, native literature composition solver, all match methods, production sharp-Maxwell integration, final TOV/tidal/Cowling pipeline, 390 regeneration, readback, Dropbox backup.

## 14. Copyable cloud-agent handoff prompt

\`\`\`text
You are the authorized cloud coding agent for cosmosapjw-quantum/NS_Cowling.

GOAL
Regenerate the historical NS_Cowling 390 crust-core EOS pair-method combinations using a single Fortran numerical backend and a Python orchestration/serialization frontend. Preserve all source and historical outputs. Do not force failed EOS models to pass.

STARTUP
1. Shallow-clone exactly:
   https://github.com/cosmosapjw-quantum/NS_Cowling.git
   branch cloud/eos390-fortran-unification-v1
2. Verify branch, HEAD, tree SHA, parent SHA and a clean worktree.
3. Read these before editing:
   data/eos390_native_fortran_unification_dag.v1.json
   docs/EOS390_NATIVE_FORTRAN_UNIFICATION_CONTRACT_20261111.md
   data/eos390_source_bibliography.v1.json
   data/receipts/eos390_m0/E390-02_C6B_adjudication_20261111.json
4. Keep main and research/eos390-fortran-unification-v1 untouched.

MANDATORY INPUT
Do not use normalized GitHub EOS tables as byte-identical provider sources. The earlier NS_Cowling_EOS390_raw_inputs_v1.tar.zst was lost; do not request it again or claim to recreate its exact bytes.
Use the surviving Dropbox archive:
  /EoS_core_crust_release_v3_2_20261106/EoS_complete_core_crust_release_v3_2.zip
Read data/eos390_reconstructed_source_inventory.v1.json. Stage the ZIP outside Git at inputs/replacement/EoS_complete_core_crust_release_v3_2.zip. Independently verify size 19,398,255 bytes, SHA-256 d671310eb630ffd13f1b278425d945741ad7e08ce14c0d9aff42309c707e6a3b, 2,516 ZIP members and CRC. The hash is a value reported by the v3.2 delivery manifest and must be recomputed from the staged ZIP.
Then run:
  python scripts/prepare_eos390_reconstructed_inputs.py --archive inputs/replacement/EoS_complete_core_crust_release_v3_2.zip --inventory data/eos390_reconstructed_source_inventory.v1.json --output-root inputs/raw_compose --output-archive inputs/NS_Cowling_EOS390_reconstructed_inputs_v2.tar.zst
Require exact 18-model/158-file/3,883,275-byte match, fresh per-file SHA256SUMS.txt, source_manifest.json, and a verified tar.zst. Preserve provider absences and source bytes. This new archive is a reconstruction, not a byte-identical restoration of the lost tar.zst. If the replacement ZIP is not staged or any check fails, stop at M0 without requesting the lost filename.

IMPLEMENTATION ORDER
A. Verify E390-02B and finish E390-02C: migrate the remaining production-path Python numerical selector(s) to Fortran; keep Python/SciPy adapter kernels diagnostic-only.\nB. Finish E390-02A: add per-row/per-species/per-group presence masks; distinguish row absence, field absence, and numerical zero. Evaluate C3 and C4 independently per row. Preserve raw unknown IDs/groups and Nadd extensions.
B. Implement E390-10: generic CompOSE T, n_B, Yq axes and 1D/2D/3D data model; no implicit flattening or slice guessing.
C. Re-run original 18-provider cohort with exact source hashes and readback. Keep unsupported dimensions explicit.
D. Resolve SPG/GPPVA abundance-field semantics from written provider documentation and independent closure. The CompOSE manual is not sufficient evidence for undocumented provider conventions. Do not guess YI/Yh/Xh.
E. Port literature composition solvers for APR/BHF2 with source provenance and correct electron/muon threshold.
F. Port the five production smooth methods and the distinct FLCHPB legacy parity method; integrate the true sharp Maxwell method under the canonical ID maxwell. The historical v1 maxwell token is not retained as an active smooth alias.
G. Integrate phase-segmented Fortran TOV, tidal jump and slow/rapid Cowling interface conditions.
H. Construct exactly 390 tasks, each with an immutable machine receipt, failure state, source/output hashes and readback.
I. Only after all predecessor gates PASS, execute the 390 run in a new attempt directory and package the output artifact.

ENGINEERING RULES
- Fortran owns all production numerical calculations. No Python/SciPy fallback.
- Add failing tests before code. Do not weaken tolerances to clear failures.
- Preserve raw EOS files, provider metadata, phases, optional fields and categorical q_extra.
- Never interpolate Iphase or categorical extensions.
- No interpolation/root/derivative stencil may cross a sharp interface or unresolved source boundary.
- Acausal source tails remain failures or diagnostic-only. Never clip them to force PASS.
- Do not infer transition order from an interpolation method or phase integer.
- Never write to historic v3.2 outputs or previous attempts.
- Use a fresh attempt ID and output directory for every run.
- Re-fetch GitHub blob SHA and branch HEAD before every mutation; no force-push.
- Inspect CI logs and all return codes. A skipped test is NOT_RUN.
- Record provenance, exact equations, units, assumptions, boundary conditions and model-paper citations.
- Run the anti-drift questionnaire in section 9 before each commit and at closeout.
- Do not claim 390 regeneration or publication until all 390 task receipts, readbacks, checksum and independent validation reports exist.

CLOSEOUT
Return exact commit/tree SHA, run ID, build/test logs, source manifest hash, all 390 status counts, failed/blocked task catalogue, readback/checksum result, artifact URI and exact next action. If any source or physical semantics remain unresolved, stop at the corresponding DAG gate and preserve that failure.
\`\`\`
