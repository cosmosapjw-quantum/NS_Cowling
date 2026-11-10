# Sharp Maxwell C6: Native Execution Verification (v1)

Status: **C6_PARTIAL**. Evidence: \`artifacts/sharp_maxwell/C6_NATIVE_EXECUTION_20261110T093714Z/\`.

## Scope

Execution check of the C5 read-only native probe on SPG_M1 and GPPVA_FSU2. This establishes executable-path and serialization behavior only where explicitly stated. It is not a Maxwell certification, pair admission, root search, TOV or Cowling run.

## Identity and provenance

- Repository: \`cosmosapjw-quantum/NS_Cowling\`
- Branch: \`research/sharp-maxwell-cohort-admission-v1\`
- Base: \`research/maxwell-native-fortran-v1\`, SHA \`a09ec30d74b363652d2d1aefe990557c468679dd\`
- C5 IO commit: \`5ecbca8eec11ce17156d360d4a02376765bf055c\`
- **Execution head:** \`43eea8bd9bee7ebb222e3c66babf5078b398d7ae\`
- The branch later advanced to \`3f9531e6796a849b6f77b193bea2ecece74348f5\` and subsequent documentation/manifest-only corrections. The native execution must not be represented as having run on that later head.
- All five C5 blobs matched the earlier handoff.
- Toolchain: gfortran 13.3.0, Python 3.12.3.

## Execution results

| Gate | Result |
|---|---|
| Native build (official \`fortran_eos/build.sh\`) | PASS |
| Committed pytest, exact test file | PASS, 6/6 |
| SPG_M1 native probe | Native reader/writer/re-reader and verifier executed; runner PASS; round-trip PASS; scientific FAIL (C3, C4) |
| GPPVA_FSU2 native probe | Native reader/writer/re-reader and verifier executed; runner PASS; round-trip PASS; scientific FAIL (C3, C4) |
| Raw composition field round-trip | PASS as serialization/field preservation only; **not** physical semantic validation |
| Raw source SHA-256 before/after (18 files) | identical |
| Legacy outputs (234 files) before/after | identical |
| 65 pairs | unchanged, all BLOCKED_BEFORE_ROOT_SEARCH |

The earlier shorthand “Yh = YI mapping PASS” means only that the committed reader/writer path preserves the selected raw field under its current mapping and reproduces it after round-trip. It does **not** establish that the raw CompOSE field \(Y_I\) has the physical meaning required by verifier variable \`Yh\`.

## Composition semantics issue: BLOCKED_REVIEW_REQUIRED

### 1. What the standard defines

The [CompOSE Reference Manual v3.01](https://doi.org/10.1140/epja/s10050-022-00847-y) defines particle number fractions as \(Y_i=n_i/n_b\) in Eq. (3.10), and distinguishes particle mass-number fractions \(X_i=B_iY_i\) in Eq. (3.12). For a group of nuclei, Eqs. (4.20)-(4.22) define the group mean mass/charge numbers and combined group number fraction \(Y_I=\sum_jY_j\). Thus the standard baryon/mass fraction of the nuclear group is \(X_h=A_{\rm av}Y_I\).

The native verifier uses its internal \`Yh\` as a baryon/mass fraction:
- C3: \(Y_p+(Z_{\rm av}/A_{\rm av})Y_h=Y_e+Y_\mu\).
- C4: \(Y_n+Y_p+Y_h=1\).

### 2. What the provider data sheets establish

Both official data sheets identify \(Iphase=0\) as outer crust and \(Iphase=7\) as inner crust:
- [SPG(M1) data sheet](https://compose.obspm.fr/download/1D/Crust/M1/eos.pdf), lines 124-128.
- [GPPVA(FSU2) data sheet](https://compose.obspm.fr/download/1D/Crust/FSU2/eos.pdf), lines 59-63. It additionally identifies \(Iphase=6\) as homogeneous matter and separately describes the extra \`eos.thermo\` geometry code (0 outer crust, 1 droplets, 2 rods, 3 slabs, 4 tubes, 5 bubbles, 6 core).

Therefore, the phase label \(Iphase=7\) is not unknown. What remains unresolved is the meaning, normalization and consistency of the composition values in those records.

### 3. Observed raw rows and native checks

An outer-crust example in both selected raw files has \((A_{\rm av},Z_{\rm av},Y_I,Y_e)=(56,26,1,0.46429)\). Under the current direct field mapping \`Yh=YI\`, C3/C4 pass on these rows. Under the standard manual conversion \(X_h=A_{\rm av}Y_I\), the inferred value is 56, which is not a normalized mass fraction.

At the first inner-crust row (row 218):
- SPG_M1: \((A_{\rm av},Z_{\rm av},Y_I,Y_e)=(82.107,24.387,0.06087,0.06087)\).
- GPPVA_FSU2: \((740,36,0.04865,0.04865)\).

The listed free-neutron, free-proton and muon fractions are zero in those rows. Under the current \`Yh=YI\` mapping, C3/C4 fail on every inner-crust row:
- SPG_M1: 148/148 \(Iphase=7\) rows fail; C3 maximum absolute residual \(0.042790704812013594\); C4 maximum absolute residual \(0.9779100000000001\).
- GPPVA_FSU2: 52/52 \(Iphase=7\) rows fail; C3 maximum absolute residual \(0.04628324324324324\); C4 maximum absolute residual \(0.98856\).
- Every \(Iphase=0\) row passes C3/C4.

The raw \`eos.compo\` Git blob identities at source-authority commit \`a09ec30d74b363652d2d1aefe990557c468679dd\` are SPG_M1 \`ab79da517271bea1f7de37294e87e019de611699\` and GPPVA_FSU2 \`fa1ad2b809eb1c95825e2ee6eb2650cf3f248693\`.

These observations do not establish a replacement conversion. Blindly using \(A_{\rm av}Y_I\) produces values far above unity, while the current direct mapping fails the native composition checks on all inner-crust rows. The underlying issue may involve provider-specific output semantics, source composition validity/alignment, or the applicability of the verifier closure to these provider tables. These possibilities remain open.

### 4. Decision and restrictions

- No code, source data, C3/C4 equation or tolerance was changed in C6.
- The C5 statement that Eq. (4.22) makes \(Y_I\) itself the combined baryon/mass fraction is superseded as a semantic error.
- Current \`Yh=YI\` is recorded as the provisional raw-field-preservation implementation only; it is not an approved physical mapping.
- The issue remains \`BLOCKED_REVIEW_REQUIRED\`. C7 is not recommended until the composition semantics are adjudicated.
- Do not infer a general Iphase-dependent conversion from these two examples.
- Preserve the 65 pair records as \`BLOCKED_BEFORE_ROOT_SEARCH\`, legacy outputs as \`LEGACY_NON_MAXWELL_READ_ONLY\`, and root/TOV/Cowling status as \`NOT_RUN\`.

## Not run

\`sharp_transition_search\`, \`tov_cowling\`, C7-C10, and the 65-pair sweep.

## Reproduce

\`\`\`bash
bash fortran_eos/build.sh /tmp/c6build
export COMPOSE_EOS_FORTRAN_EXE=/tmp/c6build/compose_eos_fortran
python -m pytest -q python_frontend/tests/test_sharp_maxwell_source_probe.py
python scripts/sharp_maxwell_source_probe.py --manifest data/sharp_maxwell_65_cohort_manifest.v1.json --source-id SPG_M1 --source-id GPPVA_FSU2 --output out.json
\`\`\`

See the exact commands, exit codes, toolchain, hashes and native outputs in the evidence directory above.
