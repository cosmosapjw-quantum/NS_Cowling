# Requested missing EOS recalculation and source-identified nested ZIP

Date: 2026-11-19 KST. Worktree: `research/eos-missing-data-recalculation-20261119`.
Scientific parent:
`7d0ff5088e87783f520a42d1f0ab7fad5ef1d380` (stage 05).
Requested scope: all seven classes of missing input and derived EOS shown in
the uploaded image headed **“❌ 없는 자료 (선배께 받아야 함)”**.

## Delivered artifact and precision of terms

The new machine-readable output is:

`data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip`

- Actual nested ZIP size: **4,413,310 bytes**.
- Actual SHA-256: **`299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04`**.
- Inner archives: **7**, each with an independent file-SHA manifest; outer/inner CRC tests passed.
- Fresh scoped regression after restoring immutable original provider EOS/PDF inputs in this isolated worktree: **153 passed, 3 skipped, exit code 0**.
- Bounds-checked GNU Fortran recompilation of the two native TOV/Cowling modules passed, with unused-argument warnings only; the 27 included source/data/ZIP hashes were unchanged over the test.
- Receipt: `data/eos_missing_data_recalc_20261119/verification/VERIFICATION_RECEIPT.json`; original ZIP bytes and PDF inputs remain untracked locally and are not separately republished to Git.

Four distinct concepts are **never conflated**:

- `original`: byte-preserved CompOSE/provider inputs or unmodified upstream raw data;
- `corrected`: a derived, explicitly energy-primary pressure/chemical-potential
  EOS on the source-selected *single inner-crust phase*, using native Fortran;
- `reconstructed`: independent candidate nuclear-theory equations solved by
  local CLD/TF or RMF Fortran generators; original author calibration need not
  have been recovered;
- `conditional_complete_stellar`: phenomenologically connected
  outer/inner/core EOS with actual 1.4 solar-mass Fortran DVODE TOV and Cowling.
  This is not microscopic functional unification.

The source and implementation origin are explicit in individual
`INTERNAL_PROVENANCE_AND_SHA256.json` manifests and in the top-level
`BUNDLE_MANIFEST.json`. All nested ZIPs were reopened in memory and subjected
to CRC and source-byte SHA-256 validation.

| Inner ZIP | Physical content | Actual action this session |
|---|---|---|
| 01_CORRECTED_13_SPG_GPPVA.zip | Original full 13 CompOSE crusts + original restricted single inner phase + 13 energy-primary repaired phases | **Regenerated all 13** using Fortran-native pressure and Euler chemical-potential equations; recovered exact original phase density grid and energy. All regenerated models were read back as CompOSE objects. |
| 02_RECONSTRUCTED_6_SPG5_FSU2.zip | SPG M1–M5 native CLD candidates + GPPVA(FSU2) TF | **Recomputed** 5 SPG 17-node stationary CLD tables and one FSU2 52-node native spatial TF table; also includes clearly labeled **previously computed 129-node SPG** comparison sources, not called fresh. |
| 03_PK1_CANONICAL_CORE_4097.zip | Published-parameter PK1 homogeneous beta-equilibrium RMF core, `0.015–1.0 fm^-3` | **Recomputed all 4097 points** with the native Fortran evaluator, physical-state validation, exact source-parameter provenance, and byte-for-byte equality against prior verified `eos.thermo`. |
| 04_BSK22_OUTER_SURFACE_LOWER.zip | Previous 217-node reconstructed BSk22 outer, 256-node newly evaluated anchored surface, analytic lower-inner 129-node source-fit and derivative variants | Re-evaluated source-anchored low-density surface and both analytic lower variants. **BSk22 microscopic ETFSI outer physics was not regenerated.** |
| 05_SPG5_PK1_10_WHOLE_STARS.zip | Five SPG source 129-node candidates + fresh PK1 4097-node core, separate Maxwell and CSS, reconstructed E BSk22 outer and causal bridge | **Newly computed 10 complete conditional stars** including `1.4 M_sun` native TOV plus slow and rapid Cowling. 9/10 pass strict `f,p1` node identity; SPG M3 CSS slow `p1` remains partial. |
| 06_NUDA_23_ORIGINAL_MODELS_PLUS_D1MS_FULL.zip | Original `nucleardatapy 1.0.2` raw crust datasets | **24 byte-preserved raw files**, representing 23 named models plus one D1M* crust+core extension; upstream license CC BY-NC-ND 4.0 and raw hashes retained. **These EDFs were not physically recalculated.** |
| 07_PROVENANCE_INPUTS_CODE_LOGS.zip | Source catalog (95 records), prior finite-pair controls, exact original Stage05 summary, selected scripts, command logs, readback receipts | Reproduction input and provenance. |

### Native full-star results at equal mass

All ten sharp-interface matched stars have `M = 1.4 M_sun` on a locally
ascending mass branch. Radii are newly computed:

| Crust + PK1 core | Maxwell R1.4 (km) | CSS R1.4 (km) | Native spectral status |
|:--|--:|--:|:--|
| SPG(M1) | 14.42042678 | 14.41581426 | Both conversions, strict f/p1 |
| SPG(M2) | 14.42041957 | 14.42030214 | Both conversions, strict f/p1 |
| SPG(M3) | 14.42043034 | 14.41569126 | **CSS slow p1 node unresolved**, other spectra accepted |
| SPG(M4) | 14.42043044 | 14.42592784 | Both conversions, strict f/p1 |
| SPG(M5) | 14.42043044 | 14.41795827 | Both conversions, strict f/p1 |

The native Fortran TOV/Cowling solver uses the previously declared
metric signature `(-,+,+,+)`, pressure and energy in MeV/fm³, density in
fm⁻³, mass in solar masses, radius in km, frequency in Hz, with distinct
slow and rapid phase-conversion junction conditions. `p1` is not
inferred from frequency order without the requisite node test.
These stellar models' crust core families are **not same-functional**;
a finite Maxwell crossing is an effective piecewise EOS construction,
not a unique microscopic phase transition.

### Input origin and missing author data

The original SPG surface coefficients `sigma0,b` were not recovered from
the source author and this archive does **not** solve the 2408-mass author
fitting process. The SPG 17-node outputs are local calibrated-candidate CLD
numerical solutions, not "original SPG" data. The generated FSU2 52-node
spatial TF output is independently computed but the CompOSE GPPVA provider
may differ; separate official source entries are delivered unaltered.

The BSk22 outer 217 nodes in the whole-star source come from a previously
reconstructed E model, **not** direct original BSk22 ETFSI calculations.
An analytic BSk22 lower 129-node model is included as **two alternatives**
`source_fit` and `energy_derivative`; they must not be merged without
an explicit thermodynamically admissible join. The attached surface and
lower bridges are explicitly phenomenological and non-original.

The 23 nominal external raw crust models from nucleardatapy 1.0.2 comprise
4 BSk 2018, 15 GMRS 2022, 3 Gogny 2020 (D1S/D1M/D1M*) and 1
Negele–Vautherin 1973. The fourth Gogny file `D1MS-full` includes the
same D1M* functional's core continuation. The raw files are covered by
the original *CC BY-NC-ND 4.0* licence. They are included unchanged with
attribution, but no Gogny/GMRS microscopic EDF solver was rerun, no missing
outer-crust continuation was invented, and no high-density extrapolation
is implicitly approved. Any public redistribution or commercial use must
respect the upstream licence.

### Reproducing this bundle

The explicit producer scripts and their Fortran/Python dependencies are
in the original repository. The original CompOSE input tables and
`nucleardatapy-1.0.2` wheel must exist at the source-hashed locations
recorded in each manifest. For each type:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend
# 13 corrected tables: Fortran source selected from each original CompOSE
python3 scripts/reconstruct_public_thermo_20261117.py INPUT_PHASE OUTPUT --model SPG_M1
# CLD candidate regeneration
python3 scripts/generate_spg_crust.py OUTPUT_DIRECTORY --points 17
# PK1 uniform RMF: see scripts/run_dense_rmf_20261118.py for the 4097 reference
# GPPVA TF: python -m compose_eos.cli generate GPPVA_FSU2 OUTPUT ...
# Exact source-anchored BSk22 analytic lower pair
python3 scripts/build_bsk22_lower_inner.py OUTPUT_DIRECTORY --points 129
# Full 10 star mode computation using fresh PK1 source:
python3 scripts/recalculate_missing_spg_pk1_stars_20261119.py
# Build all seven nested archives, manifest and master:
python3 scripts/package_nested_eos_missing_20261119.py
# Reopen and check nested CRC, source hashes, compile Fortran, test regression:
python3 scripts/verify_missing_eos_nested_20261119.py
```

Caution: the run above uses original input paths, parameter SHA and the
existing native binaries. It is an executable outline, not a self-contained
standalone build script that reconstructs unavailable original author
calibration. The detailed exact 52-density GPPVA(FSU2) invocation and input
parameter vector are in the corresponding ZIP receipts.

### Acceptance and stop condition

Scientific success for this artifact is bounded: **13 corrected + 6 newly
recalculated crust sources + native PK1 4097 + BSk22 two published analytic
alternatives + 10 native stars (9 strict full spectral, 1 partial) + 24
unmodified upstream raw data files**. Counts are not 13 + 6 + 1 + 10 + 24
different microphysical unified EOSs. The source origin and physical
uncertainties are part of the deliverable.

File-level SHA-256 and ZIP CRC check results are in the manifest, with
fresh Fortran/Python regression test receipts alongside the downloadable
master. No `main` mutation or PR is implied by creation of this branch.
