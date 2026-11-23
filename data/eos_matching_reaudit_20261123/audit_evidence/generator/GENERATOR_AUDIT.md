# EOS generator / source-projection audit — 2026-11-23

Baseline: e51dadd132125f0c5389da8c92c50b84eeb524e8 in `/home/cosmosapjw/ns_cowling_wt/cowling_c20_native_tov_pchip_20261123`; integration branch `research/eos-core-crust-all-methods-20261123`. This report covers the generator frontend, public native core domain, source census, and native raw-table projection. It does not assert that every imported source or resulting match is physically admissible. Parent integration owns the final all-pairs campaign, backup and publication. No commits were made by this audit worker.

## Confirmed defects and corrections

1. **High — default executable freshness was bypassed.** `python_frontend/compose_eos/native.py::executable()` formerly returned an existing path without calling the content-signature build gate, and accepted a non-executable existing file. Reproduced with the public API. The default source build now goes through `build.sh --ensure`, which already owns SHA/compiler/flags signatures and cross-process locking. A thread-safe in-process metadata/environment cache avoids shell/build work for unchanged repeated calls; source, compiler, executable, build script or build-environment changes invalidate it. An explicit executable override must be an executable file. This change adds orchestration only, not numerical Python.

2. **High — native core APIs accepted densities outside their own declared numerical model domains.** Original probes returned `ok` for SLy4 n=2.0 fm^-3, SLy4 bulk [1.5,1.6,1.7], and GM1 n=.001. `solve_core_state` and `generate_core_native` now enforce finite positive density within the exact `core_model_bounds`. Declared domains remain APR [.085,1.34], BBB [.08,1.5], GM1 [.02,1.0], KDE0v1 [.02,1.8], SLy4 [.02,1.4] fm^-3; default sample intervals were not confused with supported domains. All ten exact domain endpoints were tested.

3. **Test fixture defect exposed by the domain fix.** `fortran_eos/tests/core_native_probe.f90` generated a geometric SLy4 endpoint 1.4000000000000001 for a request of 1.4. It now pins requested first/last endpoints, as the production generator already did. No physical bound or test assertion was weakened.

4. **New importer normalization-argument parser hardened.** A Fortran scalar list-directed read silently accepted `939 extra`, `939,938`, `939 938`, `1*939`, and an invalid proton-mass suffix. These five observed failures now reject before any output directory exists. Normalization arguments must each contain exactly one finite positive numeric token.

## Actual verification evidence

- `RED.log`: before generator/domain corrections, 23 failed / 5 passed, exit 1.
- `GREEN.log`: dedicated generator/domain suite, 28 passed, exit 0.
- `FOCUSED_REGRESSION.log`: first broader regression, 1 failed / 48 passed; the remaining failure was the endpoint fixture defect above.
- `FOCUSED_REGRESSION_FINAL.log`: 49 passed in 3.74 s, exit 0, running the dedicated suite plus `test_core_native_backend.py` and `test_core_support_domain_audit.py`.
- `IMPORTER_RED.log`: new CLI feature tests against explicitly unimplemented Fortran entry point, 3 failed. This is new-feature RED evidence, not a pre-existing importer defect claim.
- `IMPORTER_GREEN.log`: first full importer suite, 17 passed.
- `IMPORTER_MASS_RED.log`: five scalar parser adversarial cases failed before hardening.
- `IMPORTER_GREEN_FINAL.log`: 22 passed in 0.70 s, exit 0. Test build uses Fortran bounds checks and invalid/zero/overflow traps. Three analytic layout tests independently check epsilon=m*n+100*n^2, P=100*n^2 and mu=m+200*n; malformed values, duplicate/zero density, wrong field counts, insufficient rows, invalid masses, unknown layout and output reuse are rejected.

This is implementation verification of the specified contracts; it is not a continuum-limit proof, global beta-equilibrium minimum proof, or blanket source-physics acceptance.

## Source census and identity

Starting authority `data/all_origin_maxwell_css_20261118/SOURCE_CATALOG.json` contains 95 source records: 25 core and 70 crust. It contains 24 available core snapshots plus one blocked uploaded OMEG1 source, and 44 available crust snapshots plus 18 known raw P/E/n tables, four originally unresolved BSk tables and four remaining blocked/incomplete crust records. All 68 available CompOSE source paths were present in this checkout.

The old `enumerate_pairs` grouped pairs by origin tier, and the old runner reused sealed legacy Maxwell results and a cubic prefilter. Those are historical scopes, not a full fresh Cartesian all-method campaign. The new user authorization calls for all source identities across the full Cartesian product; original, corrected and reconstructed identities must remain distinguishable, as must blocked input rows.

New distinct post-catalog core source snapshots found in `data/cowling_maxwell_native_20261119` are generated DDME2, FSU2 and TM1e plus their three original provider homogeneous-core sources. FSU2/TM1e generated snapshots have historical Python stationary-RMF oracle provenance; honest reuse as historical input must not be described as a newly implemented Fortran generator. DDME2 is Fortran-native. Newly sampled SPG17/129, PK14097 and other 4097-point core snapshots are resolution variants rather than new nuclear functionals. Parent integration plans 101 source identities after the six distinct additions.

## Native raw-table import contract

New source: `fortran_eos/app/compose_eos_import_cold_table.f90`.

```
compose_eos_import_cold_table LAYOUT RAW_FILE NEW_OUTPUT_DIRECTORY NEUTRON_OFFSET_MEV PROTON_HEADER_MEV
```

The output parent must exist, and the output directory must not exist. Directory creation is atomic and refuses duplicate-worker overwrite. Raw data must have increasing positive density and exactly the declared field count. Blank and `#` comment lines are allowed. Column guessing, interpolation, inserted zero-density points, smoothing, sorting, source energy shifts and repairs are absent.

| Layout | Raw fields | 1-based n / E-offset / P | n unit | E-offset unit | P unit | Rest-energy offset |
|---|---:|---|---|---|---|---:|
| gmrs_2022 | 11 | 1 / 8 / 9 | fm^-3 | MeV/baryon | MeV fm^-3 | 939.565346 MeV |
| gogny_2020 | 9 | 1 / 5 / 6 | fm^-3 | MeV/baryon | MeV fm^-3 | 939.0 MeV |
| bsk_2018 | 28 | 1 / 10 / 8 | fm^-3 | MeV/baryon | MeV fm^-3 | 939.5651828 MeV |

Fortran computes epsilon=n*(m_n+E_offset), mu=m_n+E_offset+P/n and the exact corresponding CompOSE normalized columns. It writes CompOSE, reads it with the existing native reader, checks density identity and normalized thermodynamic roundtrip at scaled tolerance 2e-14, and then writes native CSVs. It does not infer composition or phases. The source Euler identity used to compute mu is a normalization definition; independent differential first-law and causality checks remain necessary.

Outputs: eos.nb, eos.t, eos.yq, eos.thermo (no eos.compo), eos_data.csv, eos_data.nml, normalized_eos.csv, SOURCE_VALUES.csv and IMPORT.json. `IMPORT.json` explicitly sets `physical_acceptance_asserted=false`, `composition_available=false`, `phase_labels_reconstructed=false`.

Frozen manifest: `NUDA_NATIVE_IMPORT_LAYOUTS.json`, 22 records / 22,908 rows, each with raw-byte SHA, wheel path, explicit columns and normalization provenance. It includes 15 GMRS, 3 Gogny and 4 BSk tables. Pinned wheel SHA256 is `6d8f57cf056696f39ee45fe4d516bc4eae7a2461c9aba5d83c2012b499116d61`.

Gogny +939.0 is author-verified by Mondal et al., arXiv:2003.03338v3 Table II caption and Fig.2 text: https://arxiv.org/html/2003.03338v3. BSk +939.5651828 and proton header 938.2718440 are verified from all four primary CompOSE eos.thermo headers; examples: https://compose.obspm.fr/download/1D/NS/Skyrme/BSK24/eos.thermo and https://arxiv.org/html/1903.04981. BSk column 10 is the ETFSI energy offset; column 11 is cgs mass density, despite the misleading NUDA loader name. Pinned NUDA BSk24/25 rows differ slightly from current provider revisions and remain separate source identities. Upstream NUDA spline-derived pressures with a fabricated (0,0) point are deliberately not used.

## Reproducible importer build

```
gfortran -std=f2018 -O2 -fcheck=all -ffree-line-length-none \
  -J/tmp/ns_eos_reaudit_20261123/generator/importer_build \
  -I/tmp/ns_eos_reaudit_20261123/generator/importer_build \
  fortran_eos/src/compose_eos_types.f90 fortran_eos/src/compose_eos_io.f90 \
  fortran_eos/app/compose_eos_import_cold_table.f90 \
  -o /tmp/ns_eos_reaudit_20261123/generator/importer_build/compose_eos_import_cold_table
```

Frozen importer executable SHA256: `b6f5edd5fe9bd9847a1741656619923afdf0245ea37b2a1facab9cb5942a39b9`.
Frozen source SHA256: `f4060b976722cfa91876e0c0d41cebe0ad29467bfa264e8f3ce0d61a96ec6c32`.
Frozen layout-manifest SHA256: `320958a27f9f5f9a86b8347c92d6f8eb80ecc1a973e1a6709c40d2a09ea61320`.

## Owned change list

- python_frontend/compose_eos/native.py
- fortran_eos/src/compose_eos_core_native.f90
- fortran_eos/tests/core_native_probe.f90
- python_frontend/tests/test_eos_reaudit_generator_20261123.py
- fortran_eos/app/compose_eos_import_cold_table.f90
- python_frontend/tests/test_cold_table_importer_20261123.py

The shared working tree contains other agents' matching/verification changes. They are not claimed as this worker's edits or validation results. Raw catalog/source tables were not modified.

## Addendum: recovered user total-energy source layouts

Two additional explicitly declared layouts were added after the initial NUDA freeze. `user_core_energy_density` has 7 fields with n/epsilon/P in columns 1/2/3; `user_crust_energy_density` has 5 fields with n/epsilon/P in columns 2/3/4. Epsilon already includes all rest energy and has units MeV fm^-3. The supplied neutron/proton masses set the CompOSE encoding header only. No rest energy is added: the receipt records `rest_mass_offset_mev=0` and `neutron_header_mass_mev=MN`.

The CLI is backward-compatible with five arguments and additionally accepts seven arguments:

```
EXE LAYOUT RAW_FILE NEW_OUTPUT_DIRECTORY MN_HEADER_MEV MP_HEADER_MEV [NMIN NMAX]
```

Without explicit bounds every original row is retained. With bounds the selection is inclusive and uses unmodified source density values. NMIN=0 is permitted; both bounds must be finite, single numeric tokens and 0 <= NMIN <= NMAX. Every raw row is checked for numeric validity, exact field count, positive density and strict density order before selection; a selected projection must contain at least three rows. Selection does not sort, interpolate or repair anything. The project-adopted user domains are core n>=.079 and crust n<=.06757 fm^-3, supplied explicitly by orchestration, not hardcoded into the importer.

Receipt schema is now `ns_cowling.cold_raw_projection.v2`; status remains `NATIVE_RAW_PROJECTION_READBACK_VERIFIED`. `rows` means selected output rows. Additional fields preserve raw_rows, selected_rows, excluded_rows, excluded_lower_rows, excluded_upper_rows, original source-row span, density_selection_applied, and raw/requested/selected density domains. `SOURCE_VALUES.csv` preserves every raw n/energy/P row, including explicitly excluded rows; the native normalized EOS contains only selected rows.

Actual verification: `IMPORTER_USER_ENERGY_RED.log` has eight observed new-feature failures before implementation; `IMPORTER_USER_ENERGY_GREEN.log` has **37 passed in 1.14 s**, exit 0, including all previous 22 importer tests. New independent analytic tests use total epsilon=939*n+100*n^2 and P=100*n^2, confirm physical epsilon/P/mu remain unchanged under distinct CompOSE header masses 939 and 900, exercise both layouts, inclusive endpoints, zero lower bound, exact provenance counts/domains, invalid windows and malformed rows outside a requested selection.

Latest executable SHA256: `2daaa16b93e4704a09b3cafa75164dac96539d48893502e313ddf5359097946d`.
Latest source SHA256: `9b9d78d85fa6c3af2189a4c3bde7cb2f2047861aa44be751aacaf32412de69df`.
The previous NUDA-only executable is preserved as `importer_build/compose_eos_import_cold_table_nuda22_verified`. The 22-source NUDA layout manifest remains byte-identical, SHA256 `320958a27f9f5f9a86b8347c92d6f8eb80ecc1a973e1a6709c40d2a09ea61320`.

GMRS normalization cross-check received from the physics auditor: primary provider H2 header `939.56535 938.27201 1` at https://compose.obspm.fr/download/1D/NS/Skyrme/H2/eos.thermo agrees with the pinned NUDA constants at the provider's printed precision. The last digits of the imported constant remain explicitly NUDA provenance, not inferred extra precision from the rounded provider header.
