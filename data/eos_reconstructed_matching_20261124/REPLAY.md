# Reproducing the crust reconstruction and matching campaign

The source-code branch is based on b49d44bc0b2c774f4bd81ed6753f3795abddca6b in cosmosapjw-quantum/NS_Cowling. The code, regression tests, input contracts, source checks and HF plans are in the branch. The evidence archive holds the generated datasets and native outcomes. Paths in the commands below are relative to the repository root.

1. Generate BPS with the standard CLI. The source-contract README gives the independent CompOSE #68 and BPS Table 5 checks.

   ~~~sh
   PYTHONPATH=python_frontend python -m compose_eos.cli generate BPS_MTW1964_MS1966 /tmp/bps-new
   python data/eos_reconstructed_matching_20261124/check_bps_reconstruction.py --native /tmp/bps-new/generated --out /tmp/bps-check
   ~~~

2. Generate the two restricted profiles and the uniform beta branch with the native NV modules. The profile README describes the physical closure. Use the same command with 256 radial volumes to reproduce the mesh comparison. The 11 orbital cells, their charge neighborhoods, and the 96/160 mesh comparison have separate plans and instructions in orbital_hf/README.md.

   ~~~sh
   bash fortran_eos/run-nv-profile.sh /tmp/nv-profiles-new 128 49
   python scripts/audit_nv1972_skyrme_tables.py --help
   python scripts/audit_negele1970_radial_moments.py --help
   ~~~

3. The source-preserving outer/inner composite can be regenerated with the Fortran command in supplementary_bps_nv/README.md. It requires the 623-row BPS output and the 257-row NV1973_author_fit already present at the base commit. The controls at the first two physical NV densities and the analytic curvature and sound-speed checks are in bridge_verification.json.

4. Each lane uses the unchanged matching executable with SHA-256 d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588 and the archived core selections. See remaining_matches/README.md for full and public-only commands, physical-axis controls, prefix validation and CSV export. A reproducible public run can use the 31 public cores. The protected 32nd core requires the owner's original source catalog and selected source files. The numerical qualification concerns the finite admissible source views tested. It is not a scan over every possible interpolant or a stellar-structure calculation.

In the supplied evidence archive, data/eos_reconstructed_matching_20261124/ and data/eos_remaining_generated_20261124/ retain their repository-relative locations. The convenience paths in the Korean report identify their subdirectories. Independent evidence, source figures and the strict 672-row crosswalk are at archive root. The physical row ZIP uses one CSV per candidate and preserves its originating lane and identifier. Original journal scans and protected input tables are excluded.
