# Neutron-star nonradial oscillation follow-up: final checkpoint

> Scope note (2026-11-08): this is the historical v2 234-EOS f/p1 Cowling-spectrum checkpoint. COMPLETE applies only to that cohort, not active v4 EOS work. No new stellar-mode integrations are claimed by the v4 EOS changes.


**Status: COMPLETE for the declared calculation scope. No f/p1 integrations were rerun during resumption.**

## Source and cohort
- Input: /home/cosmosapjw/Dropbox/Neutron_Star_Oscillations_Seminar.zip
- Input SHA-256: b2c0dcaf6a5b07956cf8d33c39651d4ade32f44b4c49b33aa08fd84fb458185b
- Matching matrix: 390 candidates = 13 crusts x 5 core candidates x 6 matching prescriptions.
- Admitted: 234 source/join cases (78 each for BBB_BHF2, RG_SLY4, GM_GM1).
- Excluded: 156 source-core acausality cases (78 APR_APR and 78 RG_KDE0v1).

## Calculation
- ell=2 polar-fluid modes in the relativistic Cowling approximation.
- Cold beta-equilibrated barotropic perturbation closure; fixed spherical, nonrotating perfect-fluid background.
- Finite positive EOS pressure cutoff; no extrapolation to P=0.
- Target mass: 2.644606257493485e30 kg = 1.3300106263127078 nominal solar masses.
- Two-sided shooting/matching. f identity: W,V nodes=(0,0); p1 identity: W,V nodes=(1,1).
- All 234 unique EOS results are present and pass the declared node-identity and native projected-EOS cell audits.
- Parallel recovery already completed: 99 targeted reruns (27 previously incomplete + 72 GM_GM1 p1 identity corrections) with 12 workers.
- Representative convergence already completed: 3 EOS x 3 tolerance levels with 9 workers.

## Aggregate numerical diagnostics
- Maximum relative target-mass residual: 9.914291609902648e-13.
- Maximum absolute normalized matching residual: f 4.560284933008418e-08; p1 3.0891943864075344e-10.
- Maximum inward/outward matching-radius mismatch: f 7.275957614183426e-12 m; p1 9.094947017729282e-12 m.
- Maximum relative matching-mass mismatch: 6.386016275503321e-16.
- Density knots preserved by thermodynamic projection to roundoff; maximum reconstructed pressure change 1.2930514536147284e-3 P_max.
- All 234 projected EOS tables pass the native interpolation-cell audit.
- Representative loose-to-tight maximum frequency changes: f 9.35e-8 Hz; p1 1.71e-6 Hz. This is not a global error bound.

## Median frequencies and full ranges
- BBB_BHF2: f 2.5051 kHz [2.4883, 2.5136]; p1 7.0083 kHz [6.6076, 7.1823]; R 11.2684 km [11.2653, 11.2738].
- RG_SLY4: f 2.3698 kHz [2.3574, 2.3803]; p1 6.7177 kHz [6.3471, 6.8906]; R 11.7706 km [11.7670, 11.7745].
- GM_GM1: f 1.9202 kHz [1.8893, 1.9401]; p1 6.0090 kHz [5.6866, 6.0773]; R 13.8079 km [13.7980, 13.8215].

## Scope and caveats
- These are Cowling fluid frequencies, not full-GR quasinormal modes; no metric perturbations, gravitational-wave damping, rotation, magnetic field, elastic crust, superfluidity, frozen-composition response, or sharp phase-interface transfer modes are included.
- The maxwell suffix in the source matching matrix is a legacy smooth-matching label, not a sharp interface in the oscillation equations.
- The pressure/energy reconstruction is an explicit derived interpolation representation; original normalized EOS tables are preserved separately.
- The surface is a finite positive pressure cutoff, not an exact vacuum surface.
- The REVTeX author field remains a placeholder; replace it before submission. Two minor overfull hbox warnings are recorded in the validation report.

## Deliverables
- Release folder: /tmp/ns_oscillation_followup_20261106/DELIVERY_NS_COWLING_F_P1_v1
- Data archive: DELIVERY_NS_COWLING_F_P1_v1/NS_Cowling_f_p1_data_tables_v1.zip
- REVTeX PDF: DELIVERY_NS_COWLING_F_P1_v1/NS_Cowling_f_p1_REVTeX.pdf
- Full archive: /tmp/ns_oscillation_followup_20261106/NS_Cowling_f_p1_full_release_v1.zip
- Full archive SHA-256: 47ed76103ec510262ffa082dc8e78d9990f6511d623f1a43200bf4738ca4c83d
- Dropbox-synced folder: /home/cosmosapjw/Dropbox/NS_Cowling_f_p1_20261106
- Downloads copies: /home/cosmosapjw/Downloads/NS_Cowling_f_p1_data_tables_v1.zip, NS_Cowling_f_p1_REVTeX.pdf, NS_Cowling_f_p1_full_release_v1.zip

## Resume verification
- Original input and provenance copy have matching SHA-256.
- Archive contains 234 original normalized EOS CSVs, 234 projected EOS CSVs, and 234 projection receipts.
- Summary has 234 unique rows; 78 per core family; mode node identities and projected-EOS audit all PASS.
- Data and full-release ZIP tests pass.
- REVTeX PDF has 5 pages; pdfinfo and pdftotext checks pass.
- All 1,014 delivery-manifest entries verify.
- Dropbox local sync client reports Up to date; synced-folder and Downloads key-file checksums pass.
- BibTeX validation: 5 entries, 0 issues. Abstract logic review: 0 issues.

**Next action:** no f/p1 integrations remain for this v1 scope. Full-GR perturbations, frozen-composition closure, sharp-interface modes, or a zero-pressure surface are separate physical models and should use separately versioned outputs.
