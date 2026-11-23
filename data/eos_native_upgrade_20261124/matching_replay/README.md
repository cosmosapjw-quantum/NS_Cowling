# Native matching replay

The current baseline was built from git commit `cab027e9b5ccbc42d6509531f7a6201f2c13d184`. A second native build used the same checkout with the new matching kernel. Both builds used GNU Fortran 14.2 and the same options, `-std=f2018 -O2 -ffree-line-length-none -fimplicit-none -fcheck=all -fopenmp`. The library link arguments were `-l:liblapack.so.3 -l:libblas.so.3`. The existing Cartesian process frontend sets all native worker thread counts to one. Comparisons here test numerical and categorical reproducibility. See the sibling matching benchmark for within-matcher OpenMP scaling.

The historical run used a distinct pinned executable (`d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588`). It has remained separate. The three comparison JSONs identify the executable and frozen input catalog by SHA-256 for each lane. Counts are aggregates across all 32 original core inputs, including the owner's protected source.

- `archive_to_parent.json` compares archived numerical outputs with a fresh run of the compiled parent.
- `parent_to_candidate.json` compares two fresh native builds for the original 672, the 224 BPS plus NV equilibrium-fit composite cases, and the seven legacy control cases.
- `archive_to_candidate_all.json` compares the freshly executed candidate with the archived results for 11 lanes plus the seven controls. The remaining seven scientific lanes in this comparison use the historical reference on the baseline side.

Each surviving cell is compared on all archived CompOSE members, both internal native EOS files, the complete native matching report, gate receipts, phase transport and respective source selections. The readback directory prefixes and provenance binary identities are expected to vary and are checked separately. Rejected view attempts preserve the complete native gate, reason and return code; elapsed time is excluded. Comparison requires unmodified input and binary hashes before and after each run and reads the per-cell archives and failed trial data. A separate corruption control altered one byte in the archived thermo table and recomputed its archive digest: the comparator rejected that cell.

Use the frozen source catalog directory produced by `scripts/run_remaining_reconstruction_matching.py`, which calls the unchanged Cartesian runner. An owner's full catalog must resolve its protected input locally. Public-only catalogs may also be replayed with the corresponding 31-core reference. Each pair of binaries under comparison uses the same catalogs and source bytes. The supplemental legacy control catalog consists of the previous NV1973 author fit and public 2020-MVCD-D1MS core.

```bash
REPLAY_WORKERS=3 bash fortran_eos/tests/run_matching_replay.sh \
  "$PARENT_MATCHER" "$CANDIDATE_MATCHER" "$FROZEN_CATALOG_DIR" "$NEW_SCRATCH_DIR" \
  BPS_MTW1964_MS1966 NV_DME_printed_fixed_NZ NV_DME_D03x100_fixed_NZ \
  BPS_MTW1964_MS1966_NV1973_Hermite
```

The frontend launches the existing native matcher and compares its serialized results. No stellar solutions or matching physics are recomputed in Python. The unchanged historical per-input digest policy still applies to the original campaign generator.
