# C6 native execution evidence (C6_NATIVE_EXECUTION_20261110T093714Z)

Status: **C6_PARTIAL**. Native build, committed pytest (6/6), native probe of SPG_M1 and GPPVA_FSU2, and round-trip all executed and passed. Native C3/C4 composition closure FAILS on Iphase=7 rows of both sources; unresolved, needs review. Raw sources and 234 legacy files byte-identical before/after. 65 pairs unchanged.

Head: 43eea8bd9bee7ebb222e3c66babf5078b398d7ae. Executable sha256: 8c9153a93dac329b8bec5fafa85d7318f9a82dd5767d5b4cceea21f79a912d0f.
See execution_report.json, roundtrip_validation.json, probe_*.json. Verify with `sha256sum -c SHA256SUMS` inside this directory.
Note: native exe built outside the repo; path in environment.json. probe_combined.json is the single committed-runner output; probe_<ID>.json are per-source slices of it.
