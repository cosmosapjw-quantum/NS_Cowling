# Native quadrupolar Cowling regression

The existing C16 Fortran detector recomputed four roots on two archived, single-interface equilibrium stars from the C12 original-source fixtures. The M2 Maxwell model uses rapid conversion. The M3 constant-sound-speed model uses slow conversion. Both runs use `l=2`. The `f` and `p1` labels follow the brackets and branch identification from the prior C16 continuation. This replay evaluates signed determinants and roots on the physical profiles. Its separate node-count and overlap unit checks use synthetic profiles.

| Star | Root | Archived-control frequency (Hz) | Tighter-control frequency (Hz) | Change (Hz) |
|---|---|---:|---:|---:|
| M2 | f | 1842.311657463 | 1842.311651717 | −0.000005746 |
| M2 | p1 | 5143.719615379 | 5143.719686989 | +0.000071610 |
| M3 | f | 1564.244880395 | 1564.244862734 | −0.000017661 |
| M3 | p1 | 3736.210685077 | 3736.210549108 | −0.000135969 |

The native detector requires opposite determinant signs about each root, successful Brent termination, an absolute determinant residual at most `1e-5`, and agreement with each driver's reference frequency within `0.02 Hz`. The scan and root solves reproduced the expected source-knot counts of 39 for M2 and 40 for M3. Three repeats per configuration produced exactly equal double-precision roots within the executable. The checked and optimized builds emitted identical scientific text for both control settings. At the archived tolerances the full scientific output for each star matched the corresponding archived `FORTRAN_OPT_M*.log` byte for byte.

The tighter run sets relative tolerance to `2e-10` and absolute tolerance to `2e-13`, from `1e-9` and `1e-12`, and halves the maximum integration step. The stored TOV background, EOS interpolation, interface treatment and Brent tolerance stay fixed. The paired results measure sensitivity to these three ODE controls. No background mesh convergence is inferred. Multi-interface reconstructed stars were outside this calculation.

## Reproduction

From the repository root with GNU Fortran, BLAS development libraries, GNU coreutils, Bash and the two compressed C12 inputs available:

```sh
bash data/eos_native_upgrade_20261124/cowling/replay.sh > modes-rerun.csv
diff -u data/eos_native_upgrade_20261124/cowling/modes.csv modes-rerun.csv
```

An optional first argument accepts a directory containing the original fixture filenames. The shell frontend verifies the input hashes, changes the designated integration settings, invokes the existing Fortran build and detector, and formats the emitted frequencies. All root and determinant calculations use Fortran. The measured toolchain was GNU Fortran 14.2.0 on Linux x86-64. The builder uses `-O0 -g -fcheck=all` or `-O2` for project code, and `-O2 -std=legacy` for frozen Netlib routines in both builds. Cross-compiler reproduction should compare the numerical frequencies within a stated tolerance.

Baseline commit `cab027e9b5ccbc42d6509531f7a6201f2c13d184`. Inputs live under `data/cowling_fortran_source_exact_root_c12_20261121/`. The compressed SHA-256 values are pinned in `replay.sh`.

| Star | Git blob of compressed fixture | SHA-256 of decompressed input |
|---|---|---|
| M2 | `fb2974d4a367d80a5803efddf343f75b9baf9b28` | `eee499eb50de4300f9987e05bf46a4fee3ddcbfb2e9b9d99deb9ac8a423c35e5` |
| M3 | `40e3b48f9c291d8bb49940463c3bb682d80025b0` | `f6a7abb566ae68b8ac5f838d410771823634a98024007ca02264bd6e44a0d3b1` |

The test reads each `p1` reference from the C12 input trailer and each `f` reference from the full-precision anchor recorded in `scripts/run_c16_source_density.py`. These two references define the `difference_from_driver_reference_Hz` column of `modes.csv`. Original C16 scientific outputs are in `data/cowling_c16_source_continuation_20261121/FORTRAN_OPT_M{2,3}.log`. The earlier branch and stellar background evidence is in `docs/COWLING_C16_SOURCE_F_P1_LOCAL_CONTINUATION_20261121.md`.
