# APR support selection and full-table limitation

The parent requested a read-only decision after the regenerated APR core was rejected by the new explicit SUPPORT guard. No production code was changed in this follow-up.

The actual regenerated table at `data/research_matching_20261117/integration/real_pairs_v2/core/APR_APR/generated` has 131 rows:

| Opaque provider phase ID | Original zero-based rows | Density interval fm^-3 | Rows |
|---|---|---|---|
| 6 | 0–64 | 0.085–0.20439495489345874 | 65 |
| 8 | 65–75 | 0.20615528128087748–0.23676398610352045 | 11 |
| 7 | 76–130 | 0.23738424893920132–0.5 | 55 |

The generator declares cross-phase interpolation UNAPPROVED. Its genuine APR phase structure must not be erased to enter the matcher. In current `compose_eos_matching.f90`, matched output initializes Iphase=0 and disables the composition sidecar; retained core rows i2:n copy thermodynamic columns but never copy phase labels. Therefore a guard inspecting only the local endpoint stencil would permit later source phase boundaries to disappear from the output. The current CSS stellar adapter also replaces all rows with two crust/CSS-core labels; it does not support retaining additional core phase segmentation. This is an **implementation/input-support limitation**, not evidence that APR matching or its physical transitions are impossible.

A future relaxation would require preserving retained phase-boundary metadata through native output and CompOSE serialization, independently deciding the allowed interpolants on each segment, and a downstream representation that preserves all retained source boundaries. Label changes cannot by themselves determine whether an interface is a jump, a continuous phase boundary or a geometric change. Only after those obligations are implemented is a precise consumed-window/stencil check sufficient to replace whole-input branch preselection.

## Executed contiguous-branch route

`data/research_matching_20261117/pipeline_matching/apr_branch_route.py` selects exact rows 0–64 of the original APR table, including every per-row composition/presence array; the phase tags remain 6 and `source_dir` remains the original source. No source file or source energy is modified. It matches this bounded low-phase table to the repaired native SPG_M1 crust through all five public methods, then writes and rereads each result with native CompOSE I/O.

Final native executable SHA256: `920a4512e56d028f3f7e3e03680d5e32b915d297b3936963c94467b9ebde716f`.

| Method | Output rows | Selected-core endpoint index | Endpoint n_B fm^-3 | API construction / readback |
|---|---:|---:|---:|---|
| css | 222 | 13 | 0.10175984976746866 | verified / PASS |
| cubic_F | 223 | 11 | 0.09898108760515353 | verified / PASS |
| cubic_P | 223 | 11 | 0.09898108760515353 | verified / PASS |
| quintic_F | 223 | 11 | 0.09898108760515353 | verified / PASS |
| quintic_P | 223 | 11 | 0.09898108760515353 | verified / PASS |

All five retained-core F values agree with the selected source at max absolute error 0.0 MeV, retained density and Q1 arrays have exact numerical identity, and reported dm_crust=0. Source-file hashes are unchanged. Native consistency checks pass, so none is marked diagnostic-only. These results cover the finite low-phase interval ending at 0.20439495489345874 fm^-3; they do not produce an assembled full multiphase APR stellar EOS.

Evidence: `data/research_matching_20261117/pipeline_matching/apr_low_branch/RECEIPT.json`, sibling method directories containing actual CompOSE outputs, and `apr_branch_route.log`. The command exited 0 with one numerical thread and no rebuild:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONPATH=python_frontend:. \
COMPOSE_EOS_FORTRAN_EXE=$PWD/build/final_matching/compose_eos_fortran \
python data/research_matching_20261117/pipeline_matching/apr_branch_route.py
```

Public-API support errors should be classified as input/support or explicit implementation scope by the integration harness. The pipeline's `source_support` stage already maps to input_contract; no source-support failure should count as a physical no-transition result.
