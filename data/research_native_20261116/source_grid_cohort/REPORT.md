# Native original-density TF cohort

Status: FINITE_COHORT_NATIVE_VERIFIED; unchanged source bytes: True.

The actual public CLI regenerated all **590 inner-crust density locations** from the eight pinned GPPVA source density grids, at N16, using the native finite-family TF solver. Only `eos.nb` density locations were taken as inputs. Source thermodynamic values were neither loaded by the physics backend nor used as equality gates.

All eight requests and all 24 public native readbacks passed: 590 main rows and 1902 rows including separately regenerated outer/BSk source-fit segments. There were zero final held nodes. The retained attempt count is 3641, including 87 internal EVALUATION_FAILED trials (FSU2R:30, FSU2H:49, TM1e:8). Those failed trials were retained; bounded subsequent attempts resolved every requested family and yielded zero final held nodes.

| Model | Native inner rows | Attempts | Maximum selected KKT (MeV) |
|---|---:|---:|---:|
| GPPVA_DD2 | 66 | 396 | 0.000234243 |
| GPPVA_DDME2 | 71 | 426 | 0.000731485 |
| GPPVA_TW | 72 | 432 | 0.000533875 |
| GPPVA_FSU2 | 52 | 312 | 0.000276678 |
| GPPVA_FSU2R | 79 | 512 | 0.000440929 |
| GPPVA_FSU2H | 83 | 553 | 0.000710869 |
| GPPVA_NL3wrL55 | 81 | 486 | 0.000987124 |
| GPPVA_TM1e | 86 | 524 | 0.000327155 |

This is a regeneration and finite-grid numerical-support result. N16 is coarse. The selected geometries and their transition locations do not establish mesh convergence or a global minimum. In particular, low-density cells can select different geometry at this coarse mesh. The new results must not supersede the preserved N192/N384 research campaign as evidence of continuum physical ordering. The native SUPPORT metadata governs same-phase interpolation and keeps geometry crossings separate.

The lower/source attachment and original-source fidelity remain separate from successful generation. The regenerated outer/BSk seam has explicit pressure, chemical-potential and density mismatch; it is exported as distinct segments. No artificial energy shift, unrequested common-tangent gate or cross-gap interpolation was added.

Reproduce from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 data/research_native_20261116/source_grid_cohort/RUN.py
```

Each retry creates a new attempt directory. All requests, native profiles, attempted states, selected states, native tables, readback logs, checksums and aggregate receipts remain in this directory. The exact requested input grids and density-file identities are embedded in SUMMARY.json.
