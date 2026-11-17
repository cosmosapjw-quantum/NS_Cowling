# Final band-operator recheck of frozen crust seams

Four previously optimized N384 spherical seam profiles (DD2, TW, FSU2R,
FSU2H) were re-evaluated with spatial_tf.py SHA
2d1f1edc72b170b01194f2705e62905383a2d4292bd84cc58067bc2cc1067c02
and spatial_tf_equilibrium.py SHA
8be3c6756f1b04a38cb59fdeeacb5c0aed5c98d8337c8ed14c64e7091501bdca.

No nuclear profile or cell-radius minimization was performed. Meson/electron
fields were recomputed at the stored nuclear profiles and radius. This
verifies current-operator support of the stored finite-grid stationary
candidates; it does not establish mesh convergence, a new interface root,
or a global minimum.

Maximum changes from the original seam states:

| Quantity | Maximum absolute change |
|---|---:|
| Energy density (MeV fm^-3) | 4.440892098500626e-16 |
| Pressure (MeV fm^-3) | 5.329070518200751e-15 |
| Baryon chemical potential (MeV) | 2.3874235921539366e-12 |

Maximum chemical KKT residual was 5.7412684327573515e-5 MeV, radial residual
7.68702944882664e-7 MeV, baryon constraint error 1.3760693973186022e-15 fm^-3,
and charge constraint error 4.685150989500349e-20 fm^-3. The original physical
acceptance thresholds were retained. All four checks passed.

From the repository root:

```sh
PYTHONPATH=python_frontend OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 data/research_finish_20261116/crust_band_recheck/verify_band_profiles.py > data/research_finish_20261116/crust_band_recheck/run.log 2>&1
```

The executed verifier writes exclusively into this directory. Exact executed
source bytes, input identities, runtime versions, and old/current operator
identities are retained in BAND_PROFILE_RECHECK.json and SOURCE_SNAPSHOT.json.gz.
The earlier crust_assembly/final_operator_profile_verification.json and all
earlier sealed crust evidence were retained unchanged.
