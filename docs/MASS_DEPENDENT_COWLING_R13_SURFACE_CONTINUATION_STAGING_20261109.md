# R13 surface-continuation staging receipt

Date: 2026-11-09 (KST)

## Status

**STAGED_OFF_BRANCH: executable source package is in Dropbox, not yet integrated into this GitHub branch.**

Target repository: `cosmosapjw-quantum/NS_Cowling`  
Target branch at staging: `research/mass-dependent-cowling-v1`  
Branch HEAD at staging: `074974784b5c31afe6b24703dbc791e0bfc43307`

The local package is backed up at:

`/TOV_Cowling_Solver_Backups/ns-cowling-thermodynamic-surface-continuation-R13.zip`

Dropbox file ID: `id:BSpOijBcT10AAAAAAELDrA`  
Size: 19,216 bytes  
ZIP SHA-256: `5adc883b924b181c945f2cab42037ff9c656956e3defbb8f492ddbbc7b0bd879`

The ZIP contains:
- `python_frontend/compose_cowling/surface_continuation.py`
- `python_frontend/tests/test_surface_continuation.py`
- `scripts/run_surface_continuation.py`
- `.github/workflows/surface-continuation.yml`
- `docs/plans/2026-11-09-thermodynamic-surface-continuation.md`
- local test output, diagnostic JSON and per-file SHA-256 manifest

## Local verification (not GitHub CI)

- Focused tests: `8 passed in 0.09 s`.
- `compileall` for the module and runner passed.
- The five-family, 512-point thermodynamic diagnostic gives positive subluminal `cs2` and maximum first-law log-derivative residual approximately `2.22e-8`.

The EOS cutoff is (P_{cut}=9.000000079422\times10^{-12}\;\mathrm{MeV\,fm^{-3}}), (\epsilon_{cut}=4.354134734530827\times10^{-7}\;\mathrm{MeV\,fm^{-3}}). The module uses the exact right-sided raw Hermite (P(\mu_b)) tangent, distinguishes (\gamma_s) from (\Gamma_{eq}), integrates (d\ln n_b/d\ln\epsilon=\epsilon/(\epsilon+P)), and reconstructs (\mu_b=(\epsilon+P)/n_b).

## Scientific limits

This is a cold barotropic mathematical continuation, not a sourced or composition-resolved outer-crust EOS. The table writer emits positive tail rows with (x_{min}=10^{-8}) and (10^{-12}), because the native reader rejects literal zero primary fields. The included 42-case native TOV runner and dedicated workflow have **not** yet executed. No TOV M/R sensitivity is known, no continuation family is selected, and no radial (\omega_0^2), stability classification or new f/p1 result is claimed.

## Integration blocker and next action

The current connector write path could not safely transfer the larger executable files as Git blobs. The attempted malformed payloads were never linked to any branch or commit. This record is a documentation-only GitHub push; it does not pretend the executable source has been integrated.

Next integration must add the five package files to the target branch, verify every file digest against the ZIP manifest, run the dedicated surface-continuation Actions workflow, inspect all 42 native v5 TOV cases, and only then update the scientific status. Do not relax frozen-v5 baseline M/R checks.

PR #3 remains open/draft/unmerged and was not modified.
