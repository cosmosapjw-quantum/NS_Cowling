# Native generation entry points

The production entry point is `compose-eos generate`. Its Python code validates requests, reads immutable input coefficients, starts the native executable, records provenance, and serializes the native receipt. Physical evaluation, equilibrium, geometry choice, derivatives, and CompOSE readback execute in Fortran. There is no automatic Python fallback.

`python scripts/generate_core_eos.py OUTPUT` is a compatibility batch frontend for all five core models. Each model has its native `generated/` directory plus byte-identical `eos.*` copies at the historical batch path. `RECEIPT.json` serializes native `states.csv` and `first_law.csv`; it does not calculate nuclear states or finite differences in Python. `DISPATCH.json` preserves the exact request, native source hashes, exit status and output hashes. Failure retains its native execution evidence and marks the batch incomplete.

The narrow-density CLI regression uses APR [1.20,1.21] fm^-3, inside the declared native [0.085,1.34] support. Its original [1.95,1.96] interval exercised an unadvertised extrapolation of the historical Python oracle. The regression still checks that both native five-point derivative stencils stay inside the requested interval.

## Historical research oracles

`scripts/reoptimize_spg_surfaces.py` and `scripts/build_gppva_tf.py` retain the previously tested Python research algorithms and importable audit helpers. Direct execution requires the explicit flag `--historical-python-oracle`. These commands are not the production backend and their receipts must not be relabeled native physical generation: historical “native” wording referred to Fortran file I/O/readback.

Old reproduction documents refer to their recorded source commit and checksums. To reproduce the old campaign byte for byte, use that exact commit; the new flag identifies an explicit research-oracle run in the current checkout. Preserving the old helpers prevents old profile reconstruction and negative-case tests from being silently rewritten around the new solver.

The independent Python functional implementations also remain test oracles. Numerical agreement with them supplements conservation, stationarity, thermodynamic derivative and primary-source checks. Neither Python-oracle agreement nor CompOSE agreement alone establishes physical correctness or original-author equivalence.
