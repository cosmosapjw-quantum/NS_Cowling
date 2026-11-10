# Sharp-Maxwell Cowling mode-identity research slice

Status: PASS_SYNTHETIC_HELPER_TESTS_NOT_A_PHYSICAL_MODE_CLASSIFICATION.

This package is independent of the blocked C7X CompOSE crust composition admission. It adds exact slow/rapid interface maps, a positive profile norm and localization diagnostics, one-sided interface polarization, phase-aligned overlap, and conservative spectral-branch assignment.

It does **not** assign f-, p1-, or discontinuity-g-mode labels. Existing phase-local W/V node counts remain the responsibility of the project Python reference implementation. The helper is not integrated into the production Fortran API and has not been run on the four native fixture-root profiles.

Run from the repository root:

```bash
python -m pytest -q research/maxwell_mode_identity/tests
```

Files:
- `interface_transfer.py`: exact two-phase interface maps and residual checks.
- `diagnostics.py`: profile norm, trapping fractions, interface diagnostics, normalized overlap and branch matching.
- `profile_schema.json`: JSON contract for real-valued finite-frequency profile records.
- `MODE_IDENTITY_AND_INTERFACE_MAP_20261111.md`: derivation, literature synthesis and fixture gates.
- `MODE_IDENTITY_STATUS.json`: machine-readable status and provenance.
- `tests/`: deterministic synthetic tests.

The overlap is a chosen continuation similarity, not a physical observable or a proven canonical inner product between distinct stellar backgrounds. Slow and rapid spectra must be tracked separately. Branch identity and physical modal character remain separate data fields.
