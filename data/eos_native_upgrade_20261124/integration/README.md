# Full native build integration

The production executable was built with GNU Fortran 14.2.0. Its SHA-256 is `a92aa63e661738e09dcaf8dcf60dcadc4e98caf71d087ada790ca19bd162b1c5` and the full build's source signature is `c4a80dcb08d7e3791dd900ec5e56d780d85d0850a4111da81b451d7541c28858`. All checks below completed with exit status zero on one CPU.

| Check | Outcome |
|---|---|
| Seven legacy native matching cases on the production executable | Six successful tables with all 30 CompOSE files byte-identical, one identical failure. `legacy_comparison.json` identifies the pinned reference and catalog. |
| Production matching module | Both the physical slice conventions and exhaustive finite window search tests passed. |
| Native internal matching checks | Four archived drivers passed using the repository harness that exposes private numerical routines only in its test-local module. |
| TOV–Love | Smooth and Maxwell tests, the independent Hinderer normalisation and tidal RHS tests, and layered-interface and mass-target tests passed with production objects. The relative difference in the tidal RHS check was 1.12e-13. |
| PK1 production CLI | Four batch stars matched the separate single-star results exactly in radius, mass, Love number, tidal deformability and surface y. A negative central pressure produced a nonzero exit. |
| Native Cowling oscillations | Three standalone Fortran drivers linked to the normal build's ODEPACK objects. Weak-gravity uniform-star roots for Adams/BDF were 67.473154193/67.473154194 Hz. Four sharp-interface cases with slow/rapid conversion were each solved with Adams and BDF. The largest absolute difference from the frozen frequency was 0.002176 Hz against a 0.005 Hz gate. False sign flips and invalid-input controls passed. |

`native_controls.txt` contains the native output. The existing Cowling test drivers emitted two nonfatal GNU continuation-format warnings in diagnostic output statements. The original five-file PK1 CompOSE fixture and compressed one-interface Cowling fixture were recovered with Git from `cab027e9b5ccbc42d6509531f7a6201f2c13d184`. The focused run did not skip any fixture. The broader matching replay uses the separately identified matcher builds in `../matching_replay/`.

From the workspace containing the checkout and its baseline run, with the ordinary build and the populated legacy catalog available, the executed command was

```bash
export PATH="$PWD/work/toolchain/bin:$PATH"
export PYTHONPATH="$PWD/work/bps-testdeps"
export NATIVE_LIBS='-l:liblapack.so.3 -l:libblas.so.3'
taskset -c 8 NS_Cowling/data/eos_native_upgrade_20261124/integration/replay.sh \
  work/native_upgrade/full work/matching/catalogs work/matching/runs work/integration_regression/final3
```

Use a fresh fourth directory for each run. Python schedules jobs, decodes the archived CompOSE input for the existing PK1 CLI test, and compares native receipts. The Cowling, matching and TOV–Love solvers run in Fortran.
