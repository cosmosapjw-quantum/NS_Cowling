# Native NV profile readback

The previous verifier accepted a generated uniform-beta CSV after replacing only the first electron fraction with 0.49. The revised verifier reads the CompOSE slice and evaluates charge neutrality, electron fraction, the chemical potentials, beta affinity, Euler identity, and sound-speed and curvature conditions in Fortran. For the two constrained profile families it also compares the reported energy per baryon and pressure with both native mesh outputs and checks the finite-dilation and composition derivatives. The Python entry point parses named CSV columns, invokes the native verifier and renders its records as JSON.

The normal 128/256 radial calculation and eight deliberate corruptions or misuses appear in the control log. The legacy 128/256 run also passes; older directories lack mesh metadata, so the JSON marks the assumed grids. The generator now writes profile_run.nml for new results. Bounds on derivatives and secants apply only to the sampled values and the indicated composition.

~~~sh
bash fortran_eos/build.sh --ensure
python3 scripts/verify_nv_reconstructed_profiles.py "$COARSE_DIR" "$FINE_DIR"
bash fortran_eos/tests/run_nv_profile_audit_controls.sh "$COARSE_DIR" "$FINE_DIR" build/compose_eos_fortran/nv_profile_verify
~~~
