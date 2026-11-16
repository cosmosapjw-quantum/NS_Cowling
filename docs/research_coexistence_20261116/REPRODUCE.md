# Reproduction and identity

Base commit: `a0b158f3dc4af9825192bf6a185a6331619868a5`. Use a clean isolated checkout. The research calculations are conditional finite samples, not production EOS promotion. The TF and SPG reports specify their own actual commands, inputs, failure records and limitations.

```sh
python3 data/research_coexistence_20261116/provenance_check.py > /tmp/coexistence_provenance.json
sha256sum -c SHA256SUMS_v4.txt
bash fortran_eos/build.sh /tmp/ns_coexistence_build
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=/tmp/ns_coexistence_build/compose_eos_fortran python3 -m pytest -q python_frontend/tests pipeline/tests
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=eos_generator:eos_generator/compose_eos_pkg python3 -m pytest -q eos_generator/tests
```

The baseline checksum verification returned 8,114 OK lines; the complete console log SHA256 was `64d6e7bbe9d10a91566fd0fc4421ae5a5645bb02d6585639edea5394188f6ed8`. The build and tests are recorded in `data/research_coexistence_20261116/verification/`: frontend/pipeline 1078 passed, 9 skipped; generator 142 passed and one environment Matplotlib warning. This work changes no production code.

The archived source-grid manifest contains 42 identities. Three files differ from the current base: `fortran_eos/app/compose_eos_generate.f90`, `fortran_eos/build.sh`, and `python_frontend/compose_eos/native_generation.py`. `PROVENANCE.json` reports these byte mismatches; all three historical bytes are recovered and match their recorded hashes at `a49aaf9f82ddfda484e0f965aefdfe73aa7918f2`. The other 39 match. The historical N128 recovery receipt separately binds 15 linked sources; these same three differ now, while its other 12 match. Review of the base diff identifies additions for diagnostic generation, linking and diagnostic hashes; the mismatch itself is a byte-identity fact, not a scientific-integrity verdict. Do not attribute the historical binaries to the current complete source tree. Preserve the original failed dispatches and native readback records.

The SPG wrapper was additionally rerun from this checkout with `python3 fortran_eos/tests/spg_residual_run.py /tmp/spg_coexistence_independent_run`: exit 0, 760 inputs, 3,800 data rows, and stdout SHA256 `b9ccf639b61150f238cc8b90d8034aa0d6853b8798c2a0598b986d28921f8160`, equal to the preserved guarded output. The independent reviewer also inspected its failed-input path. This checks reproduction of the sampled derivative calculation only.
