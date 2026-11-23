# First-law pressure matching benchmark

For each candidate endpoint pair the first-law bridge still solves for the same intermediate plateau pressure by bisection. On each of its three cubic pieces, the pressure-work integral uses the same GL16 nodes, positive weights and order of summation. The Hermite basis and baryon density at a full piece's nodes depend only on its geometry. The Fortran kernel now precomputes these once per candidate. Once the plateau is selected, every output query begins with the cumulative work at its piece's left knot. The causality/extrema certificate and endpoint search order are retained. Independent output nodes use a static OpenMP loop for arrays with at least 256 queries. Both the fixed 401-node certification sample and the requested output grid count separately for this condition.

The native harness times the existing production routine including search, certificates, final table construction, and allocation, after source file reading and one warmup. The 80- and 4096-node sets have four interleaved blocks of eleven calls per configuration. Their reported values are medians in seconds. The 65536-node stress case has five processes of four timed calls each. The median of its per-process peak RSS values in KiB was 17032 (parent), 17220 (candidate one thread), and 17260 (candidate two threads). Process RSS therefore did not improve on this case. The constructor eliminated three transient arrays of length equal to the query grid in the Hermite-pressure and CSS paths. The covered parent/candidate output digests match, including every Q column, native report, and phase.

| Source pair | Bridge nodes | Parent, 1 thread (ms) | Candidate, 1 thread (ms) | Candidate, 2 threads (ms) |
|---|---:|---:|---:|---:|
| BPS + NV / D1MS | 80 | 2.067 | 1.684 | 1.713 |
| BPS + NV / D1MS | 4096 | 4.996 | 2.386 | 2.372 |
| BPS / DDMEX | 80 | 74.605 | 47.749 | 47.533 |
| BPS / DDMEX | 4096 | 79.007 | 50.795 | 50.967 |
| BPS + NV / D1MS | 65536 | 48.570 | 12.970 | 9.022 |

The shared machine and processor affinity limit the inference about other scales and hardware. High query count benefits most from parallel evaluation. The finite table checks and pressure-work unit tests use distinct controls. A second Fortran integrator using a composite GL4 rule on 32 subintervals per pressure piece agrees with GL16 within 1.31e-15 in the knot and log-density test set. Physical units are MeV and fm, with the pressure-work integral giving MeV per baryon after multiplication by the log-density span.

For a reproduction, use the two source-selected internal directories with the hashes in metadata.json. The public core inputs are found via the core source catalog and both generated reconstructed crusts can be converted to internal CSV through the existing reconstruction runner. Build each revision with the same flags and native harness. The final argument records a native binary diagnostic for byte comparison on the same compiler and platform. Timings exclude file reading and writing.

~~~sh
src="$repo/fortran_eos/src"
gfortran -std=f2018 -O2 -g -fopenmp -fcheck=all -fimplicit-none -ffree-line-length-none \
  "$src/compose_eos_types.f90" "$src/compose_eos_io.f90" "$src/compose_eos_matching.f90" \
  "$repo/fortran_eos/tests/bench_firstlaw_matching.f90" -o bench_firstlaw
OMP_NUM_THREADS=2 OMP_PROC_BIND=true OMP_PLACES=cores OMP_DYNAMIC=false \
  taskset -c 4,5 ./bench_firstlaw "$crust_internal" "$core_internal" 11 4096 case.bin
~~~
