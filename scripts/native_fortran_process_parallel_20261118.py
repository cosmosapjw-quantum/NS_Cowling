"""Independent-process parallel stress; numerical work exclusively in Fortran.

Never call VODE from concurrent threads in the same executable: DVODE and its
LINPACK workspace use static/common storage. ThreadPool only orchestrates
independent *processes*, each with isolated Fortran solver state.
"""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import hashlib, json, os, subprocess, time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/fortran_dvode_perf_20261118"
SCRATCH = Path("/home/cosmosapjw/ns_cowling_wt/scratch/fortran_dvode_perf_20261118")
BINARY = SCRATCH / "dvode/build/compose_eos_fortran"
SOURCE = SCRATCH / "baseline/PK1_native_4097.csv"
PRESSURES = (23.0, 35.0, 43.0, 55.0, 73.52922470749803, 90.0)


def native(pressure):
    env = os.environ.copy()
    env.update(
        NS_COWLING_TOV_INTEGRATOR="dvode_adams",
        NS_COWLING_TOV_RTOL="2e-9",
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    p = subprocess.run(
        [str(BINARY), "tov-love", "smooth", str(SOURCE),
         f"{pressure:.17e}"],
        env=env, capture_output=True, text=True, timeout=45,
    )
    if p.returncode != 0:
        raise RuntimeError(f"Fortran native parallel source failed: {p.stdout} {p.stderr}")
    x = json.loads(p.stdout)
    if x.get("status") != "ok" or not x.get("converged"):
        raise ValueError("not an admitted native TOV/Love solution")
    return x


def main():
    t0 = time.perf_counter_ns()
    serial = {p: native(p) for p in PRESSURES}
    serial_secs = (time.perf_counter_ns() - t0)*1e-9

    t1 = time.perf_counter_ns()
    with ThreadPoolExecutor(max_workers=6) as executor:
        futures = {executor.submit(native, p): p for p in PRESSURES}
        parallel = {futures[f]: f.result() for f in as_completed(futures)}
    parallel_secs = (time.perf_counter_ns() - t1)*1e-9

    fields = ("R_km", "M_msun", "k2", "Lambda", "y_surface", "n_steps")
    for p in PRESSURES:
        for field in fields:
            if serial[p][field] != parallel[p][field]:
                raise AssertionError(
                    f"fortran process-parallel result differs from serial: {p} {field}"
                )

    rec = dict(
        schema="ns_cowling.native_process_parallel_isolation.v1",
        status="PASS",
        source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        Fortran_executable_sha256=hashlib.sha256(BINARY.read_bytes()).hexdigest(),
        number_distinct_stars=len(PRESSURES),
        max_parallel_native_processes=6,
        num_threads_per_process=1,
        parent_Python_role="subprocess orchestration and equality comparison only",
        Fortran_state_isolation="each worker is one separate process, never DVODE in OpenMP threads",
        serial_seconds=serial_secs, parallel_seconds=parallel_secs,
        speedup_serial_over_parallel=serial_secs/parallel_secs,
        results=[dict(central_pressure=p, **serial[p]) for p in PRESSURES],
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "PARALLEL_STRESS.json").write_text(json.dumps(rec, indent=2)+"\n")
    print(
        "PARALLEL_FORTRAN_PROCESS_PASS",
        f"{len(PRESSURES)} stars",
        f"serial {serial_secs:.5f}s",
        f"parallel {parallel_secs:.5f}s",
        f"speedup {serial_secs/parallel_secs:.3f}x",
    )


if __name__ == "__main__":
    main()
