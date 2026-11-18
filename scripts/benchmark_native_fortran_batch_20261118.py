"""Measure no-Python-numerics CLI batch TOV/Love speed on same native EOS."""
from pathlib import Path
import json,subprocess,time,statistics,os,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/fortran_dvode_perf_20261118"
SCRATCH=Path("/home/cosmosapjw/ns_cowling_wt/scratch/fortran_dvode_perf_20261118")
BINARY=SCRATCH/"dvode/build/compose_eos_fortran"
SOURCE=SCRATCH/"baseline/PK1_native_4097.csv"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(args):
    env=os.environ.copy()
    env.update(NS_COWLING_TOV_INTEGRATOR="dvode_adams",
               NS_COWLING_TOV_RTOL="2e-9",
               OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
    t=time.perf_counter_ns()
    obj=subprocess.run([str(BINARY),"tov-love",*map(str,args)],
            capture_output=True,text=True,env=env,timeout=90)
    seconds=(time.perf_counter_ns()-t)/1e9
    if obj.returncode!=0:raise RuntimeError((obj.stdout+obj.stderr)[-600:])
    result=[json.loads(t) for t in obj.stdout.splitlines() if t.startswith("{")]
    return seconds,result
def benchmark():
    from compose_eos import read_compose
    original=ROOT/"data/native_dense_interpolation_20261118/rmf/PK1_CANONICAL/table/compose"
    raw=read_compose(original)
    densities=np.linspace(.20,.55,32)
    indices=np.searchsorted(raw.nb,densities)
    pressures=[float(raw.p[i]) for i in indices]
    cfg=SCRATCH/"batch/PK1_centers_32.csv"
    cfg.parent.mkdir(parents=True,exist_ok=True)
    np.savetxt(cfg,pressures,fmt="%.17e")
    stats=[]
    serial_reference=None
    numeric_fields=("R_km","M_msun","k2","Lambda","y_surface","n_steps")
    for batch in (False,True):
        times=[]
        mode="batch" if batch else "32 individual native process invocations"
        for k in range(3):
            if batch:
                elapsed,results=run(["batch-smooth",SOURCE,cfg])
            else:
                start=time.perf_counter_ns()
                results=[]
                for p in pressures:
                    _,z=run(["smooth",SOURCE,p])
                    results.extend(z)
                elapsed=(time.perf_counter_ns()-start)/1e9
            if len(results)!=len(pressures):raise AssertionError("lost batch numerical solutions")
            if not all(x["converged"] for x in results):raise AssertionError("native batch partial failure")
            if not batch and k==2:
                serial_reference=results
            if batch:
                if serial_reference is None:raise AssertionError("native serial reference not recorded")
                for row_a,row_b in zip(results,serial_reference):
                    if any(row_a[f]!=row_b[f] for f in numeric_fields):
                        raise AssertionError("native Fortran 32-star batch differs numerically from serial mode")
            times.append(elapsed)
        stats.append(dict(mode=mode,seconds=times,median_seconds=statistics.median(times),
            result_checksums=[sha(SOURCE),sha(cfg)],valid_solutions=len(results)))
        print(mode,stats[-1]["median_seconds"],flush=True)
    result=dict(status="BENCHMARK_COMPLETED",
      model="PK1_CANONICAL",source_state_count=4097,central_pressures=32,
      physics_scope="uniform PK1 homogeneous source truncated at n=.015 fm^-3; not a complete NS",
      fixed_numerical_controls={"integrator":"dvode_adams","rtol":2e-9},
      python_role="input specification and subprocess/json wall-clock only",
      serial_vs_batch_all_32_states_exact=True,
      numeric_fields_compared=list(numeric_fields),
      rows=stats,
      measured_batch_speedup=stats[0]["median_seconds"]/stats[1]["median_seconds"])
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"BATCH_BENCHMARK.json").write_text(json.dumps(result,indent=2)+"\n")
    print("BATCH_SPEEDUP",result["measured_batch_speedup"])
if __name__=="__main__":
    benchmark()
