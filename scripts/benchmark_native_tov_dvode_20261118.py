"""Python frontend ONLY: launch executable, compare Fortran numerical outputs.

All EOS interpolation, TOV and Love computations run inside the Fortran binary.
Never include the legacy uncorrected k2 or Lambda as a physics reference.
"""
from __future__ import annotations
from pathlib import Path
import hashlib,subprocess,os,json,statistics,time,random
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"data/fortran_dvode_perf_20261118"
SCRATCH=Path("/home/cosmosapjw/ns_cowling_wt/scratch/fortran_dvode_perf_20261118")
INPUT=SCRATCH/"baseline"
SOLVER=SCRATCH/"dvode"/"build"/"compose_eos_fortran"
OLD=INPUT/"build"/"compose_eos_fortran"
CORRECTED_REFERENCE=SCRATCH/"corrected"/"build"/"compose_eos_fortran"
PC=73.52922470749803
PHYSICAL_OUTPUTS=("M_msun","R_km","compactness","k2","Lambda","y_surface")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def call(bin_path,table,mode,rtol="2e-9"):
    env=os.environ.copy()
    env["OMP_NUM_THREADS"]="1"
    env["OPENBLAS_NUM_THREADS"]="1"
    env["NS_COWLING_TOV_INTEGRATOR"]=mode
    env["NS_COWLING_TOV_RTOL"]=rtol
    command=[str(bin_path),"tov-love","smooth",str(table),str(PC)]
    before=time.perf_counter_ns()
    p=subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True,timeout=45)
    elapsed=(time.perf_counter_ns()-before)/1e9
    if p.returncode!=0:raise RuntimeError(f"native exit {p.returncode}: "+p.stdout[-700:]+p.stderr[-700:])
    obj=json.loads(p.stdout)
    if obj["status"]!="ok" or not obj["converged"]:
        raise ValueError("native science solve rejected: "+str(obj))
    return dict(seconds=elapsed,reported=obj,command=command,
        environment={k:env[k] for k in ("NS_COWLING_TOV_INTEGRATOR","NS_COWLING_TOV_RTOL","OMP_NUM_THREADS")})
def benchmark():
    if not all(p.exists() for p in (SOLVER,OLD,CORRECTED_REFERENCE)):raise FileNotFoundError("native Fortran binaries required")
    OUTPUT.mkdir(parents=True,exist_ok=True)
    cases=[]
    for n in (257,1025,4097):
        table=INPUT/f"PK1_native_{n}.csv"
        target=[("corrected_dvode_adams",SOLVER,"dvode_adams"),
                ("corrected_dvode_bdf",SOLVER,"dvode_bdf"),
                ("corrected_rk45",SOLVER,"rk45"),
                ("legacy_binary",OLD,"rk45")]
        random.Random(n).shuffle(target)
        for label,binary,mode in target:
            trial=[]
            # Warm first launch avoids claiming favorable one-off cache timing.
            warm=call(binary,table,mode)
            for repeat in range(7):
                row=call(binary,table,mode)
                row["repeat"]=repeat
                trial.append(row)
            median=statistics.median(x["seconds"] for x in trial)
            case=dict(model="PK1_CANONICAL",source_density_range_fm3=[.015,1.0],
              surface_status="UNIFORM_CORE_ONLY_TRUNCATED_SOURCE; not a complete neutron star",
              source_table_rows=n,source_table_sha256=sha(table),
              backend_label=label,actual_binary_sha256=sha(binary),native_mode=mode,
              wall_median_seconds=median,
              wall_min_seconds=min(x["seconds"] for x in trial),
              wall_max_seconds=max(x["seconds"] for x in trial),
              n_steps=trial[0]["reported"]["n_steps"],
              result={k:trial[0]["reported"][k] for k in PHYSICAL_OUTPUTS},
              measurements=trial,one_warmup_seconds=warm["seconds"])
            cases.append(case)
            print("BENCHMARK_NATIVE",n,label,"median",f"{median:.6f}",
                  "steps",case["n_steps"],flush=True)
            (OUTPUT/"BENCHMARK_PARTIAL.json").write_text(json.dumps(dict(cases=cases,partial=True),indent=2)+"\n")
    summary=dict(schema="ns_cowling.native_dvode_perf_accuracy.v1",
        method="7 independent executable invocations per condition, native numerical computation only",
        cpu_affinity="unrestricted 24 logical cores; one native process per benchmark at a time",
        compiler="gfortran 13.3.0 -O2 -fcheck=all; original Netlib DVODE",
        original_physics_error="legacy k2 and Love Q invalid; compare old speed only, not old physics accuracy",
        numeric_source_identity="same direct model PK1 table density subsamples",
        cases=cases)
    for n in (257,1025,4097):
        base={z["backend_label"]:z for z in cases if z["source_table_rows"]==n}
        a=base["corrected_dvode_adams"];b=base["corrected_dvode_bdf"];rk=base["corrected_rk45"];old=base["legacy_binary"]
        a["speedups_at_same_fixed_source"]={
          "vs_corrected_rk45":rk["wall_median_seconds"]/a["wall_median_seconds"],
          "vs_legacy_unfixed_UNFAIR_PHYSICS":old["wall_median_seconds"]/a["wall_median_seconds"]}
        a["relative_difference_vs_dvode_bdf"]={k:(a["result"][k]-b["result"][k])/b["result"][k] for k in PHYSICAL_OUTPUTS}
        a["relative_difference_vs_tight_rk45"]={k:(a["result"][k]-rk["result"][k])/rk["result"][k] for k in PHYSICAL_OUTPUTS}
    (OUTPUT/"BENCHMARK.json").write_text(json.dumps(summary,indent=2)+"\n")
    return summary
if __name__=="__main__":
    data=benchmark()
    print("FINAL_OUTPUT",len(data["cases"]),"conditions",OUTPUT/"BENCHMARK.json")
