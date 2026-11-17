#!/usr/bin/env python3
"""Cross-check native PCHIP g(log P) Cowling against the exact source law."""
import csv,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from dataclasses import asdict
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python_frontend'))
from compose_eos.stellar_background import ContractEnthalpyEOS,solve_background,MEV_FM3_TO_PA
from compose_eos.stellar_cowling import CowlingSettings,solve_mode
from compose_cowling.solver import _compile_solver
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/stellar_consumers_20261117/cowling"
def dump(name,x): (OUT/name).write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def native_node_counts(prefix):
    modes={}
    for j,label in [(1,"f"),(2,"p1")]:
        counts={}
        for region in ["inner","outer"]:
            rows=list(csv.DictReader(Path(str(prefix)+f"_mode{j}_{region}.csv").open()))
            counts[region]={}
            for col in ["w_over_r3","v_over_r2"]:
                a=np.array([float(r[col]) for r in rows]);a=a[abs(a)>1e-12*max(abs(a))]
                counts[region][col]=int(np.count_nonzero(a[1:]*a[:-1]<0))
        counts["total_W_nodes"]=sum(counts[r]["w_over_r3"] for r in ["inner","outer"])
        counts["total_V_nodes"]=sum(counts[r]["v_over_r2"] for r in ["inner","outer"])
        modes[label]=counts
    return modes
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    tov=ROOT/"data/stellar_consumers_20261117/tov"
    canonical=json.loads((tov/"CANONICAL_1P4.json").read_text())
    eos=ContractEnthalpyEOS.from_path(ROOT/"data/core_support_20261117/primary_v1/REPRESENTATION.json")
    nc=canonical["central_density_fm3"];hcut=1e-4
    bg=solve_background(eos,nc,surface_h=hcut,rtol=1e-10,atol=1e-12,max_step_h=.002)
    pc=eos.state_at_enthalpy(bg.h_c)["p"]*MEV_FM3_TO_PA
    pcut=eos.state_at_enthalpy(hcut)["p"]*MEV_FM3_TO_PA
    settings=CowlingSettings(surface_h=hcut,rtol=1e-9,atol=1e-12)
    references=[]
    for label,br in [("f",[2250,2500]),("p1",[6500,6750])]:
        x=solve_mode(bg,br,settings);prof=x.pop("profile")
        with (OUT/("exact_hcut1e-4_"+label+"_profile.csv")).open("w") as f:
            w=csv.DictWriter(f,fieldnames=prof[0]);w.writeheader();w.writerows(prof)
        references.append(x)
    dump("EXACT_NATIVE_CUTOFF.json",dict(background=bg.summary(),settings=asdict(settings),modes=references))
    _,executable=_compile_solver()
    records=[]
    for sub in [1,4]:
        table=tov/("native_exact_samples_sub%d.csv"%sub)
        rows=list(csv.DictReader(table.open()))
        mutarget=eos.a*np.exp(.15*bg.h_c)
        index=1+int(np.argmin([abs(float(r["mu_b_mev"])-mutarget) for r in rows]))
        prefix=OUT/("native_sub%d"%sub)
        nml=prefix.with_suffix(".nml")
        nml.write_text("&precision\n"+
          f"table_path='{table}'\ntask='spectrum'\nmethod='dopri'\nshooting='matched'\n"+
          f"match_index={index}\nsplit_index=1\nalpha=0.0\npc_pa={pc:.17e}\np_cut_pa={pcut:.17e}\n"+
          "r0_m=1.0\ntarget_mass_kg=0.0\nrtol=1.0e-10\nmode_rtol=1.0e-9\nstep_cap=0.05\n"+
          "frequency_lo=2250.0,6500.0\nfrequency_hi=2500.0,6750.0\n"+
          "root_atol_hz=1.0e-6\nresidual_tol=2.0e-6\nnmode=2\n/\n")
        start=time.monotonic()
        run=subprocess.run([str(executable),str(nml),str(prefix)],capture_output=True,text=True)
        prefix.with_suffix(".stdout.log").write_text(run.stdout)
        prefix.with_suffix(".stderr.log").write_text(run.stderr)
        record=dict(subdivisions=sub,table_sha256=hashlib.sha256(table.read_bytes()).hexdigest(),
                    executable=str(executable),executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                    command=[str(executable),str(nml),str(prefix)],exit_code=run.returncode,
                    seconds=time.monotonic()-start,match_index=index,
                    interpolation="monotone_PCHIP_epsilon_minus_pressure_of_logP_with_firstlaw_density_v1",
                    provider_source="oscillation_source/solutions/src/eos/ns_table_eos.f90",
                    provider_source_sha256=hashlib.sha256((ROOT/"oscillation_source/solutions/src/eos/ns_table_eos.f90").read_bytes()).hexdigest())
        if prefix.with_suffix(".json").exists():
            result=json.loads(prefix.with_suffix(".json").read_text())
            record["result"]=result
            if run.returncode==0: record["native_profile_nodes"]=native_node_counts(prefix)
            record["frequency_minus_exact_cutoff_hz"]=[
                result["f_hz"]-references[0]["frequency_hz"],
                result["p1_hz"]-references[1]["frequency_hz"]]
        records.append(record);dump("NATIVE_COMPARISON.json",records)
        print(json.dumps(record),flush=True)
if __name__=="__main__":main()
