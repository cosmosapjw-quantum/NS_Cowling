"""One bounded H-source comparison with C1-side-speed-matched causal Bernstein bridges."""
from pathlib import Path
import csv,hashlib,json,sys,time,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from matching_interpolation_shapes_20261118 import SlopeMatchedBernsteinInterval
from run_matching_interpolation_crosscheck_20261118 import star_case,physics_samples,endpoint_args
from compose_eos.whole_crust import WholeCrustBarotrope,CausalExponentialInterval
H=ROOT/"data/low_density_crust_20261117/stellar/causal_extension"
OUT=ROOT/"data/matching_interpolation_20261118"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(x):
    (OUT/"BERNSTEIN_RESULTS.json").write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def main():
    t0=time.time();plan=json.loads((H/"TOV_PLAN.json").read_text())
    refs={x["id"]:x for x in json.loads((H/"TOV_RESULTS.json").read_text())["cases"]}
    result=dict(status="RUNNING",source_commit="634a586fdd0c214bf7fa9718b02ebafa3e927008",
       scope="Only lower phenomenological bridge; exact source upper rows and Maxwell interface are fixed",
       degree_selection="smallest feasible Bernstein derivative degree among 4,6,8,12,16,24,32; optional roughness minimization",
       side_derivatives="fixed from previous E outer final polytrope and first source upper Hermite interval",
       cases=[])
    source_paths=[ROOT/"scripts/matching_interpolation_shapes_20261118.py",Path(__file__),
                  H/"TOV_PLAN.json",H/"TOV_RESULTS.json"]
    source_paths += [ROOT/e["contract"] for e in plan["cases"]]
    result["source_hashes_before"]={str(x.relative_to(ROOT)):sha(x) for x in source_paths}
    save(result)
    csvrows=[]
    for entry in plan["cases"]:
        p=ROOT/entry["contract"]
        assert sha(p)==entry["contract_sha256"]
        c=json.loads(p.read_text())
        b=WholeCrustBarotrope.from_contract(c)
        base=b.bridge
        cl=float(b.outer.gamma[-1]*base.p_left/(base.epsilon_left+base.p_left))
        interval=b.upper_intervals[0]
        cr=float(interval.number(0)/(interval.mu_left*interval.susceptibility(0)))
        args=endpoint_args(c)
        ref=CausalExponentialInterval(*args)
        mu=ref.mu_left+ref.width*np.linspace(0,1,301)
        bs=SlopeMatchedBernsteinInterval(*args,left_cs2=cl,right_cs2=cr)
        assert bs.validate()["status"]=="PASS"
        b.bridge=bs
        star=star_case(b,refs[entry["id"]],plan["settings"],
                    plan["tighter_settings"] if entry["label"]=="primary" else None)
        diag=physics_samples(bs,ref.evaluate_mu(mu))
        eq=refs[entry["id"]]["background"]["radius_km"]
        delta=star["target_mass"]["radius_km"]-eq
        record=dict(id=entry["id"],source_contract=entry["contract"],input_sha256=entry["contract_sha256"],
            n_bridge_right=entry["upper_start_density"],
            degree=bs.degree,selection=bs.selection,LP_attempts=bs.lp_attempts,
            bernstein_coefficients=bs.bernstein_coefficients.tolist(),
            side_cs2_reference=[cl,cr],validation=bs.validate(),physics=diag,
            star=star,delta_R_from_H_exponential_km=float(delta),
            source_sharp_interface_preserved=True)
        result["cases"].append(record);save(result)
        csvrows.append(dict(case=entry["id"],degree=bs.degree,cs2_max=diag["max_cs2_over_c2"],
            R_km=star["target_mass"]["radius_km"],delta_R_mm=delta*1e6,
            delta_P_max=diag["max_abs_delta_P_MeV_fm3"]))
        print("BERNSTEIN",entry["id"],"D",bs.degree,
              "R",f"{star['target_mass']['radius_km']:.10f}",
              "delta_R_mm",f"{delta*1e6:.5f}","maxcs2",f"{diag['max_cs2_over_c2']:.6f}",flush=True)
    result["source_hashes_after"]={str(x.relative_to(ROOT)):sha(x) for x in source_paths}
    result["source_bytes_unchanged"]=result["source_hashes_before"]==result["source_hashes_after"]
    result["status"]="PASS_SCOPED" if result["source_bytes_unchanged"] and len(result["cases"])==8 else "FAIL"
    result["runtime_seconds"]=time.time()-t0
    save(result)
    with (OUT/"BERNSTEIN_COMPARISON.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(csvrows[0]));w.writeheader();w.writerows(csvrows)
    print("FINISHED",result["status"],flush=True)
if __name__=="__main__":main()
