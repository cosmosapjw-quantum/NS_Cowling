"""Probe interpolation sensitivity immediately adjacent to fixed Maxwell doublets.

The coexistence jump itself is NEVER interpolated. Each side's source
(n,P,epsilon) knot values are unchanged; only intervals j-1 and j+1 vary.
"""
from pathlib import Path
import hashlib,json,sys,csv,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from matching_interpolation_shapes_20261118 import CausalShapeInterval
from compose_eos.whole_crust import HermiteInterval,CausalExponentialInterval,WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background
DATA=ROOT/"data/matching_interpolation_20261118"
H=ROOT/"data/low_density_crust_20261117/stellar/causal_extension"
def save(x):
    (DATA/"SHARP_NEIGHBOR_RESULTS.json").write_text(json.dumps(x,indent=2,allow_nan=False,ensure_ascii=False)+"\n")
def endpoints(c,ii):
    u=c["upper"]
    return tuple(u[k][index] for index in [ii,ii+1] for k in ("n","p","epsilon"))
def at_side(bridge,at):
    mu=bridge.mu_left+bridge.width*at
    return float(bridge.number(at)/(mu*bridge.susceptibility(at)))
def main():
    plan=json.loads((H/"TOV_PLAN.json").read_text());T=json.loads((H/"TOV_RESULTS.json").read_text())
    refs={x["id"]:x for x in T["cases"]}
    settings=plan["settings"]
    out=dict(scope="Upper intervals immediately left/right of sharp interface only; no smoothing of Maxwell density gap",
       source_H_commit="634a586fdd0c214bf7fa9718b02ebafa3e927008",cases=[],status="RUNNING")
    start=time.time();save(out)
    rows=[]
    for c in plan["cases"]:
        contract=json.loads((ROOT/c["contract"]).read_text())
        bar=WholeCrustBarotrope.from_contract(contract)
        j=min(bar.sharp_interfaces)
        assert j>0 and j+2<len(bar.upper_n)
        p_t=float(bar.upper_p[j])
        mu_t=float(bar.upper_mu[j]);mu_plus=float(bar.upper_mu[j+1])
        record=dict(id=c["id"],sharp_index=j,source_transition=dict(
            p_t=p_t,mu_low=mu_t,mu_high=mu_plus,n_low=float(bar.upper_n[j]),
            n_high=float(bar.upper_n[j+1]),density_jump=float(bar.upper_n[j+1]-bar.upper_n[j])),
            methods={},source_sha256=c["contract_sha256"])
        base_R=refs[c["id"]]["background"]["radius_km"]
        ncenter=refs[c["id"]]["background"]["central_density_fm3"]
        for name,kind in [("hermite",HermiteInterval),("exponential",CausalExponentialInterval),
                          ("power",None),("rational",None)]:
            try:
                lb=endpoints(contract,j-1);rb=endpoints(contract,j+1)
                obj=lambda x:kind(*x) if kind is not None else CausalShapeInterval(*x,family=name)
                low=obj(lb);high=obj(rb)
                vl,vh=low.validate(),high.validate()
                if vl["status"]!="PASS" or vh["status"]!="PASS":
                    raise ValueError("adjacent interval failed finite causal gate")
                bl=bar.upper_intervals[j-1]
                bh=bar.upper_intervals[j+1]
                lows=[at_side(i,1.) for i in [bl,low]]
                highs=[at_side(i,0.) for i in [bh,high]]
                probe=np.linspace(0,1,101)
                differences={}
                for label,original,replacement in [("low",bl,low),("high",bh,high)]:
                    mu=original.mu_left+original.width*probe
                    sr=original.evaluate_mu(mu)
                    sa=replacement.evaluate_mu(mu)
                    differences[label]=dict(
                         max_dcs2=float(np.max(abs(np.asarray(sa["cs2_over_c2"])-sr["cs2_over_c2"]))),
                         max_dP=float(np.max(abs(np.asarray(sa["p"])-sr["p"]))),
                         max_dE=float(np.max(abs(np.asarray(sa["epsilon"])-sr["epsilon"]))))
                changed=WholeCrustBarotrope.from_contract(contract)
                changed.upper_intervals[j-1]=low
                changed.upper_intervals[j+1]=high
                assert changed.sharp_interfaces==bar.sharp_interfaces
                assert np.array_equal(changed.upper_n,bar.upper_n)
                assert np.array_equal(changed.upper_p,bar.upper_p)
                assert np.array_equal(changed.upper_epsilon,bar.upper_epsilon)
                eos=SharpContractEnthalpyEOS(changed)
                bg=solve_sharp_background(eos,ncenter,**settings)
                jump=bg.interfaces[0]
                assert abs(jump["p_MeV_fm3"]-p_t)<1e-11
                result=dict(status="PASS_SCOPED",cs2_low_side=lows[1],
                    cs2_high_side=highs[1],baseline_cs2_low_side=lows[0],
                    baseline_cs2_high_side=highs[0],variation=differences,
                    same_central_density=True,R_km=float(bg.radius_km),
                    mass_msun=float(bg.mass_msun),delta_R_from_H_km=float(bg.radius_km-base_R))
            except Exception as ex:
                result=dict(status="UNRESOLVED_OR_REJECTED",exception_type=type(ex).__name__,
                            diagnostic=str(ex))
            record["methods"][name]=result
            rows.append(dict(case=c["id"],method=name,status=result["status"],
                cs2_low=result.get("cs2_low_side"),cs2_high=result.get("cs2_high_side"),
                delta_R_from_H_km=result.get("delta_R_from_H_km")))
            print("INTERFACE",c["id"],name,result["status"],
                  result.get("delta_R_from_H_km",""),flush=True)
        out["cases"].append(record);save(out)
    out["status"]="COMPLETED_WITH_REJECTED_CANDIDATES" if any(
        m["status"]!="PASS_SCOPED" for c in out["cases"] for m in c["methods"].values()) else "PASS_SCOPED"
    out["duration_seconds"]=time.time()-start
    save(out)
    with (DATA/"SHARP_NEIGHBOR.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print("FINISHED",out["status"],flush=True)
if __name__=="__main__":main()
