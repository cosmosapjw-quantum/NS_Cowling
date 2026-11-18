"""Bounded read-only independent audit of direct native holdout evidence."""
from pathlib import Path
from collections import defaultdict
import json,hashlib,math,numpy as np
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"data/native_dense_interpolation_20261118"
SUMMARY=json.loads((D/"aggregated/SUMMARY.json").read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
errors=[]; checks=0
def check(label,pred):
    global checks
    checks+=1
    if not pred:errors.append(label)
check("all source counts",SUMMARY["generated_direct_native_state_count"]==25612)
check("TF density count",SUMMARY["TF_distinct_new_withheld"]==30)
check("figure method excludes linear","linear" not in SUMMARY["interpolation_methods_plotted"])
check("user method excluded",SUMMARY["interpolation_methods_plotted"]==["hermite","exponential","rational","bernstein","power","pchip","akima","cubic"])
for family in ("rmf","core"):
    for modeldir in sorted((D/family).iterdir()):
        if not modeldir.is_dir():continue
        model=modeldir.name
        if family=="rmf":
            x=json.loads((modeldir/"RESULT.json").read_text())
            output=modeldir/"table/compose"
            check(model+" raw status",x["status"]=="GENERATED_VERIFIED_DIRECT_NATIVE")
            expected=x["table_sha256"]
            check(model+" model state count",x["grid_count"]==4097)
            check(model+" native residual",x["direct_native_residual"]["physical_state_checks"])
        else:
            x=json.loads((modeldir/model/"RECEIPT.json").read_text())
            output=modeldir/model
            expected=x["output_sha256"]
            check(model+" core native readback",x["native_readback"]=="PASS")
            check(model+" valid source range",max(r["n_b"] for r in x["rows"])<=(.8+1e-12 if model=="BBB_BHF2" else 1+1e-12))
        check(model+" CompOSE file identity",all(sha(output/k)==v for k,v in expected.items()))
        cvp=D/"cv"/family/model/"RESULT.json"
        cv=json.loads(cvp.read_text())
        check(model+" cv status",cv["status"]=="COMPLETED_SCOPED_WITH_NATIVE_HOLDOUT")
        check(model+" source table exact",cv["native_source_sha256"]==sha(output/"eos.thermo"))
        for run in cv["runs"]:
            if run.get("status")!="COMPLETE_WITH_REJECTED_METHODS":continue
            check(model+" method set",set(run["metrics"])==set(SUMMARY["interpolation_methods_plotted"])|{"linear"})
            check(model+" holdout finite",run["training_size"]+run["heldout_size"]==run["source_row_range"][1]-run["source_row_range"][0])
            for method,data in run["metrics"].items():
                check(model+" "+method+" result",data["status"]=="COMPUTED_WITH_NATIVE_HOLDOUT")
                check(model+" "+method+" coverage valid",0<=data["finite_prediction_fraction"]<=1)
                check(model+" "+method+" physical valid",0<=data["physical_among_native_causal"]<=1)
                if data["number_density_relative"]["count"]:
                    check(model+" "+method+" nRMS finite",math.isfinite(data["number_density_relative"]["rms"]))
        if model=="APR_APR":
            check("APR has exactly 6/8/7",[int(z["phase"]) for z in cv["runs"][::2]]==[6,8,7])
            check("APR native high density acausality flagged",sum(z["native_physics"]["count_acausal"] for z in cv["runs"])>0)
original=json.loads((D/"tf/PARALLEL_RECEIPT.json").read_text())
reseed=json.loads((D/"tf/refined_source_bracket/PARALLEL_RECEIPT.json").read_text())
check("TF five original campaigns",len(original["results"])==5)
check("TF original run count",sum(v["requested"] for v in original["results"])==57)
check("TF new first-try accepted",sum(v["accepted"] for v in original["results"])==26)
check("TF recovered accepted",sum(v["recovered"] for v in reseed["models"])==7)
check("TF rejected preserved",sum(v["still_rejected"] for v in reseed["models"])>0)
for case in ("DDME2","DDMEX","NL3","NL3wrL55","TM1e"):
    tf=json.loads((D/"tf_holdout"/f"{case}.json").read_text())
    check(case+" only accepted holdout",tf["status"]=="SCOPED_FINITE_NATIVE_VALIDATION")
    check(case+" no out-of-source predictions",len(tf["excluded_direct_states"])==0)
    check(case+" train/holdout separate",set(tf["training_shas"]).isdisjoint(tf["new_direct_shas"]))
    for method,value in tf["metrics"].items():
        check(case+" "+method+" native status",value["status"]=="EVALUATED_AGAINST_NEW_NATIVE_DROPLET_N128")
        check(case+" "+method+" physical fraction range",0<=value["physical_among_native_causal"]<=1)
    check(case+" no rejected native promoted",tf["heldout_distinct"]==len(tf["new_direct_shas"]))
check("user-visible figure methods",all("linear" not in x.lower() for x in json.loads((D/"aggregated/FIGURES_RECEIPT.json").read_text())["methods"]))
figure=json.loads((D/"aggregated/FIGURES_RECEIPT.json").read_text())
check("all 5 charts generated",len(figure["figures"])==5)
for item in figure["figures"]:
    path=ROOT/item["relative_path"]
    check("figure bytes "+path.name,path.is_file() and sha(path)==item["sha256"])
# The model-level ranking must not include a method lacking overall physical coverage.
for family,methods in SUMMARY["family_summary"].items():
    for method,v in methods.items():
        if v["coverage_fraction"]<.995 or v["physical_fraction"]<.995:
            check("rank exclusion "+family+method,not v["qualifies_for_global_accuracy_rank"])
result=dict(status="PASS_SCOPED" if not errors else "FAIL",
   independent_evidence_scope="read-only source SHA/readback, TF acceptance, phase split, holdout eligibility; does not reproduce all native states",
   check_count=checks,failures=errors,
   source_stage="native-dense-interpolation-benchmark-20261118",
   reference_head="a0fa5fb93851a9647030f3589c24a62ddea527b1",
   limitations=["no new Maxwell root search","N128 droplet only TF","BBB baseline reconstructs published table","native spinodal/global minima unverified"],
   conclusion="sufficient for scoped fixed-model native-holdout numeric comparison only" if not errors else "do not publish before fixing findings")
(D/"aggregated/INDEPENDENT_AUDIT.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps(result,indent=2))
raise SystemExit(0 if not errors else 1)
