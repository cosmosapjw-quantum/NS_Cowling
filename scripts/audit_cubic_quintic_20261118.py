"""Bounded read-only numerical/source acceptance for cubic/quintic campaign."""
from pathlib import Path
from collections import Counter
import json,hashlib,math
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"data/cubic_quintic_matching_20261118"
SOURCE=ROOT/"data/native_dense_interpolation_20261118"
F=ROOT/"docs/cubic_quintic_matching_20261118"
load=lambda p:json.loads(p.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
s=load(D/"SUMMARY.json")
checks=[];fail=[]
def test(label,pred):
    condition=bool(pred)
    checks.append(dict(label=label,passed=condition))
    if not condition:fail.append(label)
test("10 native fixed-model direct sources",s["reference_density_models"]==10 and s["native_direct_state_count"]==25612)
test("5 N128 fixed TF native source families",s["tf_branch_models"]==5 and s["tf_holdout_count"]==30)
test("exactly two displayed methods",s["method_names"]==["cubic","quintic"])
test("58 phase/method source holdout runs",len(s["dense_model_phase_runs"])==58)
test("20 native construction outputs",len(s["native_smooth_css"])==20)
test("16 native Maxwell sparse-root results",len(s["maxwell"])==16)
test("paired stellar outputs",len(s["stellar_paired"])==10)
test("all method generation finishes",all(z["status"]=="CONSTRUCTED_VERIFIED_NATIVE" for z in s["native_smooth_css"]))
test("matched sharp CSS count",sum(z["css_sharp"] for z in s["native_smooth_css"])==4)
test("all native accuracy phase runs admitted",all(z["physical_fraction"]>=0 and z["physical_fraction"]<=1 for fam in s["family_levels"].values() for methods in fam.values() for z in methods.values()))
for family in ("rmf","core"):
    for p in sorted((D/"dense_cv_v2"/family).glob("*/RESULT.json")):
        x=load(p)
        for q in x["runs"]:
            for m in ("cubic","quintic"):
                z=q["methods"][m]
                test(x["model"]+" "+str(q["training_stride"])+" "+m+" CV",z["status"]=="CHECKED_NATIVE_HOLDOUT" and z["no_dense_truth_leakage"])
                test(x["model"]+" "+m+" Euler",z["max_abs_euler_residual_mev_fm3"]<1e-7)
        native=SOURCE/family/x["model"]
        source=native/"table/compose/eos.thermo" if family=="rmf" else native/x["model"]/"eos.thermo"
        test(x["model"]+" source digest",x["source_sha256"]==sha(source))
for p in sorted((D/"tf_cv").glob("*.json")):
    if p.name=="CAMPAIGN.json":continue
    x=load(p);test("TF "+x["case"]+" native count",x["n_holdout"]>0)
    test("TF "+x["case"]+" no train/hold overlap",set(x["training_ids"]).isdisjoint(x["heldout_ids"]))
    test("TF "+x["case"]+" cubic/quintic",all(x["methods"][m]["status"]=="CHECKED_NATIVE_TF_HOLDOUT" for m in ("cubic","quintic")))
    test("TF "+x["case"]+" native source no interpolation truth",x["direct_native_reference"].startswith("N128 same droplet branch"))
for p in sorted((D/"native_matching").glob("*/RESULT.json")):
    x=load(p)
    test(p.parent.name+" all five native methods",set(x["methods"])=={"cubic_F","quintic_F","cubic_P","quintic_P","css"})
    for method,z in x["methods"].items():
        test(p.parent.name+":"+method+" serialized check",z["status"]=="CONSTRUCTED_VERIFIED_NATIVE" and
             max(z["serialized_readback_max_abs"].values())<1e-9 and z["native_report"]["consistency"]["overall"]["pass_"])
        for name,digest in z["output_files_sha256"].items():
            test(p.parent.name+":"+method+":"+name,sha(p.parent/method/"compose"/name)==digest)
for p in sorted((D/"maxwell").glob("*/RESULT.json")):
    x=load(p);native=x["native_full_grid_maxwell"]
    should="PK1_CANONICAL" in p.parent.name
    test(p.parent.name+" native coexistence status",(native["status"]=="PASS_SCOPED")==should)
    for kind in ("cubic","quintic"):
        for level in ("sparse_4x","sparse_8x"):
            z=x["fits"][kind+"_"+level]
            test(p.parent.name+" "+kind+level+" root count",(z["crossings"]==1)==should)
            if should:
                q=z["roots"][0]
                test(p.parent.name+" "+kind+level+" positive jump",q["positive_density_jump"])
                test(p.parent.name+" "+kind+level+" Maxwell latent identity",
                     abs(q["delta_epsilon_MeV_fm3"]-q["mu_MeV"]*q["density_jump"])<1e-12)
for p in sorted((D/"stellar_bounded").glob("*/RESULT.json")):
    x=load(p)
    test(x["case"]+" source transition unchanged",x["source_transition_unchanged"])
    for kind in ("cubic","quintic"):
        q=x["cases"][kind]
        test(x["case"]+" "+kind+" same mass",abs(q["background"]["mass_msun"]-1.4)<2e-7)
        for key,mode in q["modes"].items():
            if mode["status"]=="PASS":
                expect=(0,0) if key.endswith("_f") else (1,1)
                if key.endswith("_p1"):expect=(1,1)
                test(x["case"]+" "+kind+key+" node count",
                     (mode["W_sign_changes_including_interfaces"],mode["V_sign_changes_including_interfaces"])==expect)
for z in s["stellar"]:
    if z["mode_f"]["status"]=="PASS":
        test(z["case"]+z["method"]+" positive f",z["mode_f"]["frequency_hz"]>0)
    if z["mode_p1"]["status"]=="PASS":
        test(z["case"]+z["method"]+" positive p1",z["mode_p1"]["frequency_hz"]>0)
test("report authored",(F/"RESULTS_KO.md").exists())
if (D/"FIGURES_RECEIPT.json").exists():
    fig=load(D/"FIGURES_RECEIPT.json")
    test("five plot outputs",len(fig["figure_files"])==5)
    for q in fig["figure_files"]:
        test("figure identity "+q["path"],sha(ROOT/q["path"])==q["sha256"])
        test("figure uses exactly requested methods",fig["methods"]==["cubic_hermite","quintic_hermite"])
rec=dict(status="PASS_SCOPED" if not fail else "FAIL",checks=len(checks),failed=fail,
    scope="five native physical models plus five RMF, five TF source, 20 native matching, 4 original-source Maxwell pairs; only scoped reliable modal roots accepted",
    source_head="2b32672da2b9b602942a362ff7208510708be459",
    unresolved=["quintic_F generated output not accepted by whole-star upper interval Hermite consumer on PK1",
                "not all p1 roots identified across complete spectrum",
                "finite density source and selected crust geometry, not global minimum",
                "original native input itself has limitations at high density (APR)"],
    checks_detail=checks)
(D/"INDEPENDENT_AUDIT.json").write_text(json.dumps(rec,indent=2)+"\n")
print("AUDIT",rec["status"],"checks",len(checks),"failed",len(fail),fail[:12])
raise SystemExit(0 if not fail else 1)
