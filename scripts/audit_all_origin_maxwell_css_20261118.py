"""Independent read-only audit: sealed origin identity, 635 exact pair cases, local
Hermite roots, native Maxwell/CSS checks and conditional full-star TOV.

A PASS certifies only enumerated source scopes, not all possible EOS physically
or nonexistence of roots on unadmitted/missing external sources.
"""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import hashlib,json,math,csv
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"data/all_origin_maxwell_css_20261118"
load=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
catalog=load(D/"SOURCE_CATALOG.json")
data=load(D/"CAMPAIGN.json")
roots=load(D/"HERMITE_ROOT_CROSSCHECK.json")
whole=load(D/"whole_stars/CAMPAIGN.json")
interp=load(D/"HERMITED_INPUT_2X2/CAMPAIGN.json")
assert interp["matchings"]==["maxwell","css"] and interp["interpolation_orders"]==["cubic","quintic"]
assert interp["inputs_count"]==7 and interp["variants_expected"]==28
summ=load(D/"RESULTS_SUMMARY.json")
src={x["source_id"]:x for x in catalog["sources"]}
checks=[];fails=[]
def check(label,expr):
    condition=bool(expr)
    checks.append(dict(label=label,pass_=condition))
    if not condition:fails.append(label)
check("source count",len(src)==95==summ["number_sources"])
check("source crust count",sum(x["role"]=="crust" for x in src.values())==70)
check("source core count",sum(x["role"]=="core" for x in src.values())==25)
check("candidate count",len(data["cases"])==635==summ["number_pairs_total"])
check("candidate origin purity",all(q["tier"]==src[q["crust_source_id"]]["tier"]==src[q["core_source_id"]]["tier"] for q in data["cases"]))
check("unique source pair IDs",len(set(q["pair_id"] for q in data["cases"]))==635)
check("source type labels",set(x["tier"] for x in src.values())=={"original","corrected","reconstructed"})
counts=Counter(x["both_status"] for x in data["cases"])
check("scoped result categories",counts==Counter({"NOT_BOTH":350,"BLOCKED_SOURCE":278,"FINITE_ONLY_NOT_WHOLE":7}))
check("strict whole negative",summ["true_original_full_whole_success"]==0)
check("five conditional complete-star TOV",whole["both_conditional_whole"]==5)
check("7 finite dual CSS/Maxwell",summ["finite_dual_success"]==7)
check("Hermite local roots seven",roots["eligible_dual_pair_count"]==7)
check("matching methods only Maxwell/CSS",summ["matching_methods"]==["maxwell","css"])
check("input interpolants only cubic/quintic",summ["input_interpolation_methods"]==["cubic","quintic"])
check("native 2x2 input variants total",summ["interpolation_matrix_2x2_total"]==28)
check("actual native Maxwell on interpolated inputs",summ["actual_native_maxwell_interpolated_source_accepted"]==14)
check("actual native CSS on interpolated inputs",summ["actual_native_css_interpolated_source_accepted"]==14)
check("number original cohort native Maxwell cases",len(load(ROOT/"data/maxwell_three_families_20261117/campaign_v2/CAMPAIGN.json")["cells"])==195)
check("user ZIP explicitly in census",all("user_zip:"+name in src for name in
     ("MYN13_RHF_TF:crust","OMEG1_EoS:core")))
check("strict source absence on ZIP",all(src["user_zip:"+name]["status"]=="BLOCKED_NO_REMOTE_VERIFIED_COMPOSE_PROJECTION" for name in
     ("MYN13_RHF_TF:crust","OMEG1_EoS:core")))
check("NV pressure absence correctly blocked",src["external_nuda:1973-Negele-Vautherin:crust"]["status"]=="BLOCKED_MISSING_PRESSURE")
check("nuda wheel has SHA",len(catalog["nuda_evidence"]["sha256"])==64)
nuda=Path("/home/cosmosapjw/ns_cowling_wt/scratch/all_origin_maxwell_css_20261118/nuda/nucleardatapy-1.0.2-py3-none-any.whl")
check("nuda original wheel identity",nuda.is_file() and sha(nuda)==catalog["nuda_evidence"]["sha256"])
check("source catalog exact type IDs",len(src)==len(catalog["sources"]))
# All sealed original/new/native paths, including public originals and reconstructed
# projection provenance hashes. Unknown external ZIP rows have no shell source.
source_count=0;sha_count=0
for x in catalog["sources"]:
    if x["status"]!="SOURCE_SNAPSHOT_AVAILABLE":continue
    source_count+=1
    source_dir=ROOT/x["source_path"]
    expect=x.get("source_table_sha256",{})
    if not expect:continue
    for name,digest in expect.items():
        p=source_dir/name
        check("raw source identity: "+x["source_id"]+" "+name, p.exists() and sha(p)==digest)
        sha_count+=1
dual=[q for q in data["cases"] if q["both_status"]=="FINITE_ONLY_NOT_WHOLE"]
local={(z["tier"],z["crust"],z["core"]):z for z in roots["cases"]}
check("Hermite exact seven pairing identities",len(local)==7 and
    set(local)==set((q["tier"],q["crust"],q["core"]) for q in dual))
for q in dual:
    k=(q["tier"],q["crust"],q["core"])
    r=local[k]
    check("different functional for "+str(k),not q["same_functional"])
    check("CSS actually native admitted "+str(k),
         q["css"]["status"]=="PASS" and 0<q["css"]["cs2"]<=1 and
         q["css"]["n_low"]<q["css"]["n_high"] and q["css"]["P_t"]>0)
    check("Maxwell actually native admitted "+str(k),
         q["maxwell"]["status"] in ("PASS_NEW_NATIVE_FINITE","PASS_VERIFIED_FINITE"))
    # Matching construction remains Maxwell or CSS. Hermite is the input EOS interpolation.
    p2x=(D/"HERMITED_INPUT_2X2"/q["tier"]/
         (q["crust"]+"__"+q["core"])/"RESULT.json")
    variant=load(p2x)["interpolation_pair_results"]
    check("exactly four interpolation×matching variants "+str(k),
          set(variant)=={"cubic+maxwell","quintic+maxwell","cubic+css","quintic+css"})
    for order in ("cubic","quintic"):
        qmax=variant[order+"+maxwell"]
        qcss=variant[order+"+css"]
        check("Maxwell local P(mu) interpolant root "+order+str(k),
              qmax["status"]=="PASS_SCOPED_LOCAL_MAXWELL_ROOT" and len(qmax["roots"])==1)
        check("CSS native on interpolated EOS "+order+str(k),
              qcss["status"]=="PASS_DERIVED_INTERPOLATED_SOURCE" and
              qcss["matching"]=="css" and qcss["actual_native_matching"]=="css" and
              qcss["P_t"]>0 and qcss["density_jump"]>0 and
              0<qcss["cs2_css"]<=1 and qcss["serialized_readback_abs"]<1e-9 and
              qcss["source_crust_interpolant"]["inserted_rows"]>0 and
              qcss["source_core_interpolant"]["inserted_rows"]>0)
    for m in ("cubic","quintic"):
        rr=r["methods"][m]
        check("Hermite source-local root count "+m+str(k),len(rr["roots"])==1)
        for z in rr["roots"]:
            check("Hermite local root physical "+m+str(k),
                  z["n_low"]<z["n_high"] and z["p_t"]>0
                  and len(z["c2_sides"])==2 and 0<min(z["c2_sides"])<=max(z["c2_sides"])<=1
                  and abs(z["latent_heat"]-z["mu_t"]*z["density_jump"])<1e-9)
            check("Hermite local root equality "+m+str(k),
                z["P_abs_residual"]<1e-8*max(1,z["p_t"]))
    if q["tier"]=="reconstructed" and q["crust"].startswith("SPG_"):
        check("SPG native originally no outer "+str(k),not src[q["crust_source_id"]]["outer_present"])
# Numeric TOV after explicit same-code E outer completion, not true original
w_root=D/"whole_stars"
for crust in ("SPG_M1","SPG_M2","SPG_M3","SPG_M4","SPG_M5","GPPVA_FSU2"):
    rec=load(w_root/crust/"RESULT.json")
    check(crust+" source/whole method pair",set(rec["methods"])=={"css","maxwell"})
    for method,q in rec["methods"].items():
        if crust=="GPPVA_FSU2":
            check("FSU2 phase restriction retained "+method,q["status"]=="BLOCKED_OR_REJECTED")
        else:
            check(crust+" "+method+" whole conditional",q["status"]=="CONDITIONAL_WHOLE_NS_VERIFIED_TOV")
            x=q["admitted"]
            check(crust+" "+method+" TOV 1.4",abs(x["stellar"]["mass_msun"]-1.4)<2e-7)
            check(crust+" "+method+" R physical",9<x["stellar"]["R_km"]<20)
            check(crust+" "+method+" source Maxwell/CSS jump",x["density_jump"]>0 and x["latent_heat"]>0)
            check(crust+" "+method+" same source bridge",x["lower_bridge_validation"]["status"]=="PASS")
            check(crust+" "+method+" no original source outer claim",not rec["source_low_sha256"]=="")
check("whole original provider no false PASS",all(x["both_status"]!="BOTH_WHOLE_PASS" for x in data["cases"]))
# Machine-readable CSV is complete and ordered.
for name,n in (("ALL_CANDIDATE_PAIRS.csv",635),("DUAL_FINITE_SUCCESS.csv",7),
               ("CONDITIONAL_WHOLE_TOV.csv",5),("INDEPENDENT_CRUST_INVENTORY.csv",70)):
    with (D/name).open(newline="") as f:r=list(csv.DictReader(f))
    check("CSV row identity "+name,len(r)==n)
report=(ROOT/"docs/all_origin_maxwell_css_20261118/RESULTS_KO.md").read_text()
check("report current exact counts",all(s in report for s in ("635 same-origin","95 distinct","7 pair","M1","M5","0.")))
output=dict(status="PASS_SCOPED" if not fails else "FAIL",
            checked=len(checks),failed=fails,source_sha_count=sha_count,materialized_source_count=source_count,
            captured_input_head="bc09c87010fccb280c455d765e6239d3910fb0db",
            evidence_scope="exact materialized input bytes; 635 same-origin cases; 7 native finite maxwell+CSS pairings; 7*2 local cubic/quintic roots; 5 conditional whole TOV",
            explicit_limitations=["No strict original-inner+outer fully admitted EOS",
                                  "not all external nuda raw models can be read as thermodynamically complete cold EOS",
                                  "screen negatives outside legacy Fortran campaign not formal nonexistence proofs",
                                  "CSS and Maxwell are different physical constructions, not interchangeable",
                                  "candidate code-generated reconstructed crust full-star BSk22 outer is phenomenological anchored closure, not a single unified EDF"],
            checks=checks)
(D/"INDEPENDENT_AUDIT.json").write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n")
print("AUDIT",output["status"],"CHECKED",len(checks),"FAILED",len(fails),
      "SOURCE_SHA",sha_count,"MATERIALIZED_SOURCES",source_count, "FIRST_FAILURES",fails[:12])
raise SystemExit(0 if not fails else 1)
