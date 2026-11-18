"""One independent read-only scientific audit of the frozen interpolation comparison.

This audits stored results against immutable H contracts; it does NOT rerun the
Maxwell root search or claim microscopic source regeneration.
"""
from pathlib import Path
from collections import Counter
import hashlib,json,math,numpy as np

R=Path(__file__).resolve().parents[1]
DATA=R/"data/matching_interpolation_20261118"
H=R/"data/low_density_crust_20261117/stellar/causal_extension"
get=lambda name:json.loads((DATA/name).read_text())
sh=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
main=get("RESULTS.json")
bern=get("BERNSTEIN_RESULTS.json")
adj=get("SHARP_NEIGHBOR_RESULTS.json")
sym=get("SYMBOLIC_VERIFICATION.json")
full=get("FULL_SUITE_RECEIPT.json")
hplan=json.loads((H/"TOV_PLAN.json").read_text())
original=json.loads((H/"TOV_RESULTS.json").read_text())
refs={z["id"]:z for z in original["cases"]}
shapes=("hermite","exponential","power","rational","loglog_polytrope")
cases={c["id"]:c for c in main["cases"]}
bern_cases={c["id"]:c for c in bern["cases"]}
adj_cases={c["id"]:c for c in adj["cases"]}
checks=[]
violations=[]
def check(name,condition,detail=""):
    checks.append(dict(test=name,pass_=bool(condition),detail=str(detail)))
    if not condition: violations.append(name+(": "+str(detail) if detail else ""))

check("eight frozen H contracts",len(hplan["cases"])==8==len(cases)==len(bern_cases)==len(adj_cases))
check("H source hashes unmodified",main.get("source_bytes_unchanged")==True and bern.get("source_bytes_unchanged")==True)
check("run statuses",main["status"]=="PASS_SCOPED" and bern["status"]=="PASS_SCOPED")
check("symbolic primitives",sym["all_zero"] and all(v=="0" for v in sym["identities"].values()))
check("full suite exit and code bytes",full["actual_exit_code"]==0 and full["source_bytes_unchanged"])
counts=Counter()
maxdiff=0
maxpressure=0
maxdeltaepsilon=0
for entry in hplan["cases"]:
    key=entry["id"]
    p=R/entry["contract"]
    check(key+" H SHA",sh(p)==entry["contract_sha256"])
    c=json.loads(p.read_text())
    i=c["upper"]["sharp_interfaces"][0]
    u=c["upper"]
    p0,p1=float(u["p"][i]),float(u["p"][i+1])
    mu0=(float(u["epsilon"][i])+p0)/float(u["n"][i])
    mu1=(float(u["epsilon"][i+1])+p1)/float(u["n"][i+1])
    n0,n1=float(u["n"][i]),float(u["n"][i+1])
    e0,e1=float(u["epsilon"][i]),float(u["epsilon"][i+1])
    maxpressure=max(maxpressure,abs(p1-p0))
    maxdeltaepsilon=max(maxdeltaepsilon,abs((e1-e0)-mu0*(n1-n0)))
    check(key+" Maxwell density direction",n1>n0)
    check(key+" Maxwell pressure equality",abs(p1-p0)<2e-12*max(1,abs(p0)))
    check(key+" Maxwell chemical equality",abs(mu1-mu0)<2e-12*max(1,abs(mu0)))
    check(key+" Maxwell latent heat",abs((e1-e0)-mu0*(n1-n0))<2e-9*max(1,abs(e1-e0)))
    case=cases[key]
    check(key+" four method statuses",all(m in case["methods"] for m in shapes))
    for m in shapes:
        counts[m+"_"+case["methods"][m]["status"]]+=1
    cbern=bern_cases[key]
    v=cbern["validation"]
    check(key+" Bernstein gate",v["status"]=="PASS" and cbern["degree"] in (4,6,8,12,16,24,32))
    check(key+" Bernstein monotonic coeff",min(cbern["bernstein_coefficients"])>=-1e-10)
    expected=cbern["side_cs2_reference"]
    found=v["endpoint_cs2_found"]
    check(key+" Bernstein endpoint speed match",np.allclose(expected,found,rtol=1e-8,atol=1e-12))
    for m in ("exponential","power","rational"):
        b=case["methods"][m]
        check(key+" "+m+" EOS validation",b["status"]=="PASS_SCOPED" and b["validation"]["status"]=="PASS")
        check(key+" "+m+" mass target",abs(b["star"]["target_mass"]["mass_msun"]-1.4)<1e-7)
        delta=abs(b["star"]["target_mass"]["radius_km"]-refs[key]["background"]["radius_km"])
        maxdiff=max(maxdiff,delta)
        check(key+" "+m+" first-law area",abs(b["EOS"]["independent_pressure_area_error_MeV_fm3"])<1e-9)
    check(key+" Bernstein mass target",abs(cbern["star"]["target_mass"]["mass_msun"]-1.4)<1e-7)
    check(key+" sharp neighbor H,p,n preserved",
          np.allclose([p0,n0,n1],
                      [adj_cases[key]["source_transition"][k] for k in ("p_t","n_low","n_high")],
                      rtol=0,atol=5e-14))
    check(key+" adjacent methods statuses",
          all(adj_cases[key]["methods"][m]["status"]=="PASS_SCOPED" for m in ("hermite","exponential","rational"))
          and adj_cases[key]["methods"]["power"]["status"]=="UNRESOLVED_OR_REJECTED")
check("expected lower Hermite six failures",counts["hermite_FAIL"]==6 and counts["hermite_PASS_SCOPED"]==2)
check("all other causal shape families pass",all(counts[m+"_PASS_SCOPED"]==8 for m in ("exponential","power","rational")))
check("naive log-log first-law endpoint negative controls",counts["loglog_polytrope_FAIL_ENDPOINT_ENERGY"]==8)
report=(R/"docs/matching_interpolation_20261118/RESULTS_KO.md").read_text()
check("report distinguishes fixed Maxwell roots","다시 탐색" in report or "root search" in report)
receipt=dict(status="PASS_SCOPED" if not violations else "FAIL",
              reference_base_commit="634a586fdd0c214bf7fa9718b02ebafa3e927008",
              independent_scope="frozen source/doublet identities and stored numeric results; no repeated physics solver, no independent Maxwell root discovery",
              checks_count=len(checks),failed_checks=violations,method_counts=dict(counts),
              max_abs_R_delta_km=maxdiff,max_abs_Maxwell_P_encoding_MeV_fm3=maxpressure,
              max_abs_Maxwell_latent_heat_residual_MeV_fm3=maxdeltaepsilon,
              checks=checks)
(DATA/"INDEPENDENT_AUDIT.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({k:v for k,v in receipt.items() if k!="checks"},indent=2))
raise SystemExit(0 if not violations else 1)
