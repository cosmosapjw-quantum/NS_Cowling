#!/usr/bin/env python3
"""C9 selective final readback: exact sources, 9 fine and 7 coarse outputs.
No native kernel invocation; uses raw serialized files and saved numerical evidence.
"""
from pathlib import Path
import hashlib, json, re
from collections import Counter
import numpy as np

ROOT=Path("/home/cosmosapjw/ns_cowling_wt/new_model_maxwell_20261117")
BASE=ROOT/"data/new_model_maxwell_20261117"
OUT=BASE/"independent_review"
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def path(p):return Path(p) if Path(p).is_absolute() else ROOT/p
def num(s):return float(s.replace("D","E").replace("d","e"))
def intern(p):
    p=Path(p)
    d=np.genfromtxt(p/"eos_data.csv",delimiter=",",names=True)
    meta=(p/"eos_data.nml").read_text()
    mass=num(re.search(r"(?im)^\s*mn\s*=\s*([^\s,]+)",meta).group(1))
    return d,mass,meta
def physical(d,mass):
    n=d["nb"]
    return n,n*d["Q1"],n*mass*(1+d["Q7"]),mass*(1+d["Q3"])
def scaled(a,b):return float(np.max(np.abs(a-b)/np.maximum(1,np.abs(b))))

manifest_path=BASE/"prepared/MANIFEST.json"
campaign_path=BASE/"campaign/CAMPAIGN.json"
convergence_path=BASE/"matching_validation/convergence/CONVERGENCE.json"
envelope_path=BASE/"matching_validation/ENVELOPE_AUDIT.json"
manifest,campaign,convergence,envelope=map(read,[manifest_path,campaign_path,convergence_path,envelope_path])
records={r["id"]:r for r in manifest["source_records"]}
assert len(records)==34
assert sha(manifest_path)==campaign["prepared_manifest_sha256"]==envelope["manifest_sha256"]
assert sha(path("python_frontend/compose_eos/new_model_matching.py"))==campaign["adapter_sha256"]
assert sha(BASE/"matching_validation/audit_success_envelopes.py")==envelope["script_sha256"]
assert campaign["git_base_head"]=="00bb5e1adcdbb8dfa0a0b2d5d7dd6f9037a4ddc3"

prepared={};source_rows=[];public_rows=[];checked_hashes={}
def check_hash(p,expected):
    p=path(p);observed=sha(p)
    assert observed==expected,(str(p),expected,observed)
    checked_hashes[str(p)]=observed
for rec in records.values():
    for name,h in rec["source_table_sha256"].items():check_hash(path(rec["source_path"])/name,h)
    for name,h in rec["internal_sha256"].items():check_hash(path(rec["internal_path"])/name,h)
    check_hash(rec["diagnostic_path"],rec["diagnostic_sha256"])
    d,m,meta=intern(path(rec["internal_path"]));prepared[rec["id"]]=(d,m,meta)
    diag=read(path(rec["diagnostic_path"]))
    intervals=diag["intervals"]
    assert len(intervals)==len(d)-1
    valid=np.array([v["valid"] for v in intervals],dtype=bool)
    components=[];start=None
    for i,v in enumerate(valid):
        if v and start is None:start=i
        if not v and start is not None:
            if i+1-start>=3:components.append((start,i+1))
            start=None
    if start is not None and len(d)-start>=3:components.append((start,len(d)))
    assert components==[(c["start_zero_based"],c["stop_exclusive"]) for c in rec["components"]],rec["id"]
    for a,b in components:assert len(np.unique(d["Iphase"][a:b]))==1,(rec["id"],a,b)
    source_rows.append({"id":rec["id"],"rows":len(d),"connected_components":len(components),
                        "interval_rejections":dict(Counter(i["category"] for i in intervals if not i["valid"]))})
    if rec["origin"]=="new_public":
        sel=rec["selection"];sidepath=path(sel["auxiliary_sidecar"])
        check_hash(sidepath,sel["auxiliary_sha256"]);side=read(sidepath)
        original=path(rec["source_path"])
        rawlines=(original/"eos.thermo").read_text().splitlines()[1:]
        raw=[s.split() for s in rawlines if s.strip()]
        comp=[s.split() for s in (original/"eos.compo").read_text().splitlines() if s.strip()]
        assert len(side["rows"])==len(raw)==len(rawlines)
        for i,(s,r,line,c) in enumerate(zip(side["rows"],raw,rawlines,comp)):
            assert s["q_tokens"]==r[3:10] and s["auxiliary_tokens"]==r[11:]
            assert s["raw_thermo_line"]==line and s["source_iphase"]==int(c[3])
            assert int(r[10])==sel["nadd_original"] and len(r[11:])==sel["nadd_original"]
        a,b=sel["source_start_zero_based"],sel["source_stop_exclusive"]
        for j in range(7):
            q=np.array([num(r[3+j]) for r in raw[a:b]])
            assert np.array_equal(d["Q"+str(j+1)],q),(rec["id"],j)
        selected=[int(c[3]) for c in comp[a:b]]
        assert all((p==0)==(rec["role"]=="core") for p in selected)
        public_rows.append({"id":rec["id"],"Nadd":sel["nadd_original"],"selected_rows":[a,b],
                            "q1_q7_exact":True,"all_auxiliary_tokens_and_lines_exact":True,
                            "original_phase_values":sorted(set(selected))})

cells=campaign["cells"]
assert len(cells)==campaign["completed_cells"]==campaign["requested_cells"]==120
assert len({(c["category"],c["id"]) for c in cells})==120
counts=Counter(c["category"] for c in cells)
assert counts=={"old_native_crust__new_native_core":65,"new_public_crust__old_public_GM1":5,
               "new_public_crust__new_public_core":20,"reference_assisted_hybrid":30}
summary={}
special=[]
for cat in counts:
    group=[c for c in cells if c["category"]==cat]
    summary[cat]={"model_variant_pairs":len(group),"pairs_with_output":sum(c["generated_count"]>0 for c in group),
      "generated_tables":sum(c["generated_count"] for c in group),
      "actual_kernel_invocations":sum(c["kernel_invocations"] for c in group),
      "candidate_attempts":sum(len(c["candidate_results"]) for c in group),
      "pair_statuses":dict(Counter(c["status"] for c in group)),
      "candidate_rejections":dict(Counter(r["category"] for c in group for r in c["candidate_results"] if r["status"]=="REJECTED"))}
assert summary==campaign["summary"]
for c in cells:
    assert c["prepared_manifest_sha256"]==sha(manifest_path)
    assert c["kernel_invocations"]==sum(r["kernel_invoked"] for r in c["candidate_results"])
    assert c["generated_count"]==sum(r["status"]=="GENERATED_VERIFIED_FINITE_DOMAIN" for r in c["candidate_results"])
    for i,r in enumerate(c["candidate_results"]):
        if r.get("category") in ["UNKNOWN","MULTIPLE_COEXISTENCE_CROSSINGS"]:
            special.append({"category":c["category"],"cell":c["id"],"candidate":i,"rejection":r["category"],"error":r["error"]})
assert sum(s["actual_kernel_invocations"] for s in summary.values())==4197
observables=["mu_transition","pressure_transition","n_low","n_high","epsilon_low","epsilon_high","density_jump","latent_energy_density"]

def audit_output(r,cell,coarse=False):
    out=path(r["output_path"]);home=out.parent
    assert read(home/"RESULT.json")==r,(cell["id"],coarse,"embedded result differs")
    for name,h in r["output_sha256"].items():check_hash(out/name,h)
    d,m,meta=intern(home/"readback_internal")
    expected,em,emeta=intern(home/"expected_internal")
    assert m==em and meta==emeta
    assert all(np.array_equal(d[k],expected[k]) for k in d.dtype.names)
    ntok=(out/"eos.nb").read_text().split()
    n=np.array([num(s) for s in ntok[2:]])
    lines=(out/"eos.thermo").read_text().splitlines()
    header=[num(s) for s in lines[0].split()]
    t=np.array([[num(v) for v in line.split()] for line in lines[1:] if line.strip()])
    assert header[0]==m and np.array_equal(n,d["nb"]) and len(t)==len(n)
    for j in range(7):assert np.array_equal(t[:,3+j],d["Q"+str(j+1)]),(cell["id"],coarse,"raw_Q",j)
    n,p,e,muq=physical(d,m);mue=(e+p)/n
    a,b=r["v2_report"]["transition_row_indices_zero_based"]
    assert b==a+1 and np.all(d["Iphase"][:b]==1) and np.all(d["Iphase"][b:]==2)
    assert n[b]>n[a] and e[b]>e[a] and np.all(np.diff(n)>0)
    pressure_relative=abs(p[b]-p[a])/max(abs(p[a]),1e-300)
    mu_relative=abs(mue[b]-mue[a])/max(abs(mue[a]),1e-300)
    latent_relative=abs(e[b]-e[a]-mue[a]*(n[b]-n[a]))/max(abs(e[b]-e[a]),1e-300)
    assert max(pressure_relative,mu_relative,latent_relative)<1e-10
    rawvalues={"mu_transition":mue[a],"pressure_transition":p[a],"n_low":n[a],"n_high":n[b],
               "epsilon_low":e[a],"epsilon_high":e[b],"density_jump":n[b]-n[a],"latent_energy_density":e[b]-e[a]}
    transition_max=max(abs(rawvalues[k]-r["transition"][k])/max(1,abs(r["transition"][k])) for k in observables)
    assert transition_max<5e-13
    retained={}
    for side,selector,sourceid in [("low",slice(0,a),cell["crust_source_id"]),("high",slice(b+1,len(n)),cell["core_source_id"])]:
        if coarse and side=="high":
            sourceid=sourceid.removesuffix(":fine")+":coarse"
        source,sm,_=prepared[sourceid];sn,sp,se,smu=physical(source,sm)
        ids=np.searchsorted(sn,n[selector])
        assert np.all(ids<len(sn))
        errs={"n":scaled(n[selector],sn[ids]),"pressure":scaled(p[selector],sp[ids]),"epsilon":scaled(e[selector],se[ids]),"mu_from_Q3":scaled(muq[selector],smu[ids])}
        assert max(errs.values())<5e-13,(cell["id"],coarse,side,errs)
        retained[side]=errs
    assert r["admission"]["status"]=="ADMITTED" and not r["admission"]["failed_checks"]
    assert r["serialized_native_verification"]["overall"]=="PASS"
    assert r["readback"]["numerical"]["status"]==r["readback"]["metadata"]["status"]=="PASS"
    assert r["v2_report"]["native_report"]["physical_energy_shift_mev"]==0
    assert r["v2_report"]["native_report"]["composition_scope"]=="thermo_only_phase_labels_no_abundances"
    return {"cell":cell["id"],"category":cell["category"],"level":"coarse" if coarse else "fine",
      "output_path":r["output_path"],"rows":len(n),"density_domain":[float(n[0]),float(n[-1])],
      "transition":r["transition"],"interface_rows":[a,b],
      "raw_interface_relative_residuals":{"P":float(pressure_relative),"mu_E":float(mu_relative),"latent":float(latent_relative)},
      "raw_transition_max_scaled":float(transition_max),"retained_physical_rows_max_scaled":retained,
      "raw_internal_exact":True,"output_hashes_match":True,"admission":"ADMITTED",
      "energy_shift_mev":0,"composition":"thermo_only; source compositions retained outside matched artifact"}

outputs=[];success_map={}
for c in cells:
    for i,r in enumerate(c["candidate_results"]):
        if r["status"]!="GENERATED_VERIFIED_FINITE_DOMAIN":continue
        check=audit_output(r,c);check["candidate_index"]=i;outputs.append(check)
        success_map[(c["category"],c["id"],i)]=(c,r)
assert len(outputs)==9
conv=[];coarse_outputs=[]
for row in convergence["results"]:
    key=(row["category"],row["cell"],row["fine_candidate_index"]);c,fine=success_map[key]
    assert row["status"]=="SAME_COMPONENT_TRANSITION_COMPARED"
    coarse=row["coarse_result"]
    coarse_outputs.append(audit_output(coarse,c,True))
    assert fine["transition"]==row["fine_transition"] and coarse["transition"]==row["coarse_transition"]
    absolute={k:abs(fine["transition"][k]-coarse["transition"][k]) for k in observables}
    relative={k:absolute[k]/abs(fine["transition"][k]) for k in observables}
    assert all(v==row["absolute_difference"][k] for k,v in absolute.items())
    assert all(v==row["relative_difference"][k] for k,v in relative.items())
    conv.append({"category":row["category"],"cell":row["cell"],"absolute_difference":absolute,"relative_difference":relative,
                 "maxwell_identity_absolute_difference":abs(fine["transition"]["maxwell_identity_residual"]-coarse["transition"]["maxwell_identity_residual"])})
assert len(conv)==7
assert envelope["status"]=="PASS" and len(envelope["results"])==9
for r in envelope["results"]:
    assert (r["category"],r["cell"],r["candidate_index"]) in success_map
    assert r["root_count"]==1 and r["deltaP_low_endpoint"]<0 and r["deltaP_high_endpoint"]>0 and r["slope_deltaP_at_transition_fm3"]>0
alloutputs=outputs+coarse_outputs
result={"schema":"ns_cowling.independent_matching_readback.v1","status":"IMPLEMENTATION_VERIFIED_FINITE_DOMAIN",
    "method":"One selective read of final summaries, recorded raw candidate results, prepared source/diagnostic identities, exact Q1-Q7/auxiliary sources, nine primary and seven coarse actual output tables. No campaign/native rerun.",
    "identities":{str(p):sha(p) for p in [manifest_path,campaign_path,convergence_path,envelope_path]},
    "source_records_checked":source_rows,"public_source_transport_checked":public_rows,
    "source_file_identities_checked":checked_hashes,"recomputed_summary":summary,
    "fine_outputs":outputs,"coarse_outputs":coarse_outputs,"convergence":conv,
    "max_relative_fine_coarse":{k:max(x["relative_difference"][k] for x in conv) for k in observables},
    "max_raw_interface_relative_residual":max(max(o["raw_interface_relative_residuals"].values()) for o in alloutputs),
    "envelope_evidence":{"status":envelope["status"],"outputs":len(envelope["results"]),"max_root_difference_MeV":max(r["root_difference_from_native_MeV"] for r in envelope["results"]),"scope":"Numerical polynomial root enumeration of the same finite Hermite component supports; not a global microscopic or continuum proof"},
    "special_rejections":special,
    "remaining_scope":["Public crust precision is frozen; no independent continuous TF mesh or source print precision convergence is inferred.",
      "Matched artifacts are thermo-only and have finite inner-crust/core support; no whole-star, TOV, mode or microscopic global-minimum certification.",
      "UNKNOWN rejection messages remain separate: three insufficient retained-knot failures and one C2 verification failure, not absence-of-root theorems."]}
(OUT/"MATCHING_READBACK.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
print(json.dumps({"status":result["status"],"source_records":len(source_rows),"public_transports":len(public_rows),
 "primary_outputs":len(outputs),"coarse_outputs":len(coarse_outputs),"max_interface_relative":result["max_raw_interface_relative_residual"],
 "max_relative_fine_coarse":result["max_relative_fine_coarse"],"envelope":result["envelope_evidence"]},indent=2))
