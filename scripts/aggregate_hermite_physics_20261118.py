"""Direct physical relative-error and matching/Cowling aggregation for cubic/quintic.

All percentages are defined against held-out *same-model* dense native states
or, for phase transitions, against the actual full-source native Maxwell
result. CSS is never described as an error relative to the Maxwell mechanism.
Frequencies use paired same-mass, same-closure settings; shifts are between
methods, NOT errors against true stellar observations.
"""
from pathlib import Path
import sys,json,math,csv
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from compare_hermite_orders_20261118 import predict
from benchmark_native_dense_eos_20261118 import load_native,phase_runs,native_truth_cs2
from benchmark_tf_native_holdout_20261118 import collect
OUT=ROOT/"data/cubic_quintic_matching_20261118"
DENSE=OUT/"dense_cv_v2"
TF=OUT/"tf_cv"
MAX=OUT/"maxwell"
MATCH=OUT/"native_matching"
HST=OUT/"stellar_bounded"
NAT=OUT/"native_stellar"
RESCUE=OUT/"p1_node_refinement"
TIGHT=OUT/"p1_tight_controls"
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def gather_tf(case):
    original=json.loads((ROOT/f"data/native_dense_interpolation_20261118/tf_holdout/{case}.json").read_text())
    anchors,vals,base,recovered=collect(case)
    ah={x["_source_sha256"]:x for x in anchors.values()}
    vh={x["sha"]:x["state"] for x in vals}
    rows=[]
    for key,digest in original["training_shas"].items():
        v=ah[digest];rows.append((float(key),v,"train"))
    for key,digest in original["new_direct_shas"].items():
        v=vh[digest];rows.append((float(key),v,"holdout"))
    rows.sort(key=lambda x:x[0])
    src={k:np.array([float(v[k2]) if k2!="n" else n for n,v,label in rows],float)
         for k,k2 in (("n","n"),("p","pressure"),("mu","mu_B"),("epsilon","epsilon"))}
    tr=np.array([i for i,(_,v,label) in enumerate(rows) if label=="train"],int)
    hold=np.array([i for i,(_,v,label) in enumerate(rows) if label=="holdout"],int)
    return src,tr,hold
def metrics(src,tr,hold,method):
    z=predict(method,src,tr,hold)
    original={k:np.asarray(src[k])[hold] for k in ("n","p","epsilon")}
    reference=native_truth_cs2(src)["cs2_from_P_over_e"][hold]
    c=np.asarray(z["cs2"],float)
    native_ok=np.isfinite(reference)&(reference>1e-5)&(reference<=1)
    sampled_ok=np.isfinite(c)&(c>0)&(c<=1+1e-12)
    positive=(np.asarray(z["n"])>0)&(np.asarray(z["p"])>0)&(np.asarray(z["epsilon"])>0)
    def rms_rel(field):
        pred=np.asarray(z[field])
        return float(100*np.sqrt(np.mean(((pred-original[field])/original[field])**2)))
    cr=100*np.sqrt(np.mean(((c[native_ok&sampled_ok]-reference[native_ok&sampled_ok])/reference[native_ok&sampled_ok])**2)) if np.any(native_ok&sampled_ok) else None
    return dict(n_percent=rms_rel("n"),P_percent=rms_rel("p"),epsilon_percent=rms_rel("epsilon"),
                cs2_relative_percent=float(cr) if cr is not None else None,
                cs2_absolute=float(np.sqrt(np.mean((c[native_ok&sampled_ok]-reference[native_ok&sampled_ok])**2))) if np.any(native_ok&sampled_ok) else None,
                physical_fraction=float(np.sum(native_ok&sampled_ok&positive)/max(1,np.sum(native_ok))),
                native_causal_count=int(np.sum(native_ok)),holdout_count=len(hold),
                numeric_identity_max_abs=float(np.max(abs(z["epsilon"]+z["p"]-z["mu"]*z["n"]))),
                min_predicted_cs2=float(np.nanmin(c)),max_predicted_cs2=float(np.nanmax(c)))
def collect_dense():
    rows=[]
    for family in ("rmf","core"):
        for modeldir in sorted((DENSE/family).iterdir()):
            if not modeldir.is_dir():continue
            model=modeldir.name
            evidence=json.loads((modeldir/"RESULT.json").read_text())
            src=load_native(family,model)
            for record in evidence["runs"]:
                if "methods" not in record:continue
                a,b=record["native_row_slice"]
                sub={k:(v[a:b] if isinstance(v,np.ndarray) else v) for k,v in src.items()}
                tr=np.array(record["training_indices"],int)
                hold=np.setdiff1d(np.arange(len(sub["mu"])),tr,assume_unique=True)
                assert len(hold)==record["holdout_rows"]
                for method in ("cubic","quintic"):
                    data=metrics(sub,tr,hold,method)
                    rows.append(dict(family=family,model=model,phase=record["phase"],
                                     stride=record["training_stride"],method=method,**data))
    for model in ("DDME2","DDMEX","NL3","NL3wrL55","TM1e"):
        src,tr,hold=gather_tf(model)
        for method in ("cubic","quintic"):
            d=metrics(src,tr,hold,method)
            rows.append(dict(family="tf",model=model,phase="N128 droplet",stride="H_fixed_knots",
                             method=method,**d))
    return rows
def aggregate(records):
    # One explicit training-grid level at a time. No pooling of the same native
    # model over different withheld training sparsities as independent data.
    from collections import defaultdict
    groups=defaultdict(list)
    for row in records:groups[(row["family"],str(row["stride"]),row["method"])].append(row)
    result={}
    for (family,stride,method),items in sorted(groups.items()):
        total=sum(v["holdout_count"] for v in items)
        val=lambda k:float(np.sqrt(sum(v["holdout_count"]*v[k]**2 for v in items)/total))
        denom=sum(x["native_causal_count"] for x in items)
        result.setdefault(family,{}).setdefault(stride,{})[method]=dict(
          n_percent=val("n_percent"),P_percent=val("P_percent"),epsilon_percent=val("epsilon_percent"),
          physical_fraction=sum(x["native_causal_count"]*x["physical_fraction"] for x in items)/max(1,denom),
          cs2_relative_percent=(float(np.sqrt(sum(x["native_causal_count"]*(x["cs2_relative_percent"] or 0)**2 for x in items if x["cs2_relative_percent"] is not None)/max(1,sum(x["native_causal_count"] for x in items if x["cs2_relative_percent"] is not None))))),
          holdout_count=total,native_causal_count=denom,
          number_source_models=len(set(x["model"] for x in items)))
    return result
def collect_maxwell():
    rows=[]
    for p in sorted(MAX.glob("*/RESULT.json")):
        x=json.loads(p.read_text());native=x["native_full_grid_maxwell"];pair=p.parent.name
        for level in ("sparse_4x","sparse_8x"):
            for method in ("cubic","quintic"):
                q=x["fits"][f"{method}_{level}"]
                ref=q.get("comparison_native_maxwell",{})
                rows.append(dict(pair=pair,level=level,method=method,native_status=native["status"],
                     n_roots=q.get("crossings"),approved=q.get("globally_accepted_single_root",False),
                     mu_percent=100*ref["mu_relative"] if ref else None,
                     Pt_percent=100*ref["pressure_relative"] if ref else None,
                     delta_n_percent=100*ref["density_jump_relative"] if ref else None,
                     latent_heat_percent=100*ref["latent_heat_relative"] if ref else None,
                     native_reference=native))
    return rows
def collect_matching():
    rows=[]
    for p in sorted(MATCH.glob("*/RESULT.json")):
        x=json.loads(p.read_text())
        for m,z in x["methods"].items():
            q=z.get("transition_or_bridge",{})
            rows.append(dict(pair=p.parent.name,method=m,status=z["status"],
              bridge_n1=q.get("n1"),bridge_n2=q.get("n2"),bridge_P1=q.get("P1"),bridge_P2=q.get("P2"),
              min_bridge_cs2=q.get("cs2_min_bridge"),max_bridge_cs2=q.get("cs2_max_bridge"),
              css_sharp=bool(q.get("is_exact_maxwell")),css_transition_pressure=q.get("css_transition_pressure"),
              css_delta_eps=q.get("css_transition_delta_eps"),css_jump_density=q.get("css_transition_nb_high",0)-q.get("css_transition_nb_low",0)))
    return rows
def selected_mode(scope,case,method,conv,family):
    # Stage-1 mode record, then separately admitted node-sorted p1 rescue.
    try:
        if scope=="H":
            orig=json.loads((HST/case/"RESULT.json").read_text())["cases"][method]
        else:
            orig=json.loads((NAT/case/"RESULT.json").read_text())["method_results"][method]
        d=orig.get("modes",{}).get(conv+"_"+family)
        if d and d.get("status") in ("PASS","PASS_SCOPED") and d.get("frequency_hz"):
            return dict(frequency_hz=float(d["frequency_hz"]),status="PASS",origin="original_mode_run",
                        numerical_hz=d.get("frequency_control_delta_hz"))
        if family=="p1":
            p=RESCUE/scope/case/method/conv/"RESULT.json"
            if p.exists():
                x=json.loads(p.read_text())
                if x["status"]=="PASS_SCOPED_NODE_CLASSIFIED_P1":
                    q=x["selected"]
                    return dict(frequency_hz=float(q["tighter_frequency_hz"]),status="PASS",
                         origin="independent_node_classification",numerical_hz=q.get("delta_numerical_hz"))
            p=TIGHT/scope/case/method/conv/"RESULT.json"
            if p.exists():
                x=json.loads(p.read_text())
                if x["status"]=="PASS_SCOPED_NODE_CLASSIFIED_P1":
                    q=x["selected"]
                    return dict(frequency_hz=float(q["frequency_hz"]),status="PASS",
                        origin="tightened_node_classification",numerical_hz=q["frequency_hz"]-q["earlier_frequency"])
    except Exception as exc:
        return dict(status="BLOCKED",reason=str(exc))
    return dict(status="BLOCKED",reason="no residual-and-node-admitted mode of expected family")
def collect_stellar():
    rows=[];pairs=[]
    for p in sorted(HST.glob("*/RESULT.json")):
        x=json.loads(p.read_text());case=x["case"]
        for method in ("cubic","quintic"):
            m=x["cases"].get(method)
            if not m or "background" not in m:continue
            for conv in ("slow","rapid"):
                z=dict(family="H_fixed_maxwell",case=case,method=method,conversion=conv,
                     mass_msun=m["background"]["mass_msun"],R_km=m["background"]["radius_km"],
                     central_n=m["background"]["central_density_fm3"],
                     mode_f=selected_mode("H",case,method,conv,"f"),
                     mode_p1=selected_mode("H",case,method,conv,"p1"))
                rows.append(z)
    for p in sorted(NAT.glob("*/RESULT.json")):
        x=json.loads(p.read_text());case=x["pair"]
        for method,m in x["method_results"].items():
            if "background" not in m:continue
            for conv in (("slow","rapid") if method in ("css","native_maxwell") else ("slow",)):
                z=dict(family="native_matching",case=case,method=method,conversion=conv,
                     mass_msun=m["background"]["mass_msun"],R_km=m["background"]["radius_km"],
                     central_n=m["background"]["central_density_fm3"],
                     mode_f=selected_mode("native",case,method,conv,"f"),
                     mode_p1=selected_mode("native",case,method,conv,"p1"))
                rows.append(z)
    for scope,orders in (("H_fixed_maxwell",("cubic","quintic")),("native_matching",("cubic_P","quintic_P"))):
        cases=sorted(set(q["case"] for q in rows if q["family"]==scope))
        for c in cases:
            for conv in ("slow","rapid") if scope=="H_fixed_maxwell" else ("slow",):
                matched=[q for q in rows if q["family"]==scope and q["case"]==c and q["conversion"]==conv]
                got={q["method"]:q for q in matched}
                if not all(o in got for o in orders):continue
                a,b=(got[k] for k in orders)
                rec=dict(family=scope,case=c,conversion=conv,baseline=orders[0],alternative=orders[1],
                    R_relative_percent=100*(b["R_km"]-a["R_km"])/a["R_km"],
                    center_n_relative_percent=100*(b["central_n"]-a["central_n"])/a["central_n"])
                for label in ("f","p1"):
                    aa=a["mode_"+label];bb=b["mode_"+label]
                    rec[label+"_relative_percent"]=(100*(bb["frequency_hz"]-aa["frequency_hz"])/aa["frequency_hz"]) if aa.get("status")=="PASS" and bb.get("status")=="PASS" else None
                pairs.append(rec)
    return rows,pairs
def main():
    native=collect_dense()
    if len(native)!=2*(5*2+4*2+3*2+5*1):
        # 5 RMF 2 strides, four core single-phase 2 strides, APR three phases 2 strides, TF five
        raise ValueError(f"unexpected model-phase count {len(native)}")
    aggregate_table=aggregate(native)
    maxwell=collect_maxwell()
    matching=collect_matching()
    stellar,comparisons=collect_stellar()
    out=dict(status="COMPUTED_WITH_PRESERVED_REJECTIONS",
        source_commit="2b32672da2b9b602942a362ff7208510708be459",
        scope="actual matched native sources, sparse vs dense cubic/quintic, fixed-mass TOV/Cowling, CSS distinct model",
        reference_density_models=10,tf_branch_models=5,native_direct_state_count=25612,tf_holdout_count=30,
        method_names=["cubic","quintic"],
        percent_semantics=dict(native="RMS signed model prediction errors relative to SAME native EOS at withheld mu; percent=100*(prediction-truth)/truth",
            matching="percent difference against full-source native Maxwell at same two sources; native Maxwell is an interpolated numerical construction, not continuum truth",
            mode="signed change of alternative relative to cubic at fixed gravitational mass; no experimental true frequency",
            css="CSS vs Maxwell shifts are DIFFERENT phase constructions, not method error"),
        family_levels=aggregate_table,dense_model_phase_runs=native,
        maxwell=maxwell,native_smooth_css=matching,stellar=stellar,stellar_paired=comparisons)
    save(OUT/"SUMMARY.json",out)
    def write_csv(name,records):
        if not records:return
        keys=list(records[0])
        with (OUT/name).open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=keys,extrasaction="ignore")
            w.writeheader();w.writerows(records)
    write_csv("NATIVE_ERRORS_BY_PHASE.csv",native)
    write_csv("MAXWELL_ROOT_ERRORS.csv",maxwell)
    write_csv("NATIVE_MATCHING_QUANTITIES.csv",matching)
    write_csv("STELLAR_MODES.csv",stellar)
    write_csv("STELLAR_PAIRED_ERRORS.csv",comparisons)
    print("AGGREGATED",len(native),len(maxwell),len(matching),len(stellar),len(comparisons))
    for k,x in aggregate_table.items():
        for level,metrics in x.items():
            a=metrics["cubic"];b=metrics["quintic"]
            print("FAMILY",k,level,"CUBIC_nRMS_PERCENT",a["n_percent"],
                  "QUINTIC_nRMS_PERCENT",b["n_percent"],
                  "CUBIC_cs2_PERCENT",a["cs2_relative_percent"],
                  "QUINTIC_cs2_PERCENT",b["cs2_relative_percent"])
if __name__=="__main__":main()
