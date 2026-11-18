"""Immutable-H endpoint-matched EOS interpolation crosscheck.

No refitting of physical Maxwell endpoints. Alternative lower bridges are
evaluation-only objects, not production EOS contracts or regenerated TF states.
"""
from pathlib import Path
import csv, hashlib, json, os, sys, time, traceback
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
sys.path.insert(0,str(ROOT/"scripts"))
from matching_interpolation_shapes_20261118 import CausalShapeInterval, naive_loglog_polytrope_endpoint
from compose_eos.whole_crust import WholeCrustBarotrope, CausalExponentialInterval, HermiteInterval
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS, solve_sharp_background

DATA=ROOT/"data/matching_interpolation_20261118"
H=ROOT/"data/low_density_crust_20261117/stellar/causal_extension"
PLAN=H/"TOV_PLAN.json"
H_REF=H/"TOV_RESULTS.json"
RUNNER=Path(__file__)
BASE_SHA="634a586fdd0c214bf7fa9718b02ebafa3e927008"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def native_json_value(value):
    if isinstance(value,np.generic):
        return value.item()
    if isinstance(value,np.ndarray):
        return value.tolist()
    raise TypeError(f"unsupported non-JSON diagnostic type: {type(value)!r}")
def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False,
                               default=native_json_value)+"\n")
def summary(bg):
    x=bg.summary()
    return dict(mass_msun=x["mass_msun"],radius_km=x["radius_km"],
        central_density_fm3=x["central_density_fm3"],compactness=x["compactness"],
        interfaces=x["interfaces"],nfev=x["nfev"])
def endpoint_args(c):
    l,h=c["lower_bridge"]["endpoint_states"]
    return (*[l[k] for k in ("n","p","epsilon")],
            *[h[k] for k in ("n","p","epsilon")])
def tsolve(eos,n,setting):
    return solve_sharp_background(eos,n,**setting)
def star_case(bar,reference,settings,tighter=None):
    eos=SharpContractEnthalpyEOS(bar)
    n_ref=reference["background"]["central_density_fm3"]
    fixed=tsolve(eos,n_ref,settings)
    lo,hi=0.16,0.32
    mlo,mhi=tsolve(eos,lo,settings).mass_msun,tsolve(eos,hi,settings).mass_msun
    if not mlo<1.4<mhi:
        raise ValueError(f"precommitted mass bracket lost: {mlo},{mhi}")
    nc=brentq(lambda n:tsolve(eos,n,settings).mass_msun-1.4,lo,hi,xtol=1e-10)
    target=tsolve(eos,nc,settings)
    control=None
    if tighter is not None:
        ct=tsolve(eos,nc,tighter)
        control=dict(same_central_density=True,delta_R_km=ct.radius_km-target.radius_km,
                     delta_M_msun=ct.mass_msun-target.mass_msun)
    return dict(fixed_central=summary(fixed),target_mass=summary(target),
                tighter_control=control)
def physics_samples(bridge,ref):
    x=np.linspace(0,1,301)
    mu=bridge.mu_left+bridge.width*x
    s=bridge.evaluate_mu(mu)
    c2=np.asarray(s["cs2_over_c2"])
    p=np.asarray(s["p"])
    epsilon=np.asarray(s["epsilon"])
    n=np.asarray(s["n"])
    p0=np.asarray(ref["p"]);c0=np.asarray(ref["cs2_over_c2"])
    independent_dp=quad(lambda mm:float(bridge.evaluate_mu(mm)["n"]),bridge.mu_left,
                        bridge.mu_right,epsabs=1e-13,epsrel=3e-12,limit=200)[0]
    return dict(max_abs_delta_P_MeV_fm3=float(np.max(abs(p-p0))),
        max_abs_delta_epsilon_MeV_fm3=float(np.max(abs(epsilon-ref["epsilon"]))),
        max_abs_delta_n_fm3=float(np.max(abs(n-ref["n"]))),
        max_abs_delta_cs2_over_c2=float(np.max(abs(c2-c0))),
        min_cs2_over_c2=float(np.min(c2)),max_cs2_over_c2=float(np.max(c2)),
        cs2_left=float(c2[0]),cs2_right=float(c2[-1]),
        independent_pressure_area_error_MeV_fm3=float(independent_dp-(bridge.p_right-bridge.p_left)),
        independent_positive_n=bool(np.all(np.diff(n)>0)),
        phase_boundary_values_modified=False)
def main():
    t0=time.time();data=json.loads(PLAN.read_text());ref=json.loads(H_REF.read_text())
    by_id={x["id"]:x for x in ref["cases"]}
    assert len(data["cases"])==8 and len(by_id)==8
    source_paths=[PLAN,H_REF,RUNNER,ROOT/"scripts/matching_interpolation_shapes_20261118.py",
                  ROOT/"python_frontend/compose_eos/whole_crust.py",
                  ROOT/"python_frontend/compose_eos/stellar_sharp.py"]
    source_paths += [ROOT/x["contract"] for x in data["cases"]]
    hashes_before={str(p.relative_to(ROOT)):sha(p) for p in source_paths}
    record=dict(stage="matching-interpolation-crosscheck",base_commit=BASE_SHA,
        scope="Only the external-to-TF lower bridge changes; sharp Maxwell doublets and TF/core tables fixed",
        source_hashes_before=hashes_before,source_contract_count=8,
        methods=["Hermite P(mu) cubic","Exponential monotone n/mu","Power monotone n/mu",
                 "Rational monotone n/mu","Naive log-log P(n) endpoint negative control"],
        cases=[],status="RUNNING")
    save(DATA/"RESULTS.json",record)
    setting=data["settings"];tighter=data["tighter_settings"]
    csv_rows=[]
    for i,case in enumerate(data["cases"]):
        assert sha(ROOT/case["contract"])==case["contract_sha256"]
        c=json.loads((ROOT/case["contract"]).read_text())
        args=endpoint_args(c)
        baseline=CausalExponentialInterval(*args)
        canonical_bar=WholeCrustBarotrope.from_contract(c)
        expected=by_id[case["id"]]["background"]
        sref=baseline.evaluate_mu(baseline.mu_left+baseline.width*np.linspace(0,1,301))
        caseout=dict(id=case["id"],crust=case["crust"],core=case["core"],
             join_n_fm3=case["upper_start_density"],label=case["label"],
             known_native_failures_in_bridge=case["known_failed_native_knots_in_closure"],
             contract_sha256=case["contract_sha256"],Maxwell_source_preserved=dict(
                  sharp_indices=sorted(canonical_bar.sharp_interfaces),
                  pressure_MeV_fm3=float(canonical_bar.upper_p[min(canonical_bar.sharp_interfaces)]),
                  n_low_fm3=float(canonical_bar.upper_n[min(canonical_bar.sharp_interfaces)]),
                  n_high_fm3=float(canonical_bar.upper_n[min(canonical_bar.sharp_interfaces)+1])),
             methods={})
        methods=caseout["methods"]
        hermite=HermiteInterval(*args)
        hv=hermite.validate()
        methods["hermite"]=dict(status=hv["status"],validation=hv,
            expected_admission="PASS iff stable causal cubic",star=None)
        poly=naive_loglog_polytrope_endpoint(*args)
        methods["loglog_polytrope"]=dict(status="FAIL_ENDPOINT_ENERGY" if not poly["thermodynamically_valid_for_fixed_endpoints"] else "PASS",diagnostic=poly,star=None)
        family=[("exponential",baseline),
                ("power",CausalShapeInterval(*args,family="power")),
                ("rational",CausalShapeInterval(*args,family="rational"))]
        if hv["status"]=="PASS":
            family.insert(1,("hermite",hermite))
        for name,bridge in family:
            val=bridge.validate()
            if val["status"]!="PASS":
                methods[name]=dict(status="FAIL",validation=val,star=None)
                continue
            bar=WholeCrustBarotrope.from_contract(c)
            bar.bridge=bridge
            st=star_case(bar,by_id[case["id"]],setting,
                        tighter if case["label"]=="primary" and name in ("power","rational") else None)
            ss=physics_samples(bridge,sref)
            r=st["target_mass"]["radius_km"]
            if name=="exponential":
                assert abs(r-expected["radius_km"])<3e-5,("H baseline drift",case["id"],r,expected["radius_km"])
            sharp=st["target_mass"]["interfaces"][0]
            assert abs(sharp["p_MeV_fm3"]-caseout["Maxwell_source_preserved"]["pressure_MeV_fm3"])<1e-11
            methods[name]=dict(status="PASS_SCOPED",validation=val,EOS=ss,star=st,
                delta_R_from_H_exponential_km=float(r-expected["radius_km"]),
                min_max_cs2_over_c2=[ss["min_cs2_over_c2"],ss["max_cs2_over_c2"]])
            csv_rows.append(dict(case=case["id"],method=name,status="PASS",
                join_n_fm3=case["upper_start_density"],
                R_1p4_km=r,delta_R_from_H_km=r-expected["radius_km"],
                min_cs2=ss["min_cs2_over_c2"],max_cs2=ss["max_cs2_over_c2"],
                max_delta_cs2_vs_exp=ss["max_abs_delta_cs2_over_c2"]))
            print("TOV",case["id"],name,"R",f"{r:.10f}","cs2_max",f"{ss['max_cs2_over_c2']:.5f}",flush=True)
        caseout["H_reference_radius_km"]=expected["radius_km"]
        record["cases"].append(caseout);save(DATA/"RESULTS.json",record)
        print("CASE",i+1,"/",len(data["cases"]),case["id"],"METHODS",[(k,v["status"]) for k,v in methods.items()],flush=True)
    record["source_hashes_after"]={str(p.relative_to(ROOT)):sha(p) for p in source_paths}
    record["source_bytes_unchanged"]=hashes_before==record["source_hashes_after"]
    record["status"]="PASS_SCOPED" if record["source_bytes_unchanged"] and all(
        c["methods"][m]["status"]=="PASS_SCOPED" for c in record["cases"]
        for m in ("exponential","power","rational")) else "FAIL_CHECK_RECEIPT"
    record["duration_seconds"]=time.time()-t0
    save(DATA/"RESULTS.json",record)
    with (DATA/"METHOD_COMPARISON.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(csv_rows[0]));writer.writeheader();writer.writerows(csv_rows)
    print("FINISHED",record["status"],"seconds",record["duration_seconds"],flush=True)
    return 0 if record["status"]=="PASS_SCOPED" else 1
if __name__=="__main__":
    raise SystemExit(main())
