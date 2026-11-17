#!/usr/bin/env python3
"""H opt-in causal lower closure: pre-observable source-knot admission."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,copy,hashlib,json,os,sys,time,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
BASE=ROOT/"data/low_density_crust_20261117/stellar/causal_extension"
INPUTS=ROOT/"data/low_density_crust_20261117/stellar/h_final/INPUTS.json"
REFERENCE=ROOT/"data/stellar_extension_20261117/eos/fine/REPRESENTATION.json"
GAPS=ROOT/"data/low_density_crust_20261117/low_extension/MANIFEST.json"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def rel(p):return str(Path(p).relative_to(ROOT))
def save(p,x):
 Path(p).parent.mkdir(parents=True,exist_ok=True)
 Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def stamp():return datetime.now(timezone.utc).isoformat()
def main(expected_whole_sha):
 whole=ROOT/"python_frontend/compose_eos/whole_crust.py"
 if sha(whole)!=expected_whole_sha:raise ValueError("whole_crust.py differs from explicitly frozen implementation")
 # Import only after checking the caller-provided frozen implementation identity.
 from compose_eos.api import read_compose
 from compose_eos.whole_crust import WholeCrustBarotrope,CausalExponentialInterval
 entries=load(INPUTS)
 if len(entries)>4:raise ValueError("bounded four-source inventory")
 core=[whole]+[ROOT/"python_frontend/compose_eos"/n for n in ("stellar_sharp.py","stellar_background.py","stellar_cowling.py")]
 files=[*core,Path(__file__),ROOT/"scripts/run_causal_low_density_backgrounds_20261117.py",INPUTS,REFERENCE,GAPS]
 for e in entries:
  cp=ROOT/e["compose_path"]
  files.extend(p for p in sorted(cp.glob("*")) if p.is_file())
  files.append(cp.parent/"RESULT.json")
 before={rel(p):sha(p) for p in files}
 plan={"schema":"ns_cowling.H_causal_closure_admission_plan.v1","created_at_utc":stamp(),"inputs":entries,
  "scope":"Opt-in phenomenological lower closure only; native sources and sharp Maxwell doublet unchanged.",
  "selection":"First and next causal-admissible actual source knot in ascending density, before stellar observables",
  "maximum_contracts":8,"kappa_selection":"Unique endpoint pressure integral; no observable or free parameter fit",
  "source_hashes_before":before,"expected_whole_crust_sha256":expected_whole_sha,
  "known_failed_native_density_knots":{m["model"]:m["failed_source_densities"] for m in load(GAPS)["models"]}}
 if (BASE/"ADMISSION_PLAN.json").exists():raise ValueError("Admission output already exists; immutable evidence will not be overwritten")
 save(BASE/"ADMISSION_PLAN.json",plan)
 receipt={"status":"RUNNING","started_at_utc":stamp(),"pid":os.getpid(),"actual_exit_code":None}
 save(BASE/"ADMISSION_PROCESS.json",receipt);start=time.monotonic()
 try:
  base=load(REFERENCE);bar=WholeCrustBarotrope.from_contract(base)
  left=dict(n=float(bar.outer.nb[-1]),p=float(bar.outer.p[-1]),epsilon=float(bar.outer.nb[-1]*bar.outer.F[-1]))
  records=[];cases=[]
  for entry in entries:
   cp=ROOT/entry["compose_path"];data=read_compose(cp);match=load(cp.parent/"RESULT.json")
   if match["status"]!="GENERATED_VERIFIED_FINITE_DOMAIN":raise ValueError("Unadmitted strict-match input")
   native=match.get("native_report") or match["v2_report"]["native_report"]
   interface=int(native["low_phase_rows"])-1
   record=dict(**entry,original_interface_index=interface,lower_join_trials=[],admitted=[])
   for j in range(interface):
    right=dict(n=float(data.nb[j]),p=float(data.p[j]),epsilon=float(data.eps[j]))
    trial=dict(source_row=j,**right,mu=(right["epsilon"]+right["p"])/right["n"])
    try:
     bridge=CausalExponentialInterval(left["n"],left["p"],left["epsilon"],right["n"],right["p"],right["epsilon"])
     trial["bridge_validation"]=bridge.validate()
     if trial["bridge_validation"]["status"]!="PASS":raise ValueError("causal connector finite validation rejected")
     full=dict(surface=copy.deepcopy(base["surface"]),outer=copy.deepcopy(base["outer"]),lower_bridge=bridge.contract(),
       upper=dict(n=data.nb[j:].tolist(),p=data.p[j:].tolist(),epsilon=data.eps[j:].tolist(),sharp_interfaces=[interface-j]))
     gaps=[n for n in plan["known_failed_native_density_knots"].get(entry["crust"],[]) if left["n"]<n<right["n"]]
     full["provenance"]=dict(source_compose=entry["compose_path"],source_result=rel(cp.parent/"RESULT.json"),
       source_row_start=j,source_files_sha256={p.name:sha(p) for p in sorted(cp.glob("*")) if p.is_file()},
       outer_reference=rel(REFERENCE),outer_sha256=sha(REFERENCE),
       choice=plan["selection"],source_values_preserved=True,physical_energy_shift_MeV=0.,
       lower_bridge_is_explicit_phenomenological_closure=True,
       closure_density_domain_fm3=[left["n"],right["n"]],known_failed_native_knots_in_closure=gaps,
       closure_does_not_repair_native_failures=True,disconnected_low_states_imported=False,
       microscopic_support_claim_for_bridge=False,geometry_global_minimum_proven=False)
     checked=WholeCrustBarotrope.from_contract(full)
     for key,arr in (("n",data.nb[j:]),("p",data.p[j:]),("epsilon",data.eps[j:])):
      assert np.array_equal(np.asarray(full["upper"][key]),arr)
     label="primary" if not record["admitted"] else "next_admissible_join_control"
     dest=BASE/"contracts"/entry["id"]/label/"REPRESENTATION.json"
     save(dest,full);restored=load(dest);WholeCrustBarotrope.from_contract(restored)
     for key in ("surface","outer"):assert restored[key]==base[key]
     trial["status"]="ADMITTED"
     case=dict(id=entry["id"]+"__"+label,source_id=entry["id"],crust=entry["crust"],core=entry["core"],label=label,
       source_row_start=j,retained_crust_rows=interface-j+1,upper_start_density=right["n"],
       bridge_density_width=right["n"]-left["n"],kappa=bridge.kappa,
       min_cs2_over_c2=trial["bridge_validation"]["min_cs2_over_c2"],max_cs2_over_c2=trial["bridge_validation"]["max_cs2_over_c2"],
       contract=rel(dest),contract_sha256=sha(dest),known_failed_native_knots_in_closure=gaps)
     record["admitted"].append(case);cases.append(case)
    except (ValueError,RuntimeError) as exc:
     trial.update(status="REJECTED",exception_type=type(exc).__name__,error=str(exc))
    record["lower_join_trials"].append(trial)
    if len(record["admitted"])==2:break
   records.append(record)
   print(json.dumps({"source":entry["id"],"admitted":[{"label":c["label"],"n":c["upper_start_density"],"kappa":c["kappa"]} for c in record["admitted"]]}),flush=True)
  after={rel(p):sha(p) for p in files};assert before==after
  report=dict(schema="ns_cowling.H_causal_closure_admission.v1",scope=plan["scope"],plan_sha256=sha(BASE/"ADMISSION_PLAN.json"),
    reference=rel(REFERENCE),reference_sha256=sha(REFERENCE),left_state=left,inputs=records,
    source_hashes_before=before,source_hashes_after=after,source_bytes_unchanged=True,
    admitted_contract_count=len(cases),stellar_observables_computed=False)
  save(BASE/"FULL_CONTRACT_ADMISSION.json",report)
  tov=dict(schema="ns_cowling.H_causal_background_plan.v1",scope="Finite 1.4 solar-mass TOV inventory of opt-in lower closures, with sharp source doublets retained",
    created_at_utc=stamp(),cases=cases,maximum_contracts=8,mass_msun=1.4,
    density_scan_fm3=[.08,.16,.32,.5,.8,1.0],selection="first sampled ascending target-mass bracket; no stability inference",
    settings=dict(r0_km=.001,rtol=1e-10,atol=1e-12,max_step_h=.002),
    tighter_settings=dict(r0_km=.0005,rtol=2e-11,atol=2e-13,max_step_h=.001),
    numerical_sensitivity="Primary contracts only; tighter solve at the same central density",
    mass_root_xtol=1e-10,mass_target_absolute_tolerance_msun=1e-7,
    admission=rel(BASE/"FULL_CONTRACT_ADMISSION.json"),admission_sha256=sha(BASE/"FULL_CONTRACT_ADMISSION.json"),
    source_and_code_sha256=after)
  save(BASE/"TOV_PLAN.json",tov)
  receipt.update(status="COMPLETED",actual_exit_code=0,finished_at_utc=stamp(),wall_seconds=time.monotonic()-start,
    admitted_contract_count=len(cases),source_bytes_unchanged=True)
  save(BASE/"ADMISSION_PROCESS.json",receipt)
 except BaseException as exc:
  receipt.update(status="FAILED",actual_exit_code=1,finished_at_utc=stamp(),wall_seconds=time.monotonic()-start,
    exception_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
  save(BASE/"ADMISSION_PROCESS.json",receipt);raise
if __name__=="__main__":
 a=argparse.ArgumentParser();a.add_argument("--expected-whole-sha",required=True);args=a.parse_args();main(args.expected_whole_sha)
