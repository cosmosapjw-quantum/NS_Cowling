"""Independent H review: raw native/output byte and semantic checks; no production mutations."""
from pathlib import Path
import json,hashlib,subprocess,collections,time
import numpy as np
ROOT=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117")
DATA=ROOT/"data/low_density_crust_20261117";OUT=DATA/"review"
EXE_SHA="21f00ecd9cd492e37329e5690d342842c19025d4b649f97f99e84e592522c189"
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check(condition,message):
 if not condition:errors.append(message)
def table(path):
 path=Path(path);n=np.loadtxt(path/"eos.nb",skiprows=2)
 lines=(path/"eos.thermo").read_text().splitlines()
 mass=float(lines[0].split()[0]);q=np.array([[float(v) for v in l.split()] for l in lines[1:] if l.strip()])
 return dict(n=n,p=n*q[:,3],epsilon=n*mass*(q[:,9]+1),mu=mass*(q[:,5]+1))
def compare_table(tab,states,label):
 n=np.array([s["n_b"] for s in states]);check(len(n)==len(tab["n"]),label+" rows")
 if len(n)!=len(tab["n"]):return
 for key,field in [("n","n_b"),("p","pressure"),("epsilon","epsilon"),("mu","mu_B")]:
  a=np.array([s[field] for s in states]);err=float(np.max(abs(tab[key]-a)/np.maximum(1.,abs(a))))
  max_residuals[key]=max(max_residuals[key],err);check(err<5e-13,label+" "+key+" changed")
errors=[];max_residuals=dict(n=0.,p=0.,epsilon=0.,mu=0.)
native=[*sorted((DATA/"low_extension/states").glob("**/STATE.json")),*sorted((DATA/"low_extension/repaired").glob("*/states/**/STATE.json")),*sorted((DATA/"support_refinement/states").glob("**/STATE.json")),*sorted((DATA/"support_refinement/bounded_anchor_control").glob("**/STATE.json"))]
check(len(native)==len(set(native))==272,"native invocation count")
statuses=collections.Counter();new_accepted=0;native_bytes={}
for p in native:
 s=load(p);lines=(p.parent/"stdout.txt").read_text().splitlines();v=list(map(float,lines[2].split()))
 check(s["input_sha256"]==sha(p.parent/"input.txt"),str(p)+" input hash")
 check(s["points"]==128 and len(s["profile"])==128 and s["executable_sha256"]==EXE_SHA,str(p)+" mesh/executable")
 check(lines[0]=="ok" and s["status"]==lines[1].split()[0],str(p)+" native status")
 check(s["actual_exit_code"]==0,str(p)+" actual exit")
 for i,k in [(0,"epsilon"),(1,"pressure"),(2,"mu_B"),(14,"optimizer_mode"),(17,"optimizer_success")]:check(s[k]==v[i],str(p)+" stdout "+k)
 expected=(s["status"] in ("STATIONARY","UNIFORM_STATIONARY") and s["optimizer_success"]==1 and s["muon_threshold_crossed"]==0 and s["baryon_error"]<1e-10 and s["charge_error"]<1e-10)
 check(s["accepted_stationary"]==expected,str(p)+" admission")
 check(not s["evaluation_only"],str(p)+" evaluation-only counted")
 check(abs((s["epsilon"]+s["pressure"])/s["n_b"]-s["mu_B"])<1e-9,str(p)+" Euler identity")
 statuses[s["status"]]+=1;new_accepted+=bool(s["accepted_stationary"]);native_bytes[str(p.relative_to(ROOT))]=sha(p)
by_model={};matches=[];reused=0
for model in ["DDME2","DDMEX_ORIGINAL_TABLE_RHO0","NL3_PROVIDER_PRECISION"]:
 index=DATA/"low_extension/repaired"/model/"FINAL_N128_STATES.json";states=load(index)
 for s in states:
  check(load(s["state_path"])==s,model+" index differs from retained original")
  check(s["points"]==128,model+" mixed mesh")
  if "/crust_refinement_20261117/" in s["state_path"]:reused+=1
 requested=np.array([s["n_b"] for s in states]);check(np.all(np.diff(requested)>0),model+" ordering")
 check(len(states)==len(load(DATA/"low_extension/PLAN.json")["configurations"][model]["densities"]),model+" requested node dropped")
 summary=load(DATA/"low_extension/repaired/tables_v2"/model/"N128/SUMMARY.json")
 exported=0;block_ranges=[]
 for b in summary["blocks"]:
  gen=b["generation"];src=[load(z["path"]) for z in gen["source_states"]]
  inds=[int(np.argmin(abs(requested-s["n_b"]))) for s in src]
  check(inds==list(range(inds[0],inds[-1]+1)),model+" failed requested gap bridged")
  check(all(s["accepted_stationary"] for s in src),model+" rejected state exported")
  for z in gen["source_states"]:check(sha(z["path"])==z["sha256"],model+" source hash")
  for name,h in gen["table_sha256"].items():check(sha(Path(gen["compose_path"])/name)==h,model+" source table hash")
  compare_table(table(gen["compose_path"]),src,model+" source block");exported+=len(src)
  block_ranges.append([src[0]["n_b"],src[-1]["n_b"]])
  for m in b["matches"]:
   if m["status"]=="GENERATED_VERIFIED_FINITE_DOMAIN":matches.append((model,m,table(gen["compose_path"])))
 by_model[model]=dict(requested=len(states),accepted=sum(s["accepted_stationary"] for s in states),failed=[s["n_b"] for s in states if not s["accepted_stationary"]],exported=exported,block_ranges=block_ranges)
support=load(DATA/"support_refinement/REGENERATION.json")
prepared=load(DATA/"support_refinement/PREPARED_MANIFEST.json")["source_records"]
for run in support["runs"]:
 states=[load(p["state_path"]) if p["accepted"] else None for p in run["points_results"]]
 check(all(s is None or (s["points"]==128 and s["accepted_stationary"]) for s in states),run["model"]+" support admission")
 ranges=[];exported=0
 for rec in prepared:
  if rec["model"]!=run["model"]:continue
  a,b=rec["grid_indices"];block=states[a:b]
  check(all(s is not None for s in block),run["model"]+" gap bridged")
  for name,h in rec["source_table_sha256"].items():check(sha(ROOT/rec["source_path"]/name)==h,run["model"]+" support table hash")
  compare_table(table(ROOT/rec["source_path"]),block,run["model"]+" source")
  ranges.append([block[0]["n_b"],block[-1]["n_b"]]);exported+=len(block)
 by_model[run["model"]]=dict(requested=len(states),accepted=sum(s is not None for s in states),failed=[run["density_grid_fm3"][i] for i,s in enumerate(states) if s is None],exported=exported,block_ranges=ranges)
p=DATA/"support_refinement/matching/GPPVA_NL3wrL55/continued_G_lower_adjacent_branch/N128/run0_component_000__component_000/RESULT.json"
m=load(p);src=table(ROOT/next(r["source_path"] for r in prepared if r["model"]=="GPPVA_NL3wrL55"));matches.append(("GPPVA_NL3wrL55",m,src))
match_results=[]
for model,m,source in matches:
 native_rep=m.get("native_report") or m["v2_report"]["native_report"]
 path=Path(m.get("compose_path") or m["output_path"]);path=path if path.is_absolute() else ROOT/path
 tab=table(path);i=native_rep["low_phase_rows"]-1
 check(native_rep["physical_energy_shift_mev"]==0 and native_rep["is_exact_maxwell"],model+" matching semantics")
 check(np.isclose(tab["n"][i],native_rep["n_low"],rtol=1e-13) and np.isclose(tab["n"][i+1],native_rep["n_high"],rtol=1e-13),model+" interface location")
 for row in range(i):
  j=int(np.argmin(abs(source["n"]-tab["n"][row])))
  check(abs(source["n"][j]-tab["n"][row])<1e-15,model+" inserted retained source knot")
  for key in ["p","epsilon","mu"]:check(abs(source[key][j]-tab[key][row])<5e-13*max(1,abs(source[key][j])),model+" retained row "+key)
 delta_e=tab["epsilon"][i+1]-tab["epsilon"][i];delta_n=tab["n"][i+1]-tab["n"][i]
 maxwell=abs(delta_e-native_rep["mu_transition"]*delta_n)/max(abs(delta_e),1e-300)
 check(maxwell<1e-12,model+" raw Maxwell identity")
 hashes=m.get("table_sha256") or m["output_sha256"]
 for name,h in hashes.items():check(sha(path/name)==h,model+" output hash")
 check(m["admission"]["status"]=="ADMITTED",model+" serialized admission")
 match_results.append(dict(model=model,rows=len(tab["n"]),lower_support=float(tab["n"][0]),n_low=float(tab["n"][i]),n_high=float(tab["n"][i+1]),pressure=native_rep["pressure_transition"],maxwell_relative=maxwell,output_sha256=hashes))
check(reused==9,"G reused count")
protected=["fortran_eos/src/compose_eos_tf_models.f90","fortran_eos/src/compose_eos_spatial_tf.f90","fortran_eos/tools/tf_candidate_driver.f90","python_frontend/compose_eos/tf_candidate.py","data/new_model_maxwell_20261117/sources/PARAMETER_INPUTS.json"]
changed=subprocess.check_output(["git","diff","--name-only","905b845ad0f14d2f1e7f67f9eb802691144850cd","--",*protected],cwd=ROOT,text=True).splitlines();check(not changed,"protected physics identity changed")
out=dict(scope="Independent reread of raw native stdout/input/STATE, connected request index and raw CompOSE thermodynamics; no new native solve",status="PASS" if not errors else "FAIL",native_invocations=len(native),native_status_counts=dict(statuses),accepted_new_native_states=new_accepted,reused_exact_G_states=reused,models=by_model,source_semantic_max_relative_residuals=max_residuals,verified_matches=match_results,protected_changed=changed,protected_current_sha256={p:sha(ROOT/p) for p in protected},native_state_sha256=native_bytes,errors=errors)
(OUT/"SOURCE_REVIEW.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps({k:v for k,v in out.items() if k not in ["native_state_sha256","verified_matches","protected_current_sha256"]},indent=2))
raise SystemExit(bool(errors))
