from pathlib import Path
import sys,json,copy,hashlib
ROOT=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117")
G=Path("/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117")
sys.path.insert(0,str(G/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope,HermiteInterval
out=ROOT/"data/low_density_crust_20261117/stellar/ddmex_comparison";out.mkdir(parents=True,exist_ok=True)
H=json.loads((ROOT/"data/low_density_crust_20261117/stellar/h_final/FULL_CONTRACT_ADMISSION.json").read_text())
entry=next(x for x in H["inputs"] if x["crust"]=="DDMEX_ORIGINAL_TABLE_RHO0")
primary,nextjoin=entry["admitted"]
GB=json.loads((ROOT/"data/low_density_crust_20261117/stellar/g_baseline/FULL_CONTRACT_ADMISSION.json").read_text())
gb=next(x for x in GB["inputs"] if x["crust"]=="DDMEX_ORIGINAL_TABLE_RHO0")["admitted"][0]
c=json.loads((ROOT/primary["contract"]).read_text());bar=WholeCrustBarotrope.from_contract(c)
j=next(i for i,n in enumerate(c["upper"]["n"]) if abs(n-.007)<1e-14)
left=bar.bridge.contract()["endpoint_states"][0]
right={k:c["upper"][a][j] for k,a in [("n","n"),("p","p"),("epsilon","epsilon")]}
h=HermiteInterval(left["n"],left["p"],left["epsilon"],right["n"],right["p"],right["epsilon"])
assert h.validate()["status"]=="PASS"
control=copy.deepcopy(c)
control["lower_bridge"]=h.contract()
for key in ["n","p","epsilon"]:control["upper"][key]=c["upper"][key][j:]
control["upper"]["sharp_interfaces"]=[i-j for i in c["upper"]["sharp_interfaces"]]
control["provenance"].update(choice="fixed G join density 0.007 to isolate N64-to-N128 source change before observables",source_row_start=primary["source_row_start"]+j)
WholeCrustBarotrope.from_contract(control)
path=out/"contracts/H_DDMEX_N128_G_join/REPRESENTATION.json";path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(control,indent=2)+"\n")
cases=[
 dict(id="H_DDMEX_N128_primary",contract=primary["contract"],mesh=128,lower_join=.003,role="primary actual extended source"),
 dict(id="G_DDMEX_N64_primary",contract=gb["contract"],mesh=64,lower_join=.007,role="G fixed-outer baseline"),
 dict(id="H_DDMEX_N128_G_join",contract=str(path.relative_to(ROOT)),mesh=128,lower_join=.007,role="same join density, mesh/source refinement control"),
 dict(id="H_DDMEX_N128_next_join",contract=nextjoin["contract"],mesh=128,lower_join=.004,role="next admissible source knot selected before observables")]
for c in cases:c["contract_sha256"]=hashlib.sha256((ROOT/c["contract"]).read_bytes()).hexdigest()
plan=dict(phase="FROZEN_BEFORE_OBSERVABLES",physics_root=str(G),physics_commit="905b845ad0f14d2f1e7f67f9eb802691144850cd",mass_msun=1.4,
 scope="Conditional electron-only TF-source + fixed E outer/surface + named Hermite lower closure; equilibrium perfect-fluid Cowling; finite frequency scan",
 cases=cases,frequencies_hz=[50,80,120,180,250,350,500,700,1000,1250,1500,1750,2000,2250,2500,3000,3500],
 root_selection={"slow":"first two sign-change roots detected in fixed scan","rapid":"first one sign-change root detected in fixed scan"},
 max_roots={"slow":2,"rapid":1},continuation="seven points at 88%-112% of primary frequency; require exactly one bracket",
 acceptance={"determinant_abs":1e-6,"join_relative":1e-6,"surface_relative":1e-10,"traction_relative":1e-10,"kinematic_abs":1e-10},
 bounded_repair="one tighter solve per rejected root, original retained",
 numerical_control="primary at fixed central density, r0/2, tighter tolerances and surface_h*0.3",
 nonclaims=["complete spectrum","radial stability","global geometry minimum","crust microphysics in phenomenological bridge"])
(out/"PLAN.json").write_text(json.dumps(plan,indent=2)+"\n")
print(json.dumps(plan,indent=2))
