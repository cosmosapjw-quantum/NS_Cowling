"""Independent phase-aware physical diagnostics of the direct native state tables."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
sys.path.insert(0,str(ROOT/"python_frontend"))
from benchmark_native_dense_eos_20261118 import load_native,phase_runs
D=ROOT/"data/native_dense_interpolation_20261118/aggregated"
D.mkdir(parents=True,exist_ok=True)
MODELS={
    "rmf":["NL3_PROVIDER_PRECISION","PK1_CANONICAL","TM1_CANONICAL",
           "DDMEX_PROVIDER_ROUNDED_N0","DDMEX_ORIGINAL_TABLE_RHO0"],
    "core":["APR_APR","BBB_BHF2","GM_GM1","RG_KDE0v1","RG_SLY4"]}
result=dict(schema="ns_cowling.native_thermo_physics_audit.v1",source="direct new native CompOSE tables",
            status="RUNNING",models=[])
for family,names in MODELS.items():
    for model in names:
        s=load_native(family,model)
        groups=[]
        for a,b in phase_runs(s["phase"]):
            n,p,e,mu=[s[key][a:b] for key in ("n","p","epsilon","mu")]
            euler=e+p-n*mu
            normalized=np.abs(euler)/np.maximum(1,np.abs(e+p))
            with np.errstate(divide="ignore",invalid="ignore"):
                de=np.gradient(e,n,edge_order=2)
                dp=np.gradient(p,n,edge_order=2)
                dmu=np.gradient(mu,n,edge_order=2)
                firstlaw=(de-mu)/np.maximum(1,np.abs(mu))
                gibbs=(dp-n*dmu)/np.maximum(1,np.abs(dp))
                speed=dp/de
            inside=slice(2,len(n)-2) if len(n)>6 else slice(None)
            record=dict(phase=int(s["phase"][a]),rows=len(n),n_fm3=[float(n[0]),float(n[-1])],
              max_scaled_Euler=float(np.max(normalized)),
              median_scaled_firstlaw_fd=float(np.median(abs(firstlaw[inside]))),
              p95_scaled_firstlaw_fd=float(np.quantile(abs(firstlaw[inside]),.95)),
              median_scaled_gibbs_duhem_fd=float(np.median(abs(gibbs[inside]))),
              p95_scaled_gibbs_duhem_fd=float(np.quantile(abs(gibbs[inside]),.95)),
              min_cs2_from_native_fd=float(np.nanmin(speed)),
              max_cs2_from_native_fd=float(np.nanmax(speed)),
              acausal_native_row_count=int(np.sum(speed>1)),
              source_nonmonotone_n=int(np.sum(np.diff(n)<=0)),
              source_nonmonotone_mu=int(np.sum(np.diff(mu)<=0)),
              source_nonmonotone_p=int(np.sum(np.diff(p)<=0)))
            groups.append(record)
        result["models"].append(dict(family=family,model=model,source_sha256=s["sha256"],phases=groups))
result["status"]="DERIVED_FROM_DIRECT_NATIVE_TABLES"
result["models_count"]=len(result["models"])
result["total_source_rows"]=sum(z["rows"] for item in result["models"] for z in item["phases"])
result["max_relative_euler"]=max(z["max_scaled_Euler"] for item in result["models"] for z in item["phases"])
(D/"NATIVE_THERMO_AUDIT.json").write_text(json.dumps(result,indent=2)+"\n")
print("MODEL_COUNT",result["models_count"],"ROW_COUNT",result["total_source_rows"],
      "MAX_EULER",result["max_relative_euler"])
for item in result["models"]:
    print(" ",item["family"],item["model"],"cs2_max",
          round(max(z["max_cs2_from_native_fd"] for z in item["phases"]),6),
          "p95_firstlaw",max(z["p95_scaled_firstlaw_fd"] for z in item["phases"]),
          "acausal_rows",sum(z["acausal_native_row_count"] for z in item["phases"]))
assert result["total_source_rows"]==25612
assert result["models_count"]==10
assert result["max_relative_euler"]<5e-10
