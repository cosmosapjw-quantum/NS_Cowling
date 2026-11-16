from pathlib import Path
import json,numpy as np
p=Path(__file__).resolve().parent
a=json.loads((p/"surface_fit_variants.json").read_text())
b=json.loads((p/"surface_fit_variants_321.json").read_text())
rows=[]
for model,d in a["models"].items():
    for key,v in d["variants"].items():
        w=b["models"][model]["variants"][key]
        rows.append({"model":model,"variant":key,
            "sigma0_relative_delta":w["sigma0_MeV_fm2"]/v["sigma0_MeV_fm2"]-1,
            "b_absolute_delta":w["b"]-v["b"],
            "max_sigma_probe_relative_delta":float(np.max(abs(np.array(w["sigma_probe_MeV_fm2"])/v["sigma_probe_MeV_fm2"]-1))),
            "max_binding_prediction_delta_MeV":float(np.max(abs(np.array(w["prediction_binding_MeV"])-v["prediction_binding_MeV"]))),
            "all_success":v["optimizer_success"] and w["optimizer_success"]})
summary={"grid_comparison":[161,321],"cases":rows,
    "max_abs_sigma0_relative_delta":max(abs(x["sigma0_relative_delta"]) for x in rows),
    "max_abs_b_delta":max(abs(x["b_absolute_delta"]) for x in rows),
    "max_sigma_probe_relative_delta":max(x["max_sigma_probe_relative_delta"] for x in rows),
    "max_binding_prediction_delta_MeV":max(x["max_binding_prediction_delta_MeV"] for x in rows)}
(p/"surface_fit_grid_comparison.json").write_text(json.dumps(summary,indent=2)+"\n")
print({k:v for k,v in summary.items() if k!="cases"})
print("lowYp_variant_vs_old_max_rel",max(abs(v["sigma_probe_MeV_fm2"][0]/d["previous_2408_geometric"]["sigma_probe_MeV_fm2"][0]-1) for d in a["models"].values() for v in d["variants"].values()))
print("endpoint_ratios",[v["endpoint_halving_error_ratio"] for v in json.loads((p/"two_phase_numerical_review.json").read_text())["models"].values()])
