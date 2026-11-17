import json
from pathlib import Path
r=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117")
for p in ["data/low_density_crust_20261117/low_extension/MANIFEST.json","data/low_density_crust_20261117/support_refinement/PREPARED_MANIFEST.json","data/low_density_crust_20261117/support_refinement/ACCEPTED_SOURCE_STATES.json","data/low_density_crust_20261117/low_extension/repaired/tables_v2/DDME2/N128/SUMMARY.json","data/low_density_crust_20261117/low_extension/repaired/tables_v2/DDME2/N128/block_01/match_00_00/RESULT.json"]:
 x=json.loads((r/p).read_text())
 print(p)
 if isinstance(x,dict):
  print({k: type(v).__name__+(" "+str(len(v)) if isinstance(v,(list,dict)) else " "+str(v)[:90]) for k,v in x.items()})
  if "blocks" in x: print("blocks",[{k:v for k,v in b.items() if k not in ["matches","source_states","thermo_v2"]} for b in x["blocks"]][:3])
 else:print(type(x).__name__,len(x),list(x[0]) if x and isinstance(x[0],dict) else "")
p=r/"data/low_density_crust_20261117/low_extension/repaired/DDME2/FINAL_N128_STATES.json";s=json.loads(p.read_text())[0]
print("state keys",list(s))
print("state small",{k:v for k,v in s.items() if k not in ["profile","parent","source_rows"] and not isinstance(v,(dict,list))})
