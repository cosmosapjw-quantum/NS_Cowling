"""Supplemental sequence-endpoint numerical convergence; no EOS mutation."""
import csv,json
from pathlib import Path
from compose_eos.stellar_background import ContractEnthalpyEOS,solve_background
root=Path(__file__).resolve().parents[3]
out=Path(__file__).resolve().parent
eos=ContractEnthalpyEOS.from_path(root/"data/core_support_20261117/primary_v1/REPRESENTATION.json")
rows=[]
for nc in (.1,1.):
    base=solve_background(eos,nc,rtol=2e-10,atol=1e-12,max_step_h=.002)
    for label,cfg in (
        ("tight",dict(rtol=3e-11,atol=2e-13,max_step_h=.001)),
        ("center_10m",dict(rtol=2e-10,atol=1e-12,max_step_h=.002,r0_km=.01)),
        ("surface_h1e-8",dict(rtol=2e-10,atol=1e-12,max_step_h=.002,surface_h=1e-8))):
        b=solve_background(eos,nc,**cfg)
        rows.append(dict(nc=nc,label=label,mass_msun=b.mass_msun,radius_km=b.radius_km,
                         delta_mass_msun=b.mass_msun-base.mass_msun,
                         delta_radius_m=1000*(b.radius_km-base.radius_km)))
with (out/"endpoint_convergence.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(rows,indent=2))
