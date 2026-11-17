"""Reproduce numerical failures and their final-equation-preserving repairs."""
from pathlib import Path
import importlib.util,json
import numpy as np
from compose_eos.spatial_tf import RadialGrid,ELECTRON_MASS
from compose_eos.spatial_tf_equilibrium import screened_electrons,solve_screened_cell

ROOT=Path(__file__).resolve().parent

def historical(version):
    name="compose_eos.spatial_tf_equilibrium_"+version
    spec=importlib.util.spec_from_file_location(name,ROOT/"source_versions"/("spatial_tf_equilibrium_"+version+".py"))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

results=[]
for version,operation in [
    ("v2",lambda m:m.solve_screened_cell("GPPVA_DD2",.03,points=32)),
    ("v3",lambda m:m.screened_electrons(RadialGrid(32,1),
            .003*(1+.8*np.cos(np.pi*RadialGrid(32,1).x)),35.,
            initial=(np.zeros(32),ELECTRON_MASS)))]:
    try:operation(historical(version));red={"expected_failure_reproduced":False}
    except Exception as exc:red={"expected_failure_reproduced":True,"exception":type(exc).__name__,"message":str(exc)}
    results.append({"version":version,"red":red})

cell=solve_screened_cell("GPPVA_DD2",.03,points=32)
g=RadialGrid(32,1);ne,v,mu,res=screened_electrons(g,.003*(1+.8*np.cos(np.pi*g.x)),35.,
    initial=(np.zeros(32),ELECTRON_MASS))
results[0]["green"]={"numerical_status":cell["numerical_status"],"KKT_MeV":cell["max_chemical_KKT_residual_MeV"]}
results[1]["green"]={"PDE_residual":res,"charge_residual":g.mean(.003*(1+.8*np.cos(np.pi*g.x))-ne),"mu_e":mu}
print(json.dumps(results,indent=2,allow_nan=False))
assert all(r["red"]["expected_failure_reproduced"] for r in results)
assert cell["numerical_status"]=="STATIONARY"
assert res<2e-10
