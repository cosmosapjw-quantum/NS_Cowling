"""Reproduce why an interior-only beta bracket is insufficient at tiny density.

Run from repo root with PYTHONPATH=python_frontend. This compares a clearly
specified restricted bracket against the repaired KKT solution, without
replacing or modifying the imported physics module.
"""
import json
from scipy.optimize import brentq
from compose_eos.gppva_rmf import baryons, beta_state
from compose_eos.cold_leptons import ELECTRON_MASS, free_fermion

results=[]
for name,n in [('FSU2',1e-8),('DD2',1e-10)]:
    def residual(y):
        s=baryons(name,n,y)
        mu_e=free_fermion(n*y,ELECTRON_MASS)[2]
        return s['mu_n']-s['mu_p']-mu_e
    lo,hi=1e-12,.5
    row={'model':name,'density_fm3':n,'restricted_bracket':[lo,hi],
         'restricted_endpoint_residuals_MeV':[residual(lo),residual(hi)]}
    try:
        row['restricted_solution']=brentq(residual,lo,hi)
    except ValueError as exc:
        row['restricted_error']=str(exc)
    repaired=beta_state(name,n)
    row['repaired']={key:repaired[key] for key in
        ['Y_p','composition_branch','beta_residual_mev','charge_residual','total_euler_residual']}
    results.append(row)
print(json.dumps({'scope':'Homogeneous boundary counterexamples; no claim of physical outer-crust ground state',
                  'results':results},indent=2,allow_nan=False))
