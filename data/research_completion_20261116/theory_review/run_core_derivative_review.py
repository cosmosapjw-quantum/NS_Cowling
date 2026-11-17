"""Independent off-equilibrium derivative checks, not target-table matching."""
import hashlib,json,platform
from pathlib import Path
import numpy as np
from compose_eos import gm1_uniform,bbb_uniform,apr_uniform,skyrme_uniform
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def canonical(d):
    if "epsilon_baryons" in d:return d["epsilon_baryons"],d["mu_n"],d["mu_p"],d["pressure_baryons"]
    if "epsilon" in d:return d["epsilon"],d["mu_n"],d["mu_p"],d["pressure"]
    return d["energy_density_mev_fm3"],d["mu_n_mev"],d["mu_p_mev"],d["pressure_mev_fm3"]
def generator(model):
    if model=="GM_GM1":return lambda nn,np_:canonical(gm1_uniform.baryons(nn+np_,np_/(nn+np_)))
    if model=="BBB_BHF2":return lambda nn,np_:canonical(bbb_uniform.baryons(nn+np_,np_/(nn+np_)))
    if model.startswith("APR_"):return lambda nn,np_:canonical(apr_uniform.nuclear_state(nn,np_,branch=model[4:]))
    return lambda nn,np_:canonical(skyrme_uniform.nuclear_state(model,nn,np_))
result={"claim":"OFF_EQUILIBRIUM_FUNCTIONAL_DERIVATIVE_CHECK; finite samples, no global proof",
    "model_checks":{},"software":{"python":platform.python_version(),"numpy":np.__version__}}
for model in ("GM_GM1","BBB_BHF2","APR_low","APR_high","RG_KDE0v1","RG_SLY4"):
    solve=generator(model);rows=[]
    for density in (.113,.237,.517,.913):
        for yp in (.03,.17,.41):
            x=np.array([density*(1-yp),density*yp])
            e,mn,mp,p=solve(*x);errors=[]
            for scale in (1e-4,5e-5):
                grads=[]
                for i in (0,1):
                    h=scale*x[i];q=np.zeros(2);q[i]=h
                    grad=(-solve(*(x+2*q))[0]+8*solve(*(x+q))[0]-8*solve(*(x-q))[0]+solve(*(x-2*q))[0])/(12*h)
                    grads.append(grad)
                errors.append((np.array(grads)-[mn,mp]).tolist())
            rows.append({"n_fm3":density,"yp":yp,"mu_derivative_errors_MeV":errors,
                         "euler_residual_MeV_fm3":e+p-x@[mn,mp]})
    maximum=max(abs(v) for r in rows for errors in r["mu_derivative_errors_MeV"] for v in errors)
    result["model_checks"][model]={"states":rows,"max_abs_mu_error_MeV":maximum,"pass_tolerance_MeV":2e-5,"status":"PASS_SCOPED" if maximum<2e-5 else "FAIL"}
    print(model,maximum,result["model_checks"][model]["status"],flush=True)
assert all(r["status"]=="PASS_SCOPED" for r in result["model_checks"].values())
result["input_sha256"]={name:hashlib.sha256((ROOT/"python_frontend/compose_eos"/name).read_bytes()).hexdigest() for name in ("gm1_uniform.py","bbb_uniform.py","apr_uniform.py","skyrme_uniform.py","cold_leptons.py")}
(HERE/"core_derivative_review.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
