"""Isolated species-rest-energy-restored SPG / BSk22 seam experiment.

No production module changes. OriginalcommonM kinetic and coupling parameters
and candidate surface are retained. Only explicit physical rest-energy linear
terms are added BEFORE minimization. This is a distinct candidate model.
"""
from pathlib import Path
from functools import lru_cache
import argparse,hashlib,json,math,platform
import numpy as np
from scipy.optimize import minimize,least_squares,brentq
from scipy.interpolate import PchipInterpolator
from compose_eos.candidate_gdfm_cld import reduced_cell,MASS
from compose_eos import bsk22_analytic as bsk
from compose_eos.bsk22_analytic import NEUTRON_DRIP_DENSITY
def bsk_state(n):
    if n<=.002:return bsk.state(n)
    # Explicit isolated use of the published full-crust energy fit beyond
    # the production helper's intentionally narrow lower-segment API.
    if not .002<n<=.004:raise ValueError("isolated extension restricted to.004")
    e,de,*_=bsk._energy(n);eps=n*(bsk.NEUTRON_REST_MEV+e)
    p,chi=bsk._pressure(eps*bsk.MEV_FM3_TO_G_CM3)
    return {"n_b":n,"epsilon":eps,"pressure":p,"pressure_from_energy_derivative":n*de,
        "mu_b_euler":(eps+p)/n,"mu_b_energy_derivative":bsk.NEUTRON_REST_MEV+e+de,
        "representation":"isolated explicit published full-crust fit extension above source2e-3seam",
        "pressure_fit_thermodynamic_relative_residual":(p-n*de)/p}
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
MN=939.5654;MP=938.2720;DN=MN-MASS;DP=MP-MASS;DIFF=DP-DN
def corrected_cell(model,n,z):
    energy,gradient,detail=reduced_cell(model,n,z)
    t,ni,y=z
    return energy+DN+DIFF*t*y,gradient+np.array([DIFF*y,0.,DIFF*t]),detail

@lru_cache(maxsize=20000)
def restored_state(model,n):
    low=np.array([1e-7,max(.035,n*(1+1e-6)),.001]);high=np.array([1.,.25,.499])
    def fun(z):return corrected_cell(model,n,z)[:2]
    records=[]
    for initial in ([.9,.16,.3],[.5,.15,.2],[.2,.14,.12]):
        opt=minimize(fun,initial,jac=True,method="L-BFGS-B",bounds=list(zip(low,high)),
                     options={"ftol":1e-14,"gtol":1e-10,"maxiter":500,"maxls":50})
        z=opt.x.copy();boundary=z[0]>1-1e-6 and fun(z)[1][0]<0
        if boundary:
            fit=least_squares(lambda q:fun([1.,*q])[1][1:],z[1:],bounds=(low[1:]+1e-12,high[1:]-1e-12),
                              jac="3-point",xtol=1e-12,ftol=1e-12,gtol=1e-12)
            z=np.r_[1.,fit.x]
        else:
            fit=least_squares(lambda q:fun(q)[1],np.clip(z,low+1e-10,high-1e-10),
                              bounds=(low+1e-12,high-1e-12),jac="3-point",xtol=1e-12,ftol=1e-12,gtol=1e-12)
            z=fit.x
        e,g,detail=corrected_cell(model,n,z)
        projected=g.copy()
        if z[0]>1-1e-8:projected[0]=max(g[0],0.)
        mu=detail["state"]["mu_n_gas"]+DN+z[0]*g[0]
        eps=n*(MASS+e);p=n*mu-eps
        margin=detail["state"]["gas_proton_insertion_margin_MeV"]+DIFF
        records.append({"n_B":n,"coordinates":z.tolist(),"energy_per_baryon":MASS+e,"epsilon":eps,
            "pressure":p,"mu_B":mu,"Y_e":float(z[0]*z[2]),"Y_n_free":float(1-z[0]),
            "gradient":g.tolist(),"projected_gradient_max":float(max(abs(projected))),
            "gas_proton_insertion_margin_MeV":margin,"muon_margin_MeV":detail["state"]["muon_insertion_margin_MeV"],
            "initial_minimize_success":bool(opt.success),"polish_success":bool(fit.success)})
    winner=min(records,key=lambda r:r["energy_per_baryon"])
    if winner["projected_gradient_max"]>1e-6 or winner["gas_proton_insertion_margin_MeV"]<0 or not winner["polish_success"]:
        raise RuntimeError("unresolved restored candidate "+str(winner))
    winner=dict(winner);winner["all_three_attempts"]=records
    return winner

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--models",nargs="+",default=[f"SPG_M{i}" for i in range(1,6)])
    parser.add_argument("--points",type=int,default=65);parser.add_argument("--bsk-max",type=float,default=.002);args=parser.parse_args()
    out={"claim":"SPECIES_REST_ENERGY_RESTORED_CANDIDATE_NOT_ORIGINAL_SP G_IMPLEMENTATION".replace("SP G","SPG"),
        "rest_energies_MeV":{"M0":MASS,"Mn":MN,"Mp":MP,"delta_n":DN,"delta_p":DP,"beta_difference_shift":DN-DP},
        "definition":"epsilon'=epsilon+delta_n*nn_total+delta_p*np_total before beta/composition minimization; kinetic masses/couplings/surface unchanged",
        "bsk_representation":"published analytic ETFSI fits; no independent microscopic ETFSI solve",
        "gibbs_primary_choice":"energy derivative P=n^2*d(E/A)/dn, mu=d(epsilon)/dn from BSk energyfit; sourceindependent pressurefit kept as sensitivity",
        "models":{}}
    out["bsk_density_max_fm3"]=args.bsk_max
    bsk_densities=np.unique(np.r_[np.geomspace(NEUTRON_DRIP_DENSITY,args.bsk_max,args.points),.002])
    bsk=[bsk_state(float(n)) for n in bsk_densities]
    out["bsk22_states"]=bsk
    for model in args.models:
        print("START",model,flush=True)
        ngrid=np.geomspace(.00012,max(.004,args.bsk_max*2),args.points)
        states=[restored_state(model,float(n)) for n in ngrid]
        pressures=np.array([r["pressure"] for r in states])
        assert np.all(np.diff(pressures)>0)
        inverse=PchipInterpolator(pressures,ngrid)
        def at_pressure(p):
            j=int(np.searchsorted(pressures,p));j=max(1,min(j,len(ngrid)-1))
            n=brentq(lambda n:restored_state(model,float(n))["pressure"]-p,ngrid[j-1],ngrid[j],xtol=2e-14)
            return restored_state(model,float(n))
        record={"states":states,"seam_at_n_0.002":restored_state(model,.002),"matches":{}}
        for label,pkey,mukey in (("energy_consistent","pressure_from_energy_derivative","mu_b_energy_derivative"),("source_pressure_fit","pressure","mu_b_euler")):
            samples=[]
            for b in bsk:
                pp=b[pkey]
                if pressures[0]<=pp<=pressures[-1]:
                    s=at_pressure(pp)
                    samples.append({"n_bsk":b["n_b"],"n_spg":s["n_B"],"pressure":pp,"delta_mu_spg_minus_bsk_MeV":s["mu_B"]-b[mukey]})
            roots=[]
            def gap(n):
                b=bsk_state(float(n));return at_pressure(b[pkey])["mu_B"]-b[mukey]
            for left,right in zip(samples[:-1],samples[1:]):
                if left["delta_mu_spg_minus_bsk_MeV"]*right["delta_mu_spg_minus_bsk_MeV"]<0:
                    n=brentq(gap,left["n_bsk"],right["n_bsk"],xtol=2e-14)
                    b=bsk_state(float(n));s=at_pressure(b[pkey])
                    roots.append({"n_bsk":n,"n_spg":s["n_B"],"pressure":b[pkey],"delta_mu_MeV":s["mu_B"]-b[mukey],
                        "spg_state":s,"bsk_state":b})
            record["matches"][label]={"samples":samples,"roots":roots,
                "min_delta_mu_MeV":min(r["delta_mu_spg_minus_bsk_MeV"] for r in samples),
                "max_delta_mu_MeV":max(r["delta_mu_spg_minus_bsk_MeV"] for r in samples)}
            print(model,label,"range",record["matches"][label]["min_delta_mu_MeV"],record["matches"][label]["max_delta_mu_MeV"],"roots",len(roots),flush=True)
        out["models"][model]=record
        out["code_sha256"]={name:hashlib.sha256((ROOT/"python_frontend/compose_eos"/name).read_bytes()).hexdigest() for name in ("candidate_gdfm_cld.py","gdfm_uniform.py","bsk22_analytic.py")}
        (HERE/"rest_mass_seam_experiment.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print("COMPLETE",len(out["models"]),flush=True)
if __name__=="__main__":main()
