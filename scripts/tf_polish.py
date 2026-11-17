"""Numerical restarts of the unchanged constrained TF functional."""
from __future__ import annotations
import numpy as np
from compose_eos.spatial_tf_equilibrium import solve_screened_cell


def profile_initial(old,points):
    p=old["profile"];x=(np.arange(points)+.5)/points
    values={k:np.interp(x,p["x"],p[k]) for k in ("nn","np","ne")}
    w=np.diff(np.linspace(0,1,points+1)**old["dimension"])
    factor=old["n_b"]/(w@(values["nn"]+values["np"]))
    values["nn"]*=factor;values["np"]*=factor
    values["ne"]*=float(w@values["np"])/max(float(w@values["ne"]),1e-30)
    return {**values,"radius":old["radius_fm"]}


def polish_profile(old,points=None):
    """Restart a candidate; every result must pass the original KKT criteria.

    The objective stopping tolerance is not a physical acceptance criterion.
    An energy-floor restart can recover SLSQP's success flag. A nearly absent
    species with a positive insertion cost may be initialized on its exact
    lower bound; all densities are then released in the next full minimization.
    """
    points=points or old["points"];attempts=[]
    def run(start,tol,label):
        result=solve_screened_cell(old["model_id"],old["n_b"],points=points,
            dimension=old["dimension"],seed=old["seed"],radius_bounds=(5.,180.),
            initial=start,ftol=tol)
        result["restart_kind"]=label;attempts.append(result);return result
    def valid(r):
        return r["numerical_status"] in ("STATIONARY","UNIFORM_STATIONARY")
    row=run(profile_initial(old,points),1e-13,"REFINED_PROFILE")
    if not valid(row):
        row=run(profile_initial(row,points),1e-14,"TIGHT_OBJECTIVE_STOPPING")
    if not valid(row) and row["max_chemical_KKT_residual_MeV"]<.02 and row["radius_stationarity_MeV"]<2e-4:
        row=run(profile_initial(row,points),2e-12,"ENERGY_PRECISION_FLOOR_RESTART")
    if not valid(row):
        initial=profile_initial(row,points);p=row["profile"];changed=False
        for density,chemical,offset in [("nn","mu_n",0.),("np","mu_p",row["mu_Q_lepton"])]:
            a=np.array(p[density]);gap=np.array(p[chemical])+offset-row["mu_B"]
            mask=(a<1e-7)&(gap>.005)
            changed=changed or bool(np.any(mask & (a>0)))
            initial[density][mask]=0.
        if changed:
            w=np.array(p["weights"])
            factor=row["n_b"]/(w@(initial["nn"]+initial["np"]))
            initial["nn"]*=factor;initial["np"]*=factor
            row=run(initial,1e-13,"INACTIVE_SPECIES_BOUND_SEED")
            if not valid(row) and row["max_chemical_KKT_residual_MeV"]<.02 and row["radius_stationarity_MeV"]<2e-4:
                row=run(profile_initial(row,points),2e-12,"BOUND_SEED_ENERGY_PRECISION_RESTART")
    accepted=[a for a in attempts if valid(a)]
    chosen=min(accepted or attempts,key=lambda r:r["epsilon"])
    result=dict(chosen)
    result["restart_attempts"]=attempts
    return result
