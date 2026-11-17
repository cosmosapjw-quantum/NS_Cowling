"""Two-charge Gibbs coexistence for cold uniform APR.

Each nuclear phase may carry charge. A common ideal electron/muon background
is included once. Surface and Coulomb finite-size terms are absent. This is
the bulk Gibbs limit, not a prediction of a particular pasta geometry.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import brentq,least_squares
from .apr_uniform import beta_state,nuclear_state
from .cold_leptons import leptons

PROVIDER_PHASE_SOURCE = {
    "url":"https://compose.obspm.fr/download/1D/NS/Classical/APR/eos.pdf",
    "date":"2025-03-14",
    "pdf_sha256":"f4d5f27e88d9f5752a95c39062f487cbe762e2e28aa7ee7f889196c26dd612af",
    "phase_description":"inhomogeneous mixed phase with a Gibbs construction",
    "provider_thermo_sha256":"45011d7e94607825b6f00f9e70ba48470aa21652a138aafaaa6b31bc0f2f291a",
    "provider_thermo_identity":"downloaded current file byte-identical to repository APR_APR/eos.thermo",
}


def _solve_fraction(chi: float, include_muons: bool,
                    seed=(.205,.0735,.237,.0588,138.0)) -> dict:
    if not math.isfinite(chi) or not 0 <= chi <= 1:
        raise ValueError("high phase volume fraction must be finite and in [0,1]")
    def components(z):
        nl,yl,nh,yh,mu_e=(float(v) for v in z)
        low=nuclear_state(nl*(1-yl),nl*yl,"low")
        high=nuclear_state(nh*(1-yh),nh*yh,"high")
        lepton=leptons(mu_e,include_muons)
        return low,high,lepton
    def residual(z):
        low,high,lepton=components(z)
        nl,yl,nh,yh,mu_e=z
        return [low["mu_n"]-high["mu_n"],
                low["mu_p"]-high["mu_p"],
                (low["pressure"]-high["pressure"])/.2,
                low["mu_n"]-low["mu_p"]-mu_e,
                ((1-chi)*nl*yl+chi*nh*yh-lepton["n_e"]-lepton["n_mu"])*5000]
    fit=least_squares(residual,seed,bounds=([.10,.001,.20,.001,100.],
                      [.25,.30,.36,.30,180.]),xtol=1e-13,ftol=1e-13,
                      gtol=1e-13,max_nfev=300)
    if not fit.success or max(abs(v) for v in residual(fit.x))>1e-7:
        raise RuntimeError("APR bulk Gibbs coexistence did not converge")
    nl,yl,nh,yh,mu_e=(float(v) for v in fit.x)
    if nh<=nl:
        raise RuntimeError("APR Gibbs phase densities are not ordered")
    low,high,lepton=components(fit.x)
    n=(1-chi)*nl+chi*nh
    np_total=(1-chi)*nl*yl+chi*nh*yh
    epsilon=(1-chi)*low["epsilon"]+chi*high["epsilon"]+lepton["epsilon_leptons"]
    pressure=.5*(low["pressure"]+high["pressure"])+lepton["pressure_leptons"]
    mu_b=.5*(low["mu_n"]+high["mu_n"])
    mu_p=mu_b-mu_e
    charge=np_total-lepton["n_e"]-lepton["n_mu"]
    return {"n_b":n,"epsilon":epsilon,"pressure":pressure,"mu_n":mu_b,"mu_p":mu_p,
        "mu_e":mu_e,"Y_n":1-np_total/n,"Y_p":np_total/n,"Y_e":lepton["n_e"]/n,
        "Y_mu":lepton["n_mu"]/n,"mstar_n":None,"mstar_p":None,
        "branch":"gibbs_mixed","species":"npemu" if include_muons else "npe",
        "high_phase_volume_fraction":chi,"n_low":nl,"n_high":nh,
        "Yp_low":yl,"Yp_high":yh,
        "charge_residual_fm3":charge,
        "euler_residual_mev_fm3":epsilon+pressure-n*mu_b,
        "gibbs_residuals":{"mu_n_difference_mev":low["mu_n"]-high["mu_n"],
            "mu_p_difference_mev":low["mu_p"]-high["mu_p"],
            "nuclear_pressure_difference_mev_fm3":low["pressure"]-high["pressure"],
            "beta_residual_mev":low["mu_n"]-low["mu_p"]-mu_e,
            "charge_residual_fm3":charge},
        "nuclear_pressure":.5*(low["pressure"]+high["pressure"]),
        "nuclear_phase_states":{"low":low,"high":high}}


def gibbs_state_at_volume_fraction(chi: float,include_muons: bool=True) -> dict:
    """Equilibrate both charged phases and the common lepton background."""
    return _solve_fraction(chi,include_muons)


def gibbs_coexistence(include_muons: bool=True) -> dict:
    """Find both globally neutral coexistence endpoints with multistart checks."""
    endpoints=[]
    disagreement=[]
    for chi in (0.0,1.0):
        samples=[_solve_fraction(chi,include_muons,seed) for seed in
                 ((.205,.0735,.237,.0588,138.),(.19,.06,.26,.07,140.),
                  (.22,.10,.28,.09,145.))]
        endpoints.append(samples[0])
        disagreement.extend(abs(s[k]-samples[0][k]) for s in samples[1:]
            for k in ("n_b","n_low","n_high","Yp_low","Yp_high"))
    if max(disagreement)>1e-8:
        raise RuntimeError("APR Gibbs endpoint multistart roots disagree")
    path=[gibbs_state_at_volume_fraction(float(chi),include_muons)
          for chi in np.linspace(0,1,21)]
    densities=np.array([s["n_b"] for s in path])
    if np.any(np.diff(densities)<=0):
        raise RuntimeError("APR Gibbs density is not monotone in volume fraction")
    # The grand-potential plane is checked against both nuclear branches over
    # a declared sampled domain. No all-phase/global-minimum proof is claimed.
    minimum_excess=math.inf
    for endpoint in endpoints:
        for n in np.linspace(.08,1.34,101):
            for yp in np.linspace(0,.5,51):
                for branch in ("low","high"):
                    state=nuclear_state(float(n*(1-yp)),float(n*yp),branch)
                    excess=(state["epsilon"]-endpoint["mu_n"]*n*(1-yp)
                            -endpoint["mu_p"]*n*yp+endpoint["nuclear_pressure"])
                    minimum_excess=min(minimum_excess,excess)
    if minimum_excess < -1e-7:
        raise RuntimeError("APR Gibbs tangent plane is above a sampled nuclear state")
    return {"construction":"bulk_two_charge_Gibbs","species":"npemu" if include_muons else "npe",
        "left":endpoints[0],"right":endpoints[1],
        "multistart_max_density_or_fraction_difference":max(disagreement),
        "volume_fraction_path":path,
        "sampled_tangent_plane":{"n_interval_fm3":[.08,1.34],"n_points":101,
            "Yp_interval":[0.,.5],"Yp_points":51,"branches":["low","high"],
            "endpoint_count":2,"minimum_excess_mev_fm3":minimum_excess,
            "status":"SAMPLED_MINORANT"},
        "global_phase_minimum_status":"NOT_PROVEN",
        "source_phase_prescription":PROVIDER_PHASE_SOURCE,
        "finite_size_terms":"surface and Coulomb energies omitted; bulk Gibbs limit"}


def gibbs_equilibrium_state(n: float,coexistence: dict|None=None,
                            include_muons: bool=True) -> dict:
    """Return a pure or Gibbs-mixed state at the requested mean density."""
    if not math.isfinite(n) or not .01<=n<=2.:
        raise ValueError("supported computational density range is .01 <= n <= 2 fm^-3")
    phase=gibbs_coexistence(include_muons) if coexistence is None else coexistence
    if phase["species"]!=("npemu" if include_muons else "npe"):
        raise ValueError("coexistence species differ from requested lepton sector")
    if n<=phase["left"]["n_b"]:
        return {**beta_state(n,"low",include_muons),"high_phase_volume_fraction":0.}
    if n>=phase["right"]["n_b"]:
        return {**beta_state(n,"high",include_muons),"high_phase_volume_fraction":1.}
    chi=brentq(lambda x:gibbs_state_at_volume_fraction(x,include_muons)["n_b"]-n,
               0.,1.,xtol=2e-13)
    state=gibbs_state_at_volume_fraction(chi,include_muons)
    if abs(state["n_b"]-n)>1e-11:
        raise RuntimeError("APR Gibbs mean-density inversion did not converge")
    return state
