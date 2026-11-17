"""Cold uniform APR low/high-density branches with explicit composition.

Hamiltonian and p1..p21: Akmal et al. PRC 58, 1804 (1998),
Schneider et al. PRC 100, 025803 (2019), arXiv:1901.09652.
Branch evaluation uses the published functional. The optional coexistence
calculation explicitly assumes local neutrality in each cold beta phase.
Densities fm^-3, energies/masses MeV, energy density and pressure MeV/fm^3.
"""
from __future__ import annotations

import math
from scipy.optimize import brentq
from .cold_leptons import HBARC, leptons

PARAMETERS = (337.2, -382.0, 89.8, .457, -59.0, -19.1, 214.6,
              -384.0, 6.4, 69.0, -33.0, .35, 0.0, 0.0, 287.0,
              -1.54, 175.0, -1.45, .32, .195, 0.0)
M_N = 939.5654133
M_P = 938.2720813


def nuclear_state(nn: float, np_: float, branch: str = "low") -> dict:
    """Evaluate one homogeneous phase, including nucleon rest energies."""
    if not all(math.isfinite(v) and v >= 0 for v in (nn, np_)) or nn + np_ <= 0:
        raise ValueError("require finite nonnegative species densities and n > 0")
    if branch not in ("low", "high"):
        raise ValueError("branch must be low or high")
    n = nn + np_
    if not math.isfinite(n) or not .01 <= n <= 2.0:
        raise ValueError("supported computational density range is .01 <= n <= 2 fm^-3")
    p = (0.0,) + PARAMETERS
    exp9 = math.exp(-(p[9]*n)**2)
    dexp9 = -2*p[9]**2*n*exp9
    g1 = -(p[1]+p[2]*n+p[6]*n*n+(p[10]+p[11]*n)*exp9)
    g2 = -(p[12]/n+p[7]+p[8]*n+p[13]*exp9)
    dg1 = -(p[2]+2*p[6]*n+p[11]*exp9+(p[10]+p[11]*n)*dexp9)
    dg2 = p[12]/n**2-p[8]-p[13]*dexp9
    if branch == "high":
        d1, d2 = n-p[19], n-p[20]
        e1, e2 = math.exp(p[18]*d1), math.exp(p[16]*d2)
        g1 -= (p[17]*d1+p[21]*d1*d1)*e1
        g2 -= (p[15]*d2+p[14]*d2*d2)*e2
        dg1 -= (p[17]+2*p[21]*d1+p[18]*(p[17]*d1+p[21]*d1*d1))*e1
        dg2 -= (p[15]+2*p[14]*d2+p[16]*(p[15]*d2+p[14]*d2*d2))*e2
    z = math.exp(-p[4]*n)
    kn = HBARC**2/(2*M_N)+(p[3]*n+p[5]*nn)*z
    kp = HBARC**2/(2*M_P)+(p[3]*n+p[5]*np_)*z
    if kn <= 0 or kp <= 0:
        raise ValueError("nonpositive inverse effective mass")
    taun = .6*(3*math.pi**2)**(2/3)*nn**(5/3)
    taup = .6*(3*math.pi**2)**(2/3)*np_**(5/3)
    common_n = z*(p[3]-p[4]*(p[3]*n+p[5]*nn))
    common_p = z*(p[3]-p[4]*(p[3]*n+p[5]*np_))
    delta = nn-np_
    common_u = 4*dg1*nn*np_+dg2*delta**2
    mun = (M_N+kn*(3*math.pi**2*nn)**(2/3)
           +(common_n+p[5]*z)*taun+common_p*taup
           +common_u+4*g1*np_+2*g2*delta)
    mup = (M_P+kp*(3*math.pi**2*np_)**(2/3)
           +common_n*taun+(common_p+p[5]*z)*taup
           +common_u+4*g1*nn-2*g2*delta)
    eps = M_N*nn+M_P*np_+kn*taun+kp*taup+4*g1*nn*np_+g2*delta**2
    return {"n_b":n,"epsilon":eps,"mu_n":mun,"mu_p":mup,
            "pressure":nn*mun+np_*mup-eps,"branch":branch,
            "mstar_n":HBARC**2/(2*kn),"mstar_p":HBARC**2/(2*kp)}


def beta_state(n: float, branch: str = "low", include_muons: bool = True) -> dict:
    """Charge-neutral beta equilibrium on the explicitly requested branch."""
    if not math.isfinite(n) or not .01 <= n <= 2.0:
        raise ValueError("supported computational density range is .01 <= n <= 2 fm^-3")
    def residual(yp):
        matter = nuclear_state(n*(1-yp),n*yp,branch)
        lep = leptons(matter["mu_n"]-matter["mu_p"],include_muons)
        return n*yp-lep["n_e"]-lep["n_mu"]
    yp = brentq(residual,1e-12,.5,xtol=1e-14)
    state = nuclear_state(n*(1-yp),n*yp,branch)
    lep = leptons(state["mu_n"]-state["mu_p"],include_muons)
    eps = state["epsilon"]+lep["epsilon_leptons"]
    pressure = state["pressure"]+lep["pressure_leptons"]
    return {**state,"epsilon":eps,"pressure":pressure,"Y_n":1-yp,"Y_p":yp,
            "Y_e":lep["n_e"]/n,"Y_mu":lep["n_mu"]/n,
            "mu_e":state["mu_n"]-state["mu_p"],
            "charge_residual_fm3":n*yp-lep["n_e"]-lep["n_mu"],
            "euler_residual_mev_fm3":eps+pressure-n*state["mu_n"],
            "species":"npemu" if include_muons else "npe"}

SOURCE = {
    "hamiltonian": "https://arxiv.org/html/1901.09652v2",
    "original": "https://doi.org/10.1103/PhysRevC.58.1804",
    "author_code": "https://bitbucket.org/andschn/apreos",
}


def maxwell_coexistence(include_muons: bool = True) -> dict:
    """Common tangent between locally neutral, cold beta APR branches.

    This is the one-conserved-charge Maxwell construction: each bulk phase
    is electrically neutral. It is not a globally neutral Gibbs construction
    with a shared electron chemical potential. Coulomb/interface energies
    are absent. The individual endpoint electron chemical potentials differ.
    """
    import numpy as np
    from scipy.optimize import least_squares

    def residual(densities):
        left = beta_state(float(densities[0]), "low", include_muons)
        right = beta_state(float(densities[1]), "high", include_muons)
        return [left["mu_n"]-right["mu_n"],
                (left["pressure"]-right["pressure"])/.2]
    solutions = []
    for seed in ((.19, .25), (.205, .235), (.22, .28)):
        fit = least_squares(residual, seed, bounds=([.12, .21], [.24, .36]),
                            xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=300)
        if not fit.success or max(abs(x) for x in residual(fit.x)) > 1e-7:
            raise RuntimeError("APR local-neutral Maxwell root did not converge")
        solutions.append(fit.x)
    densities = solutions[0]
    if max(float(np.max(abs(x-densities))) for x in solutions) > 1e-8:
        raise RuntimeError("APR Maxwell multistart roots disagree")
    left, right = [beta_state(float(n), branch, include_muons)
                   for n, branch in zip(densities, ("low", "high"))]
    if right["n_b"] <= left["n_b"]:
        raise RuntimeError("APR coexistence density jump is not positive")
    pressure = .5*(left["pressure"]+right["pressure"])
    chemical = .5*(left["mu_n"]+right["mu_n"])
    # A sampled minorant check supports this tangent; it is not a proof of a
    # global minimum over all density/composition phases.
    excess = []
    for density in np.linspace(.08, 1.34, 401):
        for branch in ("low", "high"):
            state = beta_state(float(density), branch, include_muons)
            excess.append(state["epsilon"]-(chemical*density-pressure))
    if min(excess) < -1e-7:
        raise RuntimeError("APR Maxwell tangent lies above a sampled pure branch")
    return {"construction": "local_neutrality_Maxwell_npemu" if include_muons
                else "local_neutrality_Maxwell_npe",
            "left": left, "right": right, "pressure": pressure, "mu_b": chemical,
            "pressure_residual_mev_fm3": left["pressure"]-right["pressure"],
            "chemical_residual_mev": left["mu_n"]-right["mu_n"],
            "multistart_max_density_difference_fm3":
                max(float(np.max(abs(x-densities))) for x in solutions),
            "tangent_check": {"n_min_fm3":.08, "n_max_fm3":1.34,
                "points_per_branch":401, "minimum_energy_excess_mev_fm3":min(excess),
                "status":"SAMPLED_MINORANT"},
            "global_phase_minimum_status": "NOT_PROVEN",
            "source_provider_phase_prescription": "GIBBS in official 2025 provider PDF; Maxwell is an explicit alternative"}


def equilibrium_state(n: float, coexistence: dict | None = None,
                      include_muons: bool = True) -> dict:
    """Pure or volume-averaged state on the explicitly local-neutral Maxwell EOS.

    Mixed states have no single local electron/proton chemical potential or
    effective mass. Those values are None; consult the endpoint phases instead.
    """
    if not math.isfinite(n) or not .01 <= n <= 2.0:
        raise ValueError("supported computational density range is .01 <= n <= 2 fm^-3")
    coexistence = maxwell_coexistence(include_muons) if coexistence is None else coexistence
    left, right = coexistence["left"], coexistence["right"]
    expected = "npemu" if include_muons else "npe"
    if left["species"] != expected or right["species"] != expected:
        raise ValueError("coexistence species differ from requested lepton sector")
    if n <= left["n_b"]:
        return {**beta_state(n, "low", include_muons), "high_phase_volume_fraction":0.0}
    if n >= right["n_b"]:
        return {**beta_state(n, "high", include_muons), "high_phase_volume_fraction":1.0}
    fraction = (n-left["n_b"])/(right["n_b"]-left["n_b"])
    eps = (1-fraction)*left["epsilon"]+fraction*right["epsilon"]
    composition = {key:((1-fraction)*left["n_b"]*left[key]
                       +fraction*right["n_b"]*right[key])/n
                   for key in ("Y_n","Y_p","Y_e","Y_mu")}
    return {"n_b":n, "epsilon":eps, "pressure":coexistence["pressure"],
            "mu_n":coexistence["mu_b"], "mu_p":None, "mu_e":None,
            "mstar_n":None, "mstar_p":None, **composition, "branch":"mixed",
            "high_phase_volume_fraction":fraction, "species":expected,
            "charge_residual_fm3":n*(composition["Y_p"]-composition["Y_e"]-composition["Y_mu"]),
            "euler_residual_mev_fm3":eps+coexistence["pressure"]-n*coexistence["mu_b"]}
