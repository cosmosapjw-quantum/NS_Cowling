"""Author-parameterized three-meson GDFM homogeneous matter.

Scurto et al. 2024, arXiv:2402.15548, Eqs. 1-20. Author Zenodo
11084025 parameter candidates are recorded separately from the historical
four-meson SPG proxy. The archived rho coefficients define Gamma_rho,
with g_rho=2*Gamma_rho in the paper's explicit g_rho/2 vertex convention.
This convention is tested against the archived nuclear matter parameters.
"""
from __future__ import annotations

from functools import lru_cache
import json
import math
from pathlib import Path

from scipy.optimize import brentq

from .cold_leptons import HBARC, free_fermion, leptons
from .gm1_uniform import scalar_density
from .native import project_root

MASS = 938.9
MODEL_IDS = tuple(f"SPG_M{i}" for i in range(1, 6))


@lru_cache(maxsize=1)
def parameters() -> dict:
    path = project_root() / "data/research_generator_20261116/spg_author_candidates.json"
    return json.loads(path.read_text())["models"]


def coupling(parameter: dict, density: float, saturation: float) -> tuple[float, float]:
    """Return Gamma and dGamma/dn in fm^3 using the model's own rho_sat."""
    x = density / saturation
    a, b, c, d = (parameter[name] for name in ("a", "b", "c", "d"))
    exponential = math.exp(-c*x)
    value = a + (b+d*x**3)*exponential
    derivative = (3*d*x*x-c*(b+d*x**3))*exponential/saturation
    return value, derivative


def baryons(model_id: str, density: float, proton_fraction: float,
            *, rho_vertex_scale: float = 1.0) -> dict[str, float]:
    """Return stationary homogeneous nuclear matter including rearrangement.

    rho_vertex_scale=1 uses the author archive Gamma_rho convention.
    Scale 0.5 is exposed for the regression comparison to the paper-literal
    interpretation of the raw archive coefficients.
    """
    if model_id not in MODEL_IDS:
        raise ValueError(f"supported source GDFM models are {MODEL_IDS}")
    if (not math.isfinite(density) or density <= 0 or
            not math.isfinite(proton_fraction) or not 0 <= proton_fraction <= 1 or
            not math.isfinite(rho_vertex_scale) or rho_vertex_scale <= 0):
        raise ValueError("finite positive density and vertex scale, 0<=Yp<=1 required")
    model = parameters()[model_id]
    coefficients, derivatives = {}, {}
    for meson in ("sigma", "omega", "rho"):
        g, dg = coupling(model["couplings"][meson], density, model["rho0_fm3"])
        if meson == "rho":
            g, dg = rho_vertex_scale*g, rho_vertex_scale*dg
        factor = (HBARC/model["meson_masses_MeV"][meson])**2
        coefficients[meson] = factor*g*g
        derivatives[meson] = 2*factor*g*dg
    cs, cw, cr = (coefficients[key] for key in ("sigma", "omega", "rho"))
    nn, np_ = density*(1-proton_fraction), density*proton_fraction
    asymmetry = np_-nn
    mass_fm = MASS/HBARC

    def residual(scalar):
        effective = mass_fm-scalar
        return scalar/cs-scalar_density(nn,effective)-scalar_density(np_,effective)

    scalar = brentq(residual, 0.0, mass_fm*(1-1e-12), xtol=1e-13)
    effective = (mass_fm-scalar)*HBARC
    en, _, mun_kin = free_fermion(nn,effective)
    ep, _, mup_kin = free_fermion(np_,effective)
    eps = en+ep+HBARC*(scalar*scalar/cs+cw*density*density+cr*asymmetry*asymmetry)/2
    scalar_density_total = scalar/cs
    rearrangement = HBARC*(derivatives["omega"]*density*density
                          +derivatives["rho"]*asymmetry*asymmetry
                          -derivatives["sigma"]*scalar_density_total**2)/2
    mun = mun_kin+HBARC*(cw*density-cr*asymmetry)+rearrangement
    mup = mup_kin+HBARC*(cw*density+cr*asymmetry)+rearrangement
    return {"epsilon_baryons":eps, "pressure_baryons":nn*mun+np_*mup-eps,
            "mu_n":mun,"mu_p":mup,"effective_mass_mev":effective,
            "rearrangement_mev":rearrangement,"scalar_residual_fm3":residual(scalar)}


def symmetry_energy(model_id: str, density: float, *, rho_vertex_scale: float = 1.0) -> float:
    """Exact quadratic symmetry energy of the three-meson functional."""
    state = baryons(model_id,density,0.5,rho_vertex_scale=rho_vertex_scale)
    model = parameters()[model_id]
    gamma, _ = coupling(model["couplings"]["rho"],density,model["rho0_fm3"])
    momentum = HBARC*(3*math.pi**2*density/2)**(1/3)
    kinetic = momentum*momentum/(6*math.hypot(momentum,state["effective_mass_mev"]))
    potential = (rho_vertex_scale*gamma)**2*HBARC**3*density/(2*model["meson_masses_MeV"]["rho"]**2)
    return kinetic+potential


def beta_state(model_id: str, density: float, *, include_muons: bool = True) -> dict[str, float]:
    """Neutral cold npe or npe-muon equilibrium of the author GDFM candidate."""
    def charge(yp):
        state = baryons(model_id,density,yp)
        lepton = leptons(state["mu_n"]-state["mu_p"],include_muons)
        return lepton["n_e"]+lepton["n_mu"]-density*yp

    yp = brentq(charge,1e-12,0.5,xtol=1e-13)
    state = baryons(model_id,density,yp)
    chemical = state["mu_n"]-state["mu_p"]
    lepton = leptons(chemical,include_muons)
    eps = state["epsilon_baryons"]+lepton["epsilon_leptons"]
    pressure = state["pressure_baryons"]+lepton["pressure_leptons"]
    return {"n_B":density,"Yp":yp,"Yn":1-yp,"Ye":lepton["n_e"]/density,
            "Ymu":lepton["n_mu"]/density,"epsilon":eps,"pressure":pressure,
            "mu_n":state["mu_n"],"mu_p":state["mu_p"],"mu_e":chemical,
            "charge_residual":density*yp-lepton["n_e"]-lepton["n_mu"],
            "euler_residual":eps+pressure-density*state["mu_n"],
            "scalar_residual_fm3":state["scalar_residual_fm3"]}
