"""Cold GM1 uniform matter from Glendenning--Moszkowski (1991), Table II.

C_sigma,C_omega,C_rho are (g/m)^2 in fm^2. Scalar potential is
S^2/(2 C_sigma)+b*m*S^3/3+c*S^4/4 in fm^-4, with S=g_sigma*sigma.
The isovector vertex is g_rho/2. A common nucleon mass 939 MeV is explicit.
CompOSE GM1 has npe matter; npemu is an optional, separately labelled extension.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import brentq
from .cold_leptons import HBARC, free_fermion, leptons

MASS = 939.0
C_SIGMA, C_OMEGA, C_RHO = 11.79, 7.149, 4.411
B, C = 0.002947, -0.001070
SOURCE = "https://escholarship.org/uc/item/6sw424k1"


def scalar_density(density: float, mass_fm: float) -> float:
    if density == 0:
        return 0.0
    momentum = (3*math.pi**2*density)**(1/3)
    energy = math.hypot(momentum, mass_fm)
    return mass_fm/(2*math.pi**2)*(momentum*energy-mass_fm**2*math.asinh(momentum/mass_fm))


def baryons(density: float, proton_fraction: float) -> dict[str, float]:
    """Stationary scalar field and analytic nucleon chemical potentials."""
    if not math.isfinite(density) or density <= 0 or not math.isfinite(proton_fraction) or not 0 <= proton_fraction <= 1:
        raise ValueError("positive finite density and 0 <= proton_fraction <= 1 required")
    nn, np_ = density*(1-proton_fraction), density*proton_fraction
    mass = MASS/HBARC

    def residual(scalar):
        effective = mass-scalar
        return scalar/C_SIGMA+B*mass*scalar**2+C*scalar**3-scalar_density(nn,effective)-scalar_density(np_,effective)

    # Positive effective-mass stationary branch. Multiple physical roots are
    # ranked by energy, rather than choosing a root by iteration history.
    grid = np.linspace(0.0, mass*(1-1e-10), 65)
    roots = []
    left, fl = grid[0], residual(grid[0])
    for right in grid[1:]:
        fr = residual(right)
        if fl*fr <= 0:
            roots.append(brentq(residual,left,right,xtol=1e-13))
        left, fl = right, fr
    if not roots:
        raise RuntimeError("no positive-effective-mass GM1 stationary scalar solution")

    def energy(scalar):
        effective = (mass-scalar)*HBARC
        kinetic = free_fermion(nn,effective)[0]+free_fermion(np_,effective)[0]
        field = scalar**2/(2*C_SIGMA)+B*mass*scalar**3/3+C*scalar**4/4
        field += C_OMEGA*density**2/2+C_RHO*(np_-nn)**2/8
        return kinetic+HBARC*field

    scalar = min(roots,key=energy)
    effective = (mass-scalar)*HBARC
    mn = free_fermion(nn,effective)[2]+HBARC*(C_OMEGA*density-C_RHO*(np_-nn)/4)
    mp = free_fermion(np_,effective)[2]+HBARC*(C_OMEGA*density+C_RHO*(np_-nn)/4)
    eps = energy(scalar)
    return {"epsilon_baryons":eps,"pressure_baryons":nn*mn+np_*mp-eps,
            "mu_n":mn,"mu_p":mp,"effective_mass_mev":effective,
            "scalar_residual_fm3":residual(scalar)}


def beta_state(density: float, *, include_muons: bool = False) -> dict[str,float]:
    """Neutral cold beta equilibrium at the requested density."""
    def charge(yp):
        state = baryons(density,yp)
        lepton = leptons(state["mu_n"]-state["mu_p"],include_muons)
        return lepton["n_e"]+lepton["n_mu"]-density*yp
    yp = brentq(charge,1e-12,0.5,xtol=1e-13)
    state = baryons(density,yp)
    chemical = state["mu_n"]-state["mu_p"]
    lepton = leptons(chemical,include_muons)
    eps = state["epsilon_baryons"]+lepton["epsilon_leptons"]
    pressure = state["pressure_baryons"]+lepton["pressure_leptons"]
    return dict(state,**lepton,n_b=density,Y_n=1-yp,Y_p=yp,
                Y_e=lepton["n_e"]/density,Y_mu=lepton["n_mu"]/density,
                mu_e=chemical,epsilon=eps,pressure=pressure,
                charge_residual_fm3=density*yp-lepton["n_e"]-lepton["n_mu"],
                euler_residual=eps+pressure-density*state["mu_n"])
