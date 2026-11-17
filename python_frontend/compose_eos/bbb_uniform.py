"""BBB2 (Paris+Urbana) published-table composition reconstruction.

Source Baldo, Bombaci & Burgio 1997, Table 1 and Eqs.16--22.
SNM/PNM energies are interpolated with PCHIP and parabolic asymmetry.
The source microscopic calculations extend to n=0.8 fm^-3; tabulated points
above this include the author's extrapolation. This implementation interpolates
rounded published values and records that approximation explicitly.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
from functools import lru_cache
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from .cold_leptons import leptons

MASS = 939.0
SOURCE = "https://arxiv.org/pdf/astro-ph/9707277"
TABLE_PATH = Path(__file__).resolve().parents[2]/"data/research_generator_20261116/bbb2_published_tables.json"


@lru_cache(maxsize=1)
def source_tables():
    table = json.loads(TABLE_PATH.resolve().read_text())
    symmetric = PchipInterpolator(table["n"],table["symmetric_EA"],extrapolate=False)
    neutron = PchipInterpolator(table["n"],table["neutron_EA"],extrapolate=False)
    return table,symmetric,neutron


def baryons(density: float, proton_fraction: float) -> dict[str,float]:
    table,symmetric,neutron = source_tables()
    if not math.isfinite(density) or not table["n"][0] <= density <= table["n"][-1]:
        raise ValueError("BBB2 published-table support is 0.08 <= n <= 1.5 fm^-3")
    if not math.isfinite(proton_fraction) or not 0 <= proton_fraction <= 1:
        raise ValueError("proton_fraction must be in [0,1]")
    e0 = float(symmetric(density))
    symmetry = float(neutron(density)-symmetric(density))
    derivative0 = float(symmetric.derivative()(density))
    derivative_symmetry = float(neutron.derivative()(density))-derivative0
    delta = 1-2*proton_fraction
    energy = e0+symmetry*delta**2
    derivative = derivative0+derivative_symmetry*delta**2
    common = MASS+energy+density*derivative
    difference = 4*symmetry*delta
    mn = common+proton_fraction*difference
    mp = common-(1-proton_fraction)*difference
    return {"epsilon_baryons":density*(MASS+energy),"pressure_baryons":density**2*derivative,
            "mu_n":mn,"mu_p":mp,"symmetry_energy_mev":symmetry}


def beta_state(density: float, *, include_muons: bool = True) -> dict[str,float]:
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
    return dict(state,**lepton,n_b=density,Y_n=1-yp,Y_p=yp,Y_e=lepton["n_e"]/density,
                Y_mu=lepton["n_mu"]/density,mu_e=chemical,epsilon=eps,pressure=pressure,
                charge_residual_fm3=density*yp-lepton["n_e"]-lepton["n_mu"],
                euler_residual=eps+pressure-density*state["mu_n"])
