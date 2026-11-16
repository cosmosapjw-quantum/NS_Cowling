"""Cold spin-saturated uniform Skyrme matter with explicit n,p,e,mu composition.

Energies include nucleon rest masses. Densities are fm^-3, energies MeV and
energy densities/pressures MeV fm^-3. HBARC carries hbar*c in MeV fm.
Uniform states omit the gradient, spin-orbit, Coulomb and pairing terms.
The mathematical uniform branch is returned even beyond its causal domain.
"""
from dataclasses import asdict, dataclass
import math

from scipy.optimize import brentq
from .cold_leptons import lepton_state as _shared_lepton_state

HBARC = 197.3269804
MN = 939.5654
MP = 938.2720
ME = 0.51099895
MMU = 105.6583755
KINETIC_COEFFICIENT = 20.73553  # hbar^2/(2m), MeV fm^2, common fitted nucleon mass
FERMI_FACTOR = (3.0 * math.pi**2)**(2.0 / 3.0)


@dataclass(frozen=True)
class SkyrmeParameters:
    model_id: str
    t0: float
    t1: float
    t2: float
    t3: float
    x0: float
    x1: float
    x2: float
    x3: float
    alpha: float
    source: str
    parameter_precision: str


PARAMETERS = {
    "RG_KDE0v1": SkyrmeParameters(
        "RG_KDE0v1", -2553.08, 411.70, -419.87, 14603.6,
        0.6483, -0.3472, -0.9268, 0.9475, 0.1673,
        "https://arxiv.org/abs/nucl-th/0505071 Table 3",
        "Published Table 3 precision; alpha=0.1673"),
    "RG_SLY4": SkyrmeParameters(
        "RG_SLY4", -2488.91, 486.82, -546.39, 13777.0,
        0.834, -0.344, -1.0, 1.354, 1.0/6.0,
        "https://arxiv.org/abs/astro-ph/0006135 Table 1",
        "Published Douchin-Haensel Table 1 rounded coefficients"),
}


def parameter_provenance(model_id: str) -> dict:
    record = asdict(PARAMETERS[model_id])
    record.update(
        kinetic_coefficient_mev_fm2=KINETIC_COEFFICIENT,
        kinetic_convention_source="https://www.fuw.edu.pl/~dobaczew/hfodd/pb208-a.out.txt",
        neutron_rest_mass_mev=MN, proton_rest_mass_mev=MP,
        electron_mass_mev=ME, muon_mass_mev=MMU, hbar_c_mev_fm=HBARC)
    return record


def nuclear_state(model_id: str, neutron_density: float, proton_density: float) -> dict:
    """Evaluate the full asymmetric functional and its analytic density gradients."""
    p = PARAMETERS[model_id]
    nn, np = float(neutron_density), float(proton_density)
    if not all(math.isfinite(x) and x >= 0.0 for x in (nn, np)) or nn + np <= 0:
        raise ValueError("finite nonnegative species densities with positive total required")
    n = nn + np
    tau_n = 0.6 * FERMI_FACTOR * nn**(5.0/3.0)
    tau_p = 0.6 * FERMI_FACTOR * np**(5.0/3.0)
    tau = tau_n + tau_p
    grad_tau = [FERMI_FACTOR * nn**(2.0/3.0), FERMI_FACTOR * np**(2.0/3.0)]
    density_square = nn*nn + np*np
    a0 = (2.0+p.x0)*n*n - (2.0*p.x0+1.0)*density_square
    a3 = (2.0+p.x3)*n*n - (2.0*p.x3+1.0)*density_square
    momentum_isoscalar = (p.t1*(2.0+p.x1)+p.t2*(2.0+p.x2))/8.0
    momentum_isovector = -(p.t1*(2.0*p.x1+1.0)-p.t2*(2.0*p.x2+1.0))/8.0
    binding = (KINETIC_COEFFICIENT*tau + p.t0*a0/4.0 + p.t3*n**p.alpha*a3/24.0
               + momentum_isoscalar*n*tau + momentum_isovector*(nn*tau_n+np*tau_p))
    chemical = []
    for density, tau_q, derivative, mass in zip((nn,np),(tau_n,tau_p),grad_tau,(MN,MP)):
        chemical.append(
            mass + KINETIC_COEFFICIENT*derivative
            + p.t0*((2+p.x0)*n-(2*p.x0+1)*density)/2
            + p.t3*(p.alpha*n**(p.alpha-1)*a3
                    + 2*n**p.alpha*((2+p.x3)*n-(2*p.x3+1)*density))/24
            + momentum_isoscalar*(tau+n*derivative)
            + momentum_isovector*(tau_q+density*derivative))
    energy = MN*nn+MP*np+binding
    pressure = nn*chemical[0]+np*chemical[1]-energy
    return dict(n_b_fm3=n, binding_energy_density_mev_fm3=binding,
                energy_density_mev_fm3=energy, pressure_mev_fm3=pressure,
                mu_n_mev=chemical[0], mu_p_mev=chemical[1])


def lepton_state(chemical_potential_mev: float, mass_mev: float) -> dict:
    """Adapt the shared exact zero-temperature lepton EOS to named outputs."""
    mu, mass = float(chemical_potential_mev), float(mass_mev)
    if not math.isfinite(mu) or not math.isfinite(mass) or mass <= 0:
        raise ValueError("finite chemical potential and positive finite mass required")
    density, energy, pressure = _shared_lepton_state(mu, mass)
    return dict(density_fm3=density, energy_density_mev_fm3=energy,
                pressure_mev_fm3=pressure)


def beta_equilibrium(model_id: str, density_fm3: float, *, include_muons: bool = True) -> dict:
    """Minimize cold uniform matter over charge-neutral beta-equilibrium candidates.

    Root candidates are compared by total energy. This provides an explicit
    check of multiple composition roots on the stated 0<=Yp<=1/2 interval.
    """
    density = float(density_fm3)
    if not math.isfinite(density) or density <= 0:
        raise ValueError("positive finite baryon density required")

    def evaluate(proton_fraction):
        state = nuclear_state(model_id, density*(1-proton_fraction), density*proton_fraction)
        mu = state["mu_n_mev"]-state["mu_p_mev"]
        electron = lepton_state(mu, ME)
        muon = lepton_state(mu, MMU) if include_muons else lepton_state(0, MMU)
        residual = density*proton_fraction-electron["density_fm3"]-muon["density_fm3"]
        return state, mu, electron, muon, residual

    grid = [0.5*i/64 for i in range(65)]
    values = [evaluate(y)[-1] for y in grid]
    candidates = []
    for i in range(64):
        if values[i] == 0:
            candidates.append(grid[i])
        if values[i]*values[i+1] < 0:
            candidates.append(brentq(lambda y:evaluate(y)[-1], grid[i],grid[i+1],
                                     xtol=1e-15, rtol=1e-14))
    if values[-1] == 0:
        candidates.append(0.5)
    if not candidates:
        raise ValueError("no neutral beta-equilibrium candidate on 0<=Yp<=1/2")
    states = []
    for yp in candidates:
        state, mu, electron, muon, charge = evaluate(yp)
        energy = state["energy_density_mev_fm3"]+electron["energy_density_mev_fm3"]+muon["energy_density_mev_fm3"]
        pressure = state["pressure_mev_fm3"]+electron["pressure_mev_fm3"]+muon["pressure_mev_fm3"]
        states.append(dict(
            model_id=model_id, n_b_fm3=density, Y_n=1-yp, Y_p=yp,
            Y_e=electron["density_fm3"]/density, Y_mu=muon["density_fm3"]/density,
            mu_n_mev=state["mu_n_mev"], mu_p_mev=state["mu_p_mev"],
            mu_e_mev=mu, mu_mu_mev=mu if include_muons else None,
            energy_density_mev_fm3=energy, pressure_mev_fm3=pressure,
            charge_residual_fm3=charge,
            euler_residual_mev_fm3=energy+pressure-density*state["mu_n_mev"],
            candidate_count=len(candidates),
            composition_regime="npemu" if include_muons else "npe"))
    return min(states,key=lambda state:state["energy_density_mev_fm3"])

def compose_table(model_id: str, densities, *, include_muons: bool = True):
    """Construct a composition-complete CompOSE container from solved uniform states."""
    import numpy as np
    from .eos_data import EOSData

    density_array = np.asarray(densities, dtype=float)
    if (density_array.ndim != 1 or density_array.size < 2 or
            not np.all(np.isfinite(density_array)) or np.any(density_array <= 0) or
            np.any(np.diff(density_array) <= 0)):
        raise ValueError("at least two finite, positive, strictly increasing densities required")
    rows = [beta_equilibrium(model_id, n, include_muons=include_muons) for n in density_array]
    column = lambda name: np.array([row[name] for row in rows])
    zero = np.zeros(len(rows))
    table = EOSData(
        nb=density_array, mn=MN, mp=MP,
        Q1=column("pressure_mev_fm3")/density_array, Q2=zero.copy(),
        Q3=column("mu_n_mev")/MN-1, Q4=-column("mu_e_mev")/MN,
        Q5=zero.copy(), Q6=column("energy_density_mev_fm3")/(density_array*MN)-1,
        Q7=column("energy_density_mev_fm3")/(density_array*MN)-1,
        q_extra=zero.copy(), Iphase=np.full(len(rows),3,dtype=int),
        Yn=column("Y_n"), Yp=column("Y_p"), Ye=column("Y_e"), Ymu=column("Y_mu"),
        Aav=np.ones(len(rows)), Zav=zero.copy(), Yh=zero.copy(),
        composition_file_present=True, composition_available=True, beta_equilibrated=True)
    return table, rows


def sound_speed_ratio(model_id: str, density_fm3: float, *, relative_step: float = 1e-4) -> float:
    """Centered pressure derivative along beta equilibrium, divided by exact mu_n.

    The result is c_s^2/c^2. Step refinement is required for a claimed crossing.
    """
    n = float(density_fm3)
    if not math.isfinite(relative_step) or not 0 < relative_step < 0.1:
        raise ValueError("relative derivative step must be in (0,0.1)")
    h = n*relative_step
    center = beta_equilibrium(model_id, n)
    plus = beta_equilibrium(model_id, n+h)["pressure_mev_fm3"]
    minus = beta_equilibrium(model_id, n-h)["pressure_mev_fm3"]
    return (plus-minus)/(2*h*center["mu_n_mev"])
