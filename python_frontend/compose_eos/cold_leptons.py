"""Zero-temperature free leptons; energies in MeV and densities in fm^-3."""
from __future__ import annotations
import math
from scipy.special import hyp2f1

HBARC = 197.3269804
ELECTRON_MASS = 0.51099895
MUON_MASS = 105.6583755


def free_fermion(density: float, mass: float) -> tuple[float, float, float]:
    """Return total energy density, pressure and chemical potential (spin 1/2)."""
    if not math.isfinite(density) or density < 0 or not math.isfinite(mass) or mass <= 0:
        raise ValueError("finite nonnegative density and positive mass required")
    if density == 0:
        return 0.0, 0.0, mass
    momentum = HBARC * (3 * math.pi**2 * density)**(1/3)
    chemical = math.hypot(momentum, mass)
    # Integral form avoids cancellation near a species threshold.
    energy = mass * density * float(hyp2f1(-0.5, 1.5, 2.5, -(momentum/mass)**2))
    pressure = density * chemical - energy
    return energy, pressure, chemical


def lepton_state(chemical: float, mass: float) -> tuple[float, float, float]:
    """Return density, energy and pressure; an unpopulated species has zero values."""
    if not math.isfinite(chemical) or not math.isfinite(mass) or mass <= 0:
        raise ValueError("finite chemical potential and positive finite mass required")
    if chemical <= mass:
        return 0.0, 0.0, 0.0
    momentum = math.sqrt((chemical-mass)*(chemical+mass))
    density = (momentum/HBARC)**3/(3*math.pi**2)
    energy, pressure, _ = free_fermion(density, mass)
    return density, energy, pressure


def leptons(chemical: float, include_muons: bool) -> dict[str, float]:
    ne, ee, pe = lepton_state(chemical, ELECTRON_MASS)
    nm, em, pm = lepton_state(chemical, MUON_MASS) if include_muons else (0., 0., 0.)
    return {"n_e": ne, "n_mu": nm, "epsilon_leptons": ee+em, "pressure_leptons": pe+pm}
