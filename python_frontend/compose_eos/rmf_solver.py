"""RMF homogeneous-matter solver for the EOS390 CLD reconstruction.

Implements the density-dependent relativistic mean-field (DD-RMF) energy
density for homogeneous asymmetric nuclear matter, following:
  - Scurto, Pais, Gulminelli PRD 109 103015 (2024), eq. 9-13
  - Gögelein, van Dalen, Fuchs, Muther PRC 77 025802 (2008), eq. 40 + Table I

Two coupling functional forms are supported:
  - DD2_DDME2_TYYPEL: g(rho) = a + (b + d*x)*exp(-c*x),  x = rho/rho_sat
  - SPG_GOGELEIN: g(rho_B) = a + (b + d*x^3)*exp(-c*x),  x = rho_B

UNIT CONVENTION (critical): all internal calculations are in NATURAL UNITS
with hbar = c = 1, so:
  - masses, momenta, energies, fields in MeV
  - densities in MeV^3 (converted from fm^-3 input by * HBARC^3)
  - the kinetic/scalar-density integrals are evaluated with kF in MeV and
    M* in MeV, yielding densities in MeV^3 and energy densities in MeV^4
  - the FINAL energy density and pressure are converted to MeV/fm^3 by
    dividing by HBARC^3, so they are in the standard EOS units (MeV/fm^3).

Public input/output:
  - rho_B, rho_p, rho_sat: fm^-3
  - meson masses: MeV
  - energy density, pressure: MeV/fm^3
  - chemical potentials, fields, effective masses: MeV
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

HBARC = 197.3269804  # MeV*fm;  rho[MeV^3] = rho[fm^-3] * HBARC^3
M_N = 939.565
M_P = 938.272
M_E = 0.511
M_MU = 105.658
PI = math.pi


def _kF_MeV(rho_MeV3: float) -> float:
    """Fermi momentum (MeV) from density in MeV^3, spin-1/2: kF^3 = 3 pi^2 rho."""
    if rho_MeV3 <= 0:
        return 0.0
    return (3.0 * PI * PI * rho_MeV3) ** (1.0 / 3.0)


def _scalar_density_MeV3(kF: float, Mstar: float) -> float:
    """rho_s = (M*/(2 pi^2)) [ kF E - M*^2 ln((kF+E)/M*) ], in MeV^3.
    E = sqrt(kF^2 + M*^2)."""
    if kF <= 0:
        return 0.0
    E = math.sqrt(kF * kF + Mstar * Mstar)
    return (Mstar / (2.0 * PI * PI)) * (kF * E - Mstar * Mstar * math.log((kF + E) / Mstar))


def _kinetic_MeV4(kF: float, Mstar: float, g: int) -> float:
    """Kinetic energy density (MeV^4): g/(2 pi^2) ∫_0^kF k^2 sqrt(k^2+M*^2) dk.
    Analytic: (g/(16 pi^2)) [ kF(2kF^2+M*^2)E - M*^4 ln((kF+E)/M*) ]."""
    if kF <= 0:
        return 0.0
    E = math.sqrt(kF * kF + Mstar * Mstar)
    I = 0.125 * (kF * (2 * kF * kF + Mstar * Mstar) * E - Mstar ** 4 * math.log((kF + E) / Mstar))
    return (g / (2.0 * PI * PI)) * I


@dataclass
class RMFCouplings:
    a: float
    b: float
    c: float
    d: float = 0.0

    def value(self, x: float, form: str) -> float:
        if form == "DD2_DDME2_TYYPEL":
            return self.a + (self.b + self.d * x) * math.exp(-self.c * x)
        if form == "SPG_GOGELEIN":
            return self.a + (self.b + self.d * x ** 3) * math.exp(-self.c * x)
        raise ValueError(f"unknown functional form {form}")

    def dvalue_dx(self, x: float, form: str) -> float:
        if form == "DD2_DDME2_TYYPEL":
            return (self.d - self.c * (self.b + self.d * x)) * math.exp(-self.c * x)
        if form == "SPG_GOGELEIN":
            return (3.0 * self.d * x * x - self.c * (self.b + self.d * x ** 3)) * math.exp(-self.c * x)
        raise ValueError(f"unknown functional form {form}")


@dataclass
class RMFModel:
    key: str
    family: str
    functional_form: str
    meson_masses: dict
    rho_sat: float            # fm^-3
    couplings: dict
    saturation_targets: dict
    reference: str = ""

    @property
    def rho_sat_MeV3(self) -> float:
        return self.rho_sat * HBARC ** 3

    def x_of(self, rho_B_MeV3: float) -> float:
        if self.functional_form == "SPG_GOGELEIN":
            # Gögelein uses x = rho_B; rho_B in fm^-3, so x = rho_B_MeV3 / HBARC^3
            return rho_B_MeV3 / HBARC ** 3
        return rho_B_MeV3 / self.rho_sat_MeV3

    def g(self, meson: str, rho_B_MeV3: float) -> float:
        return self.couplings[meson].value(self.x_of(rho_B_MeV3), self.functional_form)

    def dg_drho_fm(self, meson: str, rho_B_MeV3: float) -> float:
        """dg/d(rho_B[fm^-3]). For Typel x=rho/rho_sat => dg/drho_fm = (1/rho_sat) dg/dx.
        For Gögelein x=rho_B[fm^-3] => dg/drho_fm = dg/dx."""
        x = self.x_of(rho_B_MeV3)
        dgdx = self.couplings[meson].dvalue_dx(x, self.functional_form)
        if self.functional_form == "SPG_GOGELEIN":
            return dgdx  # x = rho_B directly in fm^-3
        return dgdx / self.rho_sat  # x = rho/rho_sat

    def includes_delta(self) -> bool:
        return "delta" in self.couplings


@dataclass
class HomogeneousResult:
    rho_B: float; rho_3: float; rho_p: float; rho_n: float
    sigma0: float; omega0: float; rho03: float; delta0: float
    Mstar_n: float; Mstar_p: float
    energy_density: float; pressure: float
    mu_n: float; mu_p: float; rearrangement: float
    converged: bool; n_iter: int


def solve_homogeneous(model: RMFModel, rho_B: float, rho_p: float,
                      tol: float = 1e-8, max_iter: int = 1000) -> HomogeneousResult:
    """Self-consistent DD-RMF for homogeneous matter at (rho_B[fm^-3], rho_p[fm^-3])."""
    if rho_B <= 0:
        return HomogeneousResult(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, True, 0)
    rho_B_M = rho_B * HBARC ** 3
    rho_p_M = rho_p * HBARC ** 3
    rho_n_M = rho_B_M - rho_p_M
    has_delta = model.includes_delta()
    m_s = model.meson_masses["sigma"]
    m_w = model.meson_masses["omega"]
    m_r = model.meson_masses["rho"]
    m_d = model.meson_masses.get("delta", 0.0)
    g_s = model.g("sigma", rho_B_M)
    g_w = model.g("omega", rho_B_M)
    g_r = model.g("rho", rho_B_M)
    g_d = model.g("delta", rho_B_M) if has_delta else 0.0

    Mstar_n = M_N
    Mstar_p = M_P
    converged = False
    it = 0
    for it in range(max_iter):
        kF_n = _kF_MeV(rho_n_M)
        kF_p = _kF_MeV(rho_p_M)
        rs_n = _scalar_density_MeV3(kF_n, Mstar_n)
        rs_p = _scalar_density_MeV3(kF_p, Mstar_p)
        rs_total = rs_n + rs_p
        rs_3 = rs_p - rs_n
        sigma0 = g_s * rs_total / (m_s ** 2)
        omega0 = g_w * rho_B_M / (m_w ** 2)
        rho03 = g_r * (rho_p_M - rho_n_M) / (m_r ** 2)
        delta0 = g_d * rs_3 / (m_d ** 2) if has_delta else 0.0
        Mstar_n_new = M_N - g_s * sigma0 - (g_d * delta0 if has_delta else 0.0)
        Mstar_p_new = M_P - g_s * sigma0 + (g_d * delta0 if has_delta else 0.0)
        if Mstar_n_new <= 10 or Mstar_p_new <= 10:
            break
        dn = abs(Mstar_n_new - Mstar_n)
        dp = abs(Mstar_p_new - Mstar_p)
        alpha = 0.4
        Mstar_n = alpha * Mstar_n_new + (1 - alpha) * Mstar_n
        Mstar_p = alpha * Mstar_p_new + (1 - alpha) * Mstar_p
        if max(dn, dp) / M_N < tol:
            Mstar_n, Mstar_p = Mstar_n_new, Mstar_p_new
            converged = True
            break

    # final fields with converged masses
    kF_n = _kF_MeV(rho_n_M)
    kF_p = _kF_MeV(rho_p_M)
    rs_n = _scalar_density_MeV3(kF_n, Mstar_n)
    rs_p = _scalar_density_MeV3(kF_p, Mstar_p)
    rs_total = rs_n + rs_p
    rs_3 = rs_p - rs_n
    sigma0 = g_s * rs_total / (m_s ** 2)
    omega0 = g_w * rho_B_M / (m_w ** 2)
    rho03 = g_r * (rho_p_M - rho_n_M) / (m_r ** 2)
    delta0 = g_d * rs_3 / (m_d ** 2) if has_delta else 0.0

    Ekin_n = _kinetic_MeV4(kF_n, Mstar_n, 2)
    Ekin_p = _kinetic_MeV4(kF_p, Mstar_p, 2)
    Ekin = Ekin_n + Ekin_p

    # field energy (MeV^4): +1/2 m_w^2 w0^2 +1/2 m_r^2 rho03^2 -1/2 m_s^2 sigma0^2 -1/2 m_d^2 delta0^2
    Efield = (0.5 * m_w ** 2 * omega0 ** 2 + 0.5 * m_r ** 2 * rho03 ** 2
              - 0.5 * m_s ** 2 * sigma0 ** 2
              - (0.5 * m_d ** 2 * delta0 ** 2 if has_delta else 0.0))
    # interaction energy (MeV^4): scalar -g_s sigma0 rs_total (-g_d delta0 rs_3 if delta)
    # plus vector -g_w w0 rho_B - g_r rho03 rho_3
    Eint = (-(g_s * sigma0 * rs_total)
            - (g_d * delta0 * rs_3 if has_delta else 0.0)
            - g_w * omega0 * rho_B_M
            - g_r * rho03 * (rho_p_M - rho_n_M))

    E_MeV4 = Ekin + Efield + Eint
    E_MeV_fm3 = E_MeV4 / HBARC ** 3  # convert to MeV/fm^3

    # rearrangement Sigma^R (MeV): sum_kappa (dg_kappa/drho_B[fm]) * field * density
    rho_n = rho_B - rho_p
    dgs = model.dg_drho_fm("sigma", rho_B_M)
    dgw = model.dg_drho_fm("omega", rho_B_M)
    dgr = model.dg_drho_fm("rho", rho_B_M)
    sigma_R = (dgw * omega0 * rho_B + dgs * sigma0 * (rs_n + rs_p) / HBARC ** 3
               + dgr * rho03 * (rho_p - rho_n))
    if has_delta:
        dgd = model.dg_drho_fm("delta", rho_B_M)
        sigma_R += dgd * delta0 * (rs_p - rs_n) / HBARC ** 3

    En_n = math.sqrt(kF_n ** 2 + Mstar_n ** 2)
    En_p = math.sqrt(kF_p ** 2 + Mstar_p ** 2)
    # chemical potentials: mu_i = E_F,i + g_w omega0 +/- (g_r rho03)/2 - Sigma^R
    # (the vector omega enters POSITIVELY — it is repulsive; the rho isovector
    # is + for protons, - for neutrons; the rearrangement Sigma^R subtracts).
    # Sign convention follows Scurto 2024 eq. 13 with E_F = sqrt(kF^2+M*^2).
    mu_n = En_n + g_w * omega0 - 0.5 * g_r * rho03 - sigma_R
    mu_p = En_p + g_w * omega0 + 0.5 * g_r * rho03 - sigma_R

    # pressure (MeV^4 -> MeV/fm^3): P = sum mu_i rho_i - E
    P_MeV4 = (mu_n * rho_n_M + mu_p * rho_p_M) - E_MeV4
    P_MeV_fm3 = P_MeV4 / HBARC ** 3

    return HomogeneousResult(
        rho_B=rho_B, rho_3=rho_p - rho_n, rho_p=rho_p, rho_n=rho_B - rho_p,
        sigma0=sigma0, omega0=omega0, rho03=rho03, delta0=delta0,
        Mstar_n=Mstar_n, Mstar_p=Mstar_p,
        energy_density=E_MeV_fm3, pressure=P_MeV_fm3,
        mu_n=mu_n, mu_p=mu_p, rearrangement=sigma_R,
        converged=converged, n_iter=it + 1,
    )


def load_models(path: Path | str) -> dict[str, RMFModel]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    models = {}
    for key, m in data["models"].items():
        couplings = {meson: RMFCouplings(a=p["a"], b=p["b"], c=p["c"], d=p.get("d", 0.0))
                     for meson, p in m["coupling_params"].items()}
        models[key] = RMFModel(
            key=key, family=m["family"], functional_form=m["functional_form"],
            meson_masses=dict(m["meson_masses_MeV"]),
            rho_sat=m["rho_sat_fm3"], couplings=couplings,
            saturation_targets=m.get("saturation_targets", {}),
            reference=m.get("reference", ""),
        )
    return models


def energy_per_baryon(r: HomogeneousResult) -> float:
    return r.energy_density / r.rho_B if r.rho_B > 0 else 0.0


def symmetry_energy(model: RMFModel, rho_B: float) -> float:
    """S(rho) via parabolic approx: S = (1/2) d^2(E/A)/d(delta^2)|_{delta=0}."""
    delta = 0.02
    rp = solve_homogeneous(model, rho_B, rho_B * (1 + delta) / 2)
    rm = solve_homogeneous(model, rho_B, rho_B * (1 - delta) / 2)
    return (energy_per_baryon(rm) - energy_per_baryon(rp)) / (2 * delta)
