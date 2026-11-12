"""RMF homogeneous-matter solver for the EOS390 CLD reconstruction.

Implements the density-dependent relativistic mean-field (DD-RMF) energy
density for homogeneous asymmetric nuclear matter. Multiple coupling
functional forms are supported, each with its EXACT published form:

  1. DD2_TYYPEL_RATIONAL (Typel 2010, arXiv:0908.2344, eq. 75-77):
       isoscalar (omega, sigma): rational function
         f_i(x) = a_i * [1 + b_i*(x + d_i)^2] / [1 + c_i*(x + d_i)^2],  x = n/n_sat
         g_i(n) = Gamma_i(n_sat) * f_i(x)
       isovector (rho): exponential
         f_rho(x) = exp[-a_rho * (x - 1)]
       This is the DD2 functional form. The DD-ME2, TW99, FSU2, FSU2H,
       NL3wrL55, TM1e models use related (but check-per-model) forms.

  2. SPG_GOGELEIN (Gögelein 2008, eq. 40 + Table I):
       g_kappa(rho_B) = a_kappa + (b_kappa + d_kappa*x^3)*exp(-c_kappa*x), x = rho_B
       4 mesons: sigma, omega, delta, rho. SPG_M1..M5 use this form.

UNIT CONVENTION: natural units (hbar=c=1), MeV internally, MeV/fm^3 output.
Densities input in fm^-3, converted to MeV^3 by * HBARC^3 internally.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

HBARC = 197.3269804
M_N = 939.56536
M_P = 938.27203
M_E = 0.511
M_MU = 105.658
PI = math.pi
# Gögelein 2008 PRC 77 025802 eq. 40: SPG functional form uses x = rho_B / rho_0
# with the FIXED reference saturation density rho_0 = 0.16 fm^-3 (NOT the model's
# actual rho_sat, which is 0.178 fm^-3 after the omega correction). This is the
# density at which the DBHF fit was anchored.
SPG_RHO_REF_FM3 = 0.16


def _kF_MeV(rho_MeV3: float) -> float:
    if rho_MeV3 <= 0:
        return 0.0
    return (3.0 * PI * PI * rho_MeV3) ** (1.0 / 3.0)


def _scalar_density_MeV3(kF: float, Mstar: float) -> float:
    if kF <= 0:
        return 0.0
    E = math.sqrt(kF * kF + Mstar * Mstar)
    return (Mstar / (2.0 * PI * PI)) * (kF * E - Mstar * Mstar * math.log((kF + E) / Mstar))


def _kinetic_MeV4(kF: float, Mstar: float, g: int) -> float:
    if kF <= 0:
        return 0.0
    E = math.sqrt(kF * kF + Mstar * Mstar)
    I = 0.125 * (kF * (2 * kF * kF + Mstar * Mstar) * E - Mstar ** 4 * math.log((kF + E) / Mstar))
    return (g / (2.0 * PI * PI)) * I


def _solve_sigma_nonlinear(rhs: float, m_s: float, g2: float, g3: float,
                             sigma_init: float = 0.0,
                             tol: float = 1e-12, max_iter: int = 200) -> float:
    """Solve the NL3_WALECKA_NONLINEAR sigma field equation:
        m_s^2 sigma = rhs + g_2 sigma^2 + g_3 sigma^3
        where rhs = g_sigma * rho_s, and g_2, g_3 are the nonlinear self-couplings.

    Convention (Lalazissis 1997 PRC 55 540, NL3 standard):
        L ⊃ -U(sigma) = -(1/2) m_s^2 sigma^2 - (g_2/3) sigma^3 - (g_3/4) sigma^4
        Field eq (dU/dsigma = source):  m_s^2 sigma + g_2 sigma^2 + g_3 sigma^3 = g_sigma rho_s
        i.e. m_s^2 sigma = rhs - g_2 sigma^2 - g_3 sigma^3 = rhs + g_2' sigma^2 + g_3' sigma^3
        where g_2' = -g_2, g_3' = -g_3 (the "self-coupling in field eq" coefficients).

        In this solver we use the convention: g_2, g_3 passed in are the LAGRANGIAN
        coefficients (g_2 = -10.4310, g_3 = -28.8856 for NL3). The field equation
        becomes: m_s^2 sigma = rhs - g_2 sigma^2 - g_3 sigma^3
        => f(sigma) = m_s^2 sigma + g_2 sigma^2 + g_3 sigma^3 - rhs = 0  (lagrangian form)
        which for NL3 (g_2<0, g_3<0) means the nonlinear terms REDUCE the rhs needed,
        keeping sigma physical (~40 MeV at saturation).

        Energy density (Hamiltonian) uses the LAGRANGIAN coefficients directly:
        E_sigma = (1/2) m_s^2 sigma^2 + (g_2/3) sigma^3 + (g_3/4) sigma^4
        For NL3 g_2,g_3<0, the cubic and quartic terms contribute NEGATIVELY (binding).

    f(sigma) = m_s^2 sigma + g_2 sigma^2 + g_3 sigma^3 - rhs = 0
    f'(sigma) = m_s^2 + 2 g_2 sigma + 3 g_3 sigma^2
    Returns sigma >= 0 (positive mean field).

    Robust hybrid: Newton with damping + bisection fallback, since the cubic
    can make f'(sigma) change sign for large sigma.
    """
    if rhs <= 0:
        return 0.0
    if g2 == 0.0 and g3 == 0.0:
        return rhs / (m_s * m_s)
    sigma_lin = rhs / (m_s * m_s)
    # For NL3 (g_2,g_3<0) the nonlinear terms reduce sigma, so the physical
    # root is in (0, sigma_lin). Start Newton from a modest fraction.
    sigma = sigma_init if (0 < sigma_init < sigma_lin) else sigma_lin * 0.3
    for it in range(max_iter):
        f = m_s * m_s * sigma + g2 * sigma * sigma + g3 * sigma ** 3 - rhs
        fp = m_s * m_s + 2.0 * g2 * sigma + 3.0 * g3 * sigma * sigma
        if abs(f) < tol * max(1.0, abs(rhs)):
            break
        if fp > 1.0:
            ds = -f / fp
            max_ds = 0.3 * sigma + 1.0
            if abs(ds) > max_ds:
                ds = ds / abs(ds) * max_ds
            sigma_new = sigma + ds
            if sigma_new < 0 or sigma_new > sigma_lin:
                sigma_new = max(0.1 * sigma, min(sigma_lin, sigma_new))
        else:
            # bisection fallback on [0, sigma_lin]
            lo, hi = 0.0, sigma_lin
            f_lo = -rhs
            for trial in [0.5, 0.3, 0.2, 0.15, 0.1, 0.05]:
                hi = sigma_lin * trial
                f_hi = m_s * m_s * hi + g2 * hi * hi + g3 * hi ** 3 - rhs
                if f_hi >= 0:
                    break
            if f_hi < 0:
                sigma = sigma_lin * 0.3
                break
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                f_mid = m_s * m_s * mid + g2 * mid * mid + g3 * mid ** 3 - rhs
                if abs(f_mid) < tol * max(1.0, abs(rhs)):
                    sigma = mid
                    break
                if (f_lo < 0) != (f_mid < 0):
                    hi = mid
                else:
                    lo = mid
                    f_lo = f_mid
            sigma = 0.5 * (lo + hi)
            break
        if abs(sigma_new - sigma) < tol * max(1.0, abs(sigma)):
            sigma = sigma_new
            break
        sigma = sigma_new
    return sigma


@dataclass
class RMFCouplings:
    """Generic coupling parameter set. The interpretation of (a,b,c,d,Gamma)
    depends on the functional form of the parent model.

    For NL3_WALECKA_NONLINEAR: a = g_sigma (constant), Gamma_sat = g_i(rho) for
    ALL densities (no density dependence). The nonlinear sigma self-interaction
    coefficients g_2 (cubic, MeV) and g_3 (quartic, dimensionless) and the omega-rho
    cross-coupling Lambda_v are stored in the RMFModel, not here, because they are
    model-level (not per-meson) parameters. b, c, d are unused for Walecka."""
    a: float
    b: float = 0.0
    c: float = 0.0
    d: float = 0.0
    Gamma_sat: float = 0.0  # Gamma_i(n_sat), the saturation-point coupling

    def value_rational(self, x: float) -> float:
        """DD2 isoscalar: f = a*(1 + b*(x+d)^2)/(1 + c*(x+d)^2); g = Gamma*f."""
        num = 1.0 + self.b * (x + self.d) ** 2
        den = 1.0 + self.c * (x + self.d) ** 2
        return self.Gamma_sat * self.a * num / den

    def value_exp(self, x: float) -> float:
        """DD2 isovector (rho): f = exp(-a*(x-1)); g = Gamma*f."""
        return self.Gamma_sat * math.exp(-self.a * (x - 1.0))

    def value_gogelein(self, x: float) -> float:
        """SPG Gögelein: g = a + (b + d*x^3)*exp(-c*x)."""
        return self.a + (self.b + self.d * x ** 3) * math.exp(-self.c * x)

    def value_walecka(self, x: float) -> float:
        """NL3_WALECKA_NONLINEAR: constant coupling (no density dependence).
        For the meson-nucleon coupling g_i, returns Gamma_sat. x is unused.
        The sigma self-interaction g_2, g_3 and omega-rho Lambda_v are handled
        in solve_homogeneous via model-level fields."""
        return self.Gamma_sat

    def dvalue_rational_dx(self, x: float) -> float:
        """d/dx of a*(1+b(x+d)^2)/(1+c(x+d)^2)."""
        u = x + self.d
        num = 1.0 + self.b * u * u
        den = 1.0 + self.c * u * u
        dnum = 2.0 * self.b * u
        dden = 2.0 * self.c * u
        f = self.a * num / den
        df = self.a * (dnum * den - num * dden) / (den * den)
        return self.Gamma_sat * df

    def dvalue_exp_dx(self, x: float) -> float:
        f = math.exp(-self.a * (x - 1.0))
        return self.Gamma_sat * (-self.a) * f

    def dvalue_gogelein_dx(self, x: float) -> float:
        return (3.0 * self.d * x * x - self.c * (self.b + self.d * x ** 3)) * math.exp(-self.c * x)

    def dvalue_walecka_dx(self, x: float) -> float:
        """NL3_WALECKA_NONLINEAR: density-independent coupling => derivative = 0.
        Rearrangement term vanishes for Walecka models."""
        return 0.0


@dataclass
class RMFModel:
    key: str
    family: str
    functional_form: str   # "DD2_TYYPEL_RATIONAL", "SPG_GOGELEIN", or "NL3_WALECKA_NONLINEAR"
    meson_masses: dict
    rho_sat: float
    couplings: dict        # meson -> RMFCouplings
    saturation_targets: dict
    reference: str = ""
    # which mesons use the rational (isoscalar) vs exponential (isovector) form:
    isoscalar_mesons: tuple = ("omega", "sigma")
    isovector_mesons: tuple = ("rho",)
    # NL3_WALECKA_NONLINEAR nonlinear terms (zero for DD2/SPG models):
    g2: float = 0.0   # sigma cubic self-coupling, convention: -g_2*sigma^3/3 in L;
                      # field eq: m_sigma^2 sigma + g_2 sigma^2 + g_3 sigma^3 = g_sigma rho_s.
                      # Units: MeV (sigma is in MeV), g_2 in MeV, g_3 dimensionless.
                      # NL3 canonical: g_2 = -10.4310 MeV, g_3 = -28.8856.
    g3: float = 0.0   # sigma quartic self-coupling, dimensionless.
    Lambda_v: float = 0.0  # omega-rho cross-coupling (FSU/IU convention):
                           # L += -Lambda_v * (g_omega^2 omega^2)(g_rho^2 rho_3^2);
                           # modifies field eqs:
                           #   m_omega^2 omega_0 (1 + Lambda_v g_rho^2 rho_3^2 / m_rho^2) = g_omega rho_B
                           #   m_rho^2 rho_03 (1 + Lambda_v g_omega^2 omega_0^2 / m_omega^2) = g_rho rho_3
                           # NL3/TM1 have Lambda_v = 0; FSU2/IU-FSU/FSUGarnet have Lambda_v != 0.

    @property
    def rho_sat_MeV3(self) -> float:
        return self.rho_sat * HBARC ** 3

    def x_of(self, rho_B_MeV3: float) -> float:
        if self.functional_form == "SPG_GOGELEIN":
            # Gögelein eq. 40: x = rho_B / rho_0 with rho_0 = 0.16 fm^-3 (fixed)
            return rho_B_MeV3 / HBARC ** 3 / SPG_RHO_REF_FM3
        return rho_B_MeV3 / self.rho_sat_MeV3

    def g(self, meson: str, rho_B_MeV3: float) -> float:
        x = self.x_of(rho_B_MeV3)
        c = self.couplings[meson]
        if self.functional_form == "SPG_GOGELEIN":
            return c.value_gogelein(x)
        if self.functional_form == "NL3_WALECKA_NONLINEAR":
            return c.value_walecka(x)
        # DD2_TYYPEL_RATIONAL
        if meson in self.isovector_mesons:
            return c.value_exp(x)
        return c.value_rational(x)

    def dg_drho_fm(self, meson: str, rho_B_MeV3: float) -> float:
        """dg/d(rho_B[fm^-3])."""
        x = self.x_of(rho_B_MeV3)
        c = self.couplings[meson]
        if self.functional_form == "SPG_GOGELEIN":
            # x = rho_B / rho_0 => dg/drho_B = (dg/dx) * (1/rho_0)
            return c.dvalue_gogelein_dx(x) / SPG_RHO_REF_FM3
        if self.functional_form == "NL3_WALECKA_NONLINEAR":
            # constant coupling => dg/drho = 0 (no rearrangement term)
            return 0.0
        dgdx = c.dvalue_exp_dx(x) if meson in self.isovector_mesons else c.dvalue_rational_dx(x)
        return dgdx / self.rho_sat  # x = rho/rho_sat

    def includes_delta(self) -> bool:
        return "delta" in self.couplings

    def is_walecka(self) -> bool:
        return self.functional_form == "NL3_WALECKA_NONLINEAR"


@dataclass
class HomogeneousResult:
    rho_B: float; rho_3: float; rho_p: float; rho_n: float
    sigma0: float; omega0: float; rho03: float; delta0: float
    Mstar_n: float; Mstar_p: float
    energy_density: float; binding_energy: float; pressure: float
    mu_n: float; mu_p: float; rearrangement: float
    converged: bool; n_iter: int


def solve_homogeneous(model: RMFModel, rho_B: float, rho_p: float,
                      tol: float = 1e-9, max_iter: int = 1000) -> HomogeneousResult:
    """Self-consistent DD-RMF for homogeneous matter at (rho_B, rho_p) in fm^-3."""
    if rho_B <= 0:
        return HomogeneousResult(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, True, 0)
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
    sigma0_prev = 0.0
    is_walecka = model.is_walecka()
    g2_w = getattr(model, "g2", 0.0)
    g3_w = getattr(model, "g3", 0.0)
    Lam_v = getattr(model, "Lambda_v", 0.0)
    for it in range(max_iter):
        kF_n = _kF_MeV(rho_n_M)
        kF_p = _kF_MeV(rho_p_M)
        rs_n = _scalar_density_MeV3(kF_n, Mstar_n)
        rs_p = _scalar_density_MeV3(kF_p, Mstar_p)
        rs_total = rs_n + rs_p
        rs_3 = rs_p - rs_n
        if is_walecka:
            # NL3_WALECKA_NONLINEAR: nonlinear sigma field eq + omega-rho cross-coupling
            # sigma: m_s^2 sigma + g2 sigma^2 + g3 sigma^3 = g_s * rho_s
            sigma0 = _solve_sigma_nonlinear(g_s * rs_total, m_s, g2_w, g3_w,
                                            sigma_init=sigma0_prev)
            sigma0_prev = sigma0
            # omega-rho cross-coupling: m_w^2 w0 (1 + Lam_v g_r^2 rho_3^2/m_r^2) = g_w rho_B
            # rho_3 here in MeV^3 = (rho_p - rho_n)
            rho3_M = rho_p_M - rho_n_M
            factor_w = 1.0 + Lam_v * (g_r * g_r) * (rho3_M * rho3_M) / (m_r * m_r)
            omega0 = g_w * rho_B_M / (m_w * m_w * factor_w)
            # rho field: m_r^2 b3 (1 + Lam_v g_w^2 w0^2/m_w^2) = g_r * rho_3
            factor_r = 1.0 + Lam_v * (g_w * g_w) * (omega0 * omega0) / (m_w * m_w)
            rho03 = g_r * rho3_M / (m_r * m_r * factor_r)
            delta0 = 0.0  # NL3/TM1/FSU2 do not include delta meson
        else:
            sigma0 = g_s * rs_total / (m_s ** 2)
            omega0 = g_w * rho_B_M / (m_w ** 2)
            rho03 = g_r * (rho_p_M - rho_n_M) / (m_r ** 2)
            delta0 = g_d * rs_3 / (m_d ** 2) if has_delta else 0.0
        # Gögelein eq. 5: m*_i = M + g_sigma*sigma + g_delta*delta*tau_3, with
        # delta = -g_delta*rho_3s/m_delta^2 (KG eq. 8 has minus sign). Our convention
        # sigma0>0, delta0 = g_d*rs_3/m_d^2 (sign follows rs_3). To match Gögelein's
        # m*_n = M + g_sigma*sigma - g_delta*delta (neutron lighter in neutron-rich
        # matter, delta>0 there), we use M*_n = M - g_s*sigma0 + g_d*delta0 and
        # M*_p = M - g_s*sigma0 - g_d*delta0 (delta0<0 in neutron-rich -> M*_n drops).
        Mstar_n_new = M_N - g_s * sigma0 + (g_d * delta0 if has_delta else 0.0)
        Mstar_p_new = M_P - g_s * sigma0 - (g_d * delta0 if has_delta else 0.0)
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

    kF_n = _kF_MeV(rho_n_M)
    kF_p = _kF_MeV(rho_p_M)
    rs_n = _scalar_density_MeV3(kF_n, Mstar_n)
    rs_p = _scalar_density_MeV3(kF_p, Mstar_p)
    rs_total = rs_n + rs_p
    rs_3 = rs_p - rs_n
    if is_walecka:
        sigma0 = _solve_sigma_nonlinear(g_s * rs_total, m_s, g2_w, g3_w,
                                        sigma_init=sigma0_prev)
        rho3_M = rho_p_M - rho_n_M
        factor_w = 1.0 + Lam_v * (g_r * g_r) * (rho3_M * rho3_M) / (m_r * m_r)
        omega0 = g_w * rho_B_M / (m_w * m_w * factor_w)
        factor_r = 1.0 + Lam_v * (g_w * g_w) * (omega0 * omega0) / (m_w * m_w)
        rho03 = g_r * rho3_M / (m_r * m_r * factor_r)
        delta0 = 0.0
    else:
        sigma0 = g_s * rs_total / (m_s ** 2)
        omega0 = g_w * rho_B_M / (m_w ** 2)
        rho03 = g_r * (rho_p_M - rho_n_M) / (m_r ** 2)
        delta0 = g_d * rs_3 / (m_d ** 2) if has_delta else 0.0

    # Energy density (MeV^4), Scurto 2024 eq. 26 (Esat) generalized to asymmetry:
    #   E = sum_i Ekin,i + (1/2)(m_s^2 sigma0^2 + m_w^2 w0^2 + m_r^2 b3^2)
    # This is the no-sea RMF form where ALL meson field self-energies are POSITIVE
    # and there is NO separate vector-interaction term (the omega/rho vector
    # coupling is absorbed: w0 = g_w rho/m_w^2 so 1/2 m_w^2 w0^2 = g_w^2 rho^2/2m_w^2
    # already encodes the repulsion; the scalar 1/2 m_s^2 sigma^2 with sigma = g_s rho_s/m_s^2
    # encodes the binding via M* = M - g_s sigma). The kinetic Ekin,i is the FULL
    # integral including M* rest mass. Binding E/A = E/rho - M (Scurto common M).
    # For SPG with delta meson, add +1/2 m_d^2 delta0^2 to the field terms.
    Ekin_n = _kinetic_MeV4(kF_n, Mstar_n, 2)
    Ekin_p = _kinetic_MeV4(kF_p, Mstar_p, 2)
    Ekin = Ekin_n + Ekin_p

    Efield = (0.5 * m_s ** 2 * sigma0 ** 2
              + 0.5 * m_w ** 2 * omega0 ** 2
              + 0.5 * m_r ** 2 * rho03 ** 2)
    if has_delta:
        Efield += 0.5 * m_d ** 2 * delta0 ** 2

    # NL3_WALECKA_NONLINEAR nonlinear self-interaction energy terms.
    # Lagrangian L ⊃ -U(sigma) = -(g_2/3) sigma^3 - (g_3/4) sigma^4 (Lalazissis 1997).
    # The contribution to the energy density is +U(sigma) (the potential energy
    # stored in the sigma field), i.e. +(g_2/3) sigma^3 + (g_3/4) sigma^4.
    # For NL3 with g_2<0, g_3<0, sigma>0, these terms are negative (binding contribution).
    # The omega-rho cross-coupling L ⊃ -Lam_v g_w^2 g_r^2 w0^2 b3^2 contributes
    # +Lam_v g_w^2 g_r^2 w0^2 b3^2 to the energy density (Todd-Rutel 2005 convention).
    if is_walecka:
        Efield += (g2_w / 3.0) * sigma0 ** 3 + (g3_w / 4.0) * sigma0 ** 4
        if Lam_v != 0.0:
            Efield += Lam_v * (g_w * g_w) * (g_r * g_r) * (omega0 * omega0) * (rho03 * rho03)

    E_MeV4 = Ekin + Efield
    E_MeV_fm3 = E_MeV4 / HBARC ** 3

    # binding energy density (Typel eq. 69 / Scurto common-M convention): subtract
    # free rest masses. Use individual m_n, m_p (Typel DD2) for best reproduction.
    E_bind_MeV4 = E_MeV4 - M_N * rho_n_M - M_P * rho_p_M
    E_bind_MeV_fm3 = E_bind_MeV4 / HBARC ** 3

    # rearrangement term Sigma^R (MeV), Scurto 2024 eq. 14:
    #   Sigma^R = (dg_w/drho) w0 rho + (1/2)(dg_r/drho) b3 rho_3 - (dg_s/drho) phi0 rho_s
    # (the sigma term is MINUS because the scalar field lowers M*; the rho term
    # has a 1/2 from the isospin convention). For SPG delta models add the
    # analogous delta term with the isovector scalar 1/2.
    rho_n = rho_B - rho_p
    dgs = model.dg_drho_fm("sigma", rho_B_M)
    dgw = model.dg_drho_fm("omega", rho_B_M)
    dgr = model.dg_drho_fm("rho", rho_B_M)
    sigma_R = (dgw * omega0 * rho_B
               + 0.5 * dgr * rho03 * (rho_p - rho_n)
               - dgs * sigma0 * (rs_n + rs_p) / HBARC ** 3)
    if has_delta:
        dgd = model.dg_drho_fm("delta", rho_B_M)
        sigma_R += 0.5 * dgd * delta0 * (rs_p - rs_n) / HBARC ** 3

    En_n = math.sqrt(kF_n ** 2 + Mstar_n ** 2)
    En_p = math.sqrt(kF_p ** 2 + Mstar_p ** 2)
    mu_n = En_n + g_w * omega0 - 0.5 * g_r * rho03 - sigma_R
    mu_p = En_p + g_w * omega0 + 0.5 * g_r * rho03 - sigma_R

    # Pressure (Typel eq. 66): the kinetic pressure + rearrangement contribution.
    # In the no-sea mean-field at T=0, the pressure can be computed from the
    # binding energy density via the thermodynamic identity using the chemical
    # potentials. But for consistency with the binding-energy convention, we use
    # P = sum_i mu_i rho_i - E_total, where E_total is the FULL energy density
    # (including rest mass) and mu_i are the quasiparticle chemical potentials
    # (which include the vector self-energy and rearrangement). At saturation
    # this gives P ~ 0 by construction.
    P_MeV4 = (mu_n * rho_n_M + mu_p * rho_p_M) - E_MeV4
    P_MeV_fm3 = P_MeV4 / HBARC ** 3

    return HomogeneousResult(
        rho_B=rho_B, rho_3=rho_p - rho_n, rho_p=rho_p, rho_n=rho_n,
        sigma0=sigma0, omega0=omega0, rho03=rho03, delta0=delta0,
        Mstar_n=Mstar_n, Mstar_p=Mstar_p,
        energy_density=E_MeV_fm3, binding_energy=E_bind_MeV_fm3, pressure=P_MeV_fm3,
        mu_n=mu_n, mu_p=mu_p, rearrangement=sigma_R,
        converged=converged, n_iter=it + 1,
    )


def load_models(path: Path | str) -> dict[str, RMFModel]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    models = {}
    for key, m in data["models"].items():
        couplings = {}
        for meson, p in m["coupling_params"].items():
            couplings[meson] = RMFCouplings(
                a=p["a"], b=p.get("b", 0.0), c=p.get("c", 0.0),
                d=p.get("d", 0.0), Gamma_sat=p.get("Gamma_sat", 0.0))
        models[key] = RMFModel(
            key=key, family=m["family"], functional_form=m["functional_form"],
            meson_masses=dict(m["meson_masses_MeV"]),
            rho_sat=m["rho_sat_fm3"], couplings=couplings,
            saturation_targets=m.get("saturation_targets", {}),
            reference=m.get("reference", ""),
            isoscalar_mesons=tuple(m.get("isoscalar_mesons", ("omega", "sigma"))),
            isovector_mesons=tuple(m.get("isovector_mesons", ("rho",))),
            g2=m.get("g2", 0.0),
            g3=m.get("g3", 0.0),
            Lambda_v=m.get("Lambda_v", 0.0),
        )
    return models


def energy_per_baryon(r: HomogeneousResult) -> float:
    """Binding energy per baryon (MeV), i.e. eps/n - M (free rest mass subtracted).
    This is the E/A of nuclear physics, target ~ -16 MeV at saturation."""
    return r.binding_energy / r.rho_B if r.rho_B > 0 else 0.0


def symmetry_energy(model: RMFModel, rho_B: float) -> float:
    """Symmetry energy S(rho) via the parabolic approximation:
    S = (1/2) d^2(E/A)/d(delta^2)|_{delta=0}  where delta = (rho_n - rho_p)/rho_B.
    Numerically: S = (E(+d) - 2*E(0) + E(-d)) / (2*d^2) with E(delta) the E/A at
    asymmetry delta. Since E is symmetric in delta (E(+d)=E(-d)), this simplifies
    to S = (E(d) - E(0)) / d^2. Using d=0.05 gives a stable estimate.
    """
    d = 0.05
    r0 = solve_homogeneous(model, rho_B, rho_B / 2)              # delta=0 (symmetric)
    rp = solve_homogeneous(model, rho_B, rho_B * (1 + d) / 2)    # delta=+d
    rm = solve_homogeneous(model, rho_B, rho_B * (1 - d) / 2)    # delta=-d
    e0 = energy_per_baryon(r0)
    return (energy_per_baryon(rp) - 2 * e0 + energy_per_baryon(rm)) / (2 * d ** 2)
