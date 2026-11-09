"""Thermodynamically consistent low-density continuation of a projected cold EOS.

The extension is a diagnostic surface closure, not a sourced crust EOS.  It
continues P(epsilon) below the first projected knot and obtains n_b from the
cold first law, then defines mu_b=(epsilon+P)/n_b.  This preserves
n_b*d(mu_b)=dP analytically and keeps every normalized table state positive.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

MEV_FM3_TO_PA = 1.602176634e32
MEV_TO_J = 1.602176634e-13
FM3_TO_M3 = 1.0e45
NORMALIZED_FIELDS = (
    "source_row", "i_t", "i_nb", "i_yq", "n_b_fm3", "n_b_m3",
    "pressure_mev_fm3", "pressure_pa", "energy_mev_fm3", "energy_j_m3",
    "mu_b_mev", "mu_b_j", "q1_p_over_nb_mev", "q2_entropy_per_baryon",
    "q3_mu_b_over_mn_minus_1", "q4_mu_q_over_mn", "q5_mu_l_over_mn",
    "q6_f_over_nb_mn_minus_1", "q7_e_over_nb_mn_minus_1",
)
_GL_X, _GL_W = np.polynomial.legendre.leggauss(96)


def _number(row: Mapping[str, object], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"projected EOS row is missing numeric field {field!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"projected EOS field {field!r} is nonfinite")
    return value


@dataclass(frozen=True)
class ThermodynamicSurfaceContinuation:
    """A C1 curvature family matched to the first raw Hermite P(mu) interval.

    Let x=epsilon/epsilon_cut.  The pressure family is

        P(x)=P_cut*x**gamma_s*exp(alpha*(1-x)**2).

    Since the curvature factor has zero first derivative at x=1, all alpha
    values share the same cutoff pressure and dP/depsilon.  The baryon density
    is integrated from d epsilon=(epsilon+P) d n_b/n_b.
    """

    epsilon_cut_mev_fm3: float
    pressure_cut_mev_fm3: float
    number_density_cut_fm3: float
    chemical_potential_cut_mev: float
    cs2_cut: float
    gamma_s: float
    alpha: float
    family: str
    neutron_mass_mev: float
    auxiliary_cutoff: Mapping[str, str]

    @classmethod
    def from_projected_rows(
        cls,
        rows: Sequence[Mapping[str, object]],
        *,
        alpha: float = 0.0,
        family: str = "power_law",
    ) -> "ThermodynamicSurfaceContinuation":
        if len(rows) < 2:
            raise ValueError("at least two raw projected EOS knots are required")
        if not math.isfinite(alpha):
            raise ValueError("alpha must be finite")
        if family not in {"power_law", "exp_quadratic", "poly_quadratic"}:
            raise ValueError(f"unsupported continuation family: {family}")
        if family == "power_law" and alpha != 0.0:
            raise ValueError("power_law family has no curvature parameter; alpha must be zero")
        row0, row1 = rows[0], rows[1]
        mu0 = _number(row0, "mu_b_mev")
        mu1 = _number(row1, "mu_b_mev")
        p0 = _number(row0, "pressure_mev_fm3")
        p1 = _number(row1, "pressure_mev_fm3")
        n0 = _number(row0, "n_b_fm3")
        n1 = _number(row1, "n_b_fm3")
        eps0 = _number(row0, "energy_mev_fm3")
        if not (mu1 > mu0 and p1 > p0 and n1 > n0 and eps0 > 0.0 and p0 > 0.0 and n0 > 0.0):
            raise ValueError("the first two EOS knots must be positive and strictly increasing")

        # Right-sided derivative of the exact native cubic-Hermite P(mu) cell.
        h = mu1 - mu0
        a2 = 3.0 * (p1 - p0) - h * (2.0 * n0 + n1)
        dn_dmu = 2.0 * a2 / (h * h)
        cs2 = n0 / (mu0 * dn_dmu) if dn_dmu > 0.0 else math.nan
        gamma = eps0 / p0 * cs2
        mn = mu0 / (1.0 + _number(row0, "q3_mu_b_over_mn_minus_1"))
        closure = (mu0 * n0) / (eps0 + p0) - 1.0
        if not math.isfinite(cs2) or not (0.0 < cs2 < 1.0):
            raise ValueError(f"cutoff Hermite sound speed is not causal and positive: {cs2}")
        # Gamma_s>6/5 is the finite-radius Newtonian-polytrope condition for
        # the asymptotic P~epsilon^Gamma_s tail.
        if not math.isfinite(gamma) or gamma <= 1.2:
            raise ValueError(f"surface exponent does not give a finite-radius tail: {gamma}")
        if not math.isfinite(mn) or mn <= 0.0 or abs(closure) > 2.0e-6:
            raise ValueError("cutoff state fails mu_b*n_b=epsilon+P thermodynamic closure")

        model = cls(
            epsilon_cut_mev_fm3=eps0,
            pressure_cut_mev_fm3=p0,
            number_density_cut_fm3=n0,
            chemical_potential_cut_mev=mu0,
            cs2_cut=cs2,
            gamma_s=gamma,
            alpha=float(alpha),
            family=family,
            neutron_mass_mev=mn,
            auxiliary_cutoff={
                key: str(row0[key])
                for key in (
                    "i_t", "i_yq", "q2_entropy_per_baryon",
                    "q4_mu_q_over_mn", "q5_mu_l_over_mn",
                )
            },
        )
        model._validate_shape(min_energy_ratio=1.0e-12)
        return model

    @classmethod
    def from_projected_csv(
        cls, path: str | Path, *, alpha: float = 0.0, family: str = "power_law"
    ) -> "ThermodynamicSurfaceContinuation":
        with Path(path).open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        return cls.from_projected_rows(rows, alpha=alpha, family=family)

    def _shape_factor(self, x: np.ndarray) -> np.ndarray:
        u2 = np.square(1.0 - x)
        if self.family == "power_law":
            return np.ones_like(x)
        if self.family == "exp_quadratic":
            return np.exp(self.alpha * u2)
        return 1.0 + self.alpha * u2

    def _slope_factor(self, x: np.ndarray) -> np.ndarray:
        if self.family == "power_law":
            return np.full_like(x, self.gamma_s)
        if self.family == "exp_quadratic":
            return self.gamma_s - 2.0 * self.alpha * x * (1.0 - x)
        factor = self._shape_factor(x)
        return self.gamma_s - 2.0 * self.alpha * x * (1.0 - x) / factor

    def _pressure_ratio(self, x: np.ndarray) -> np.ndarray:
        return (self.pressure_cut_mev_fm3 / self.epsilon_cut_mev_fm3) * np.power(
            x, self.gamma_s - 1.0
        ) * self._shape_factor(x)

    def _log_chemical_potential(self, x: np.ndarray) -> np.ndarray:
        """Integrate the cold first law for log(mu_b) in log energy density.

        From dP=n_b dmu_b and c_s^2=dP/depsilon,

            d ln(mu_b) / d ln(epsilon) = c_s^2 / (1 + P/epsilon).

        Integrating this small logarithmic increment is better conditioned in
        the dilute tail than integrating ln(n_b) and subtracting it from
        ln(epsilon+P) to reconstruct an almost constant chemical potential.
        """
        u = np.log(x)
        t = 0.5 * u[:, None] * (_GL_X[None, :] + 1.0)
        xq = np.exp(t)
        ratio = self._pressure_ratio(xq)
        cs2 = ratio * self._slope_factor(xq)
        integrand = cs2 / (1.0 + ratio)
        integral = 0.5 * u * np.sum(integrand * _GL_W[None, :], axis=1)
        return math.log(self.chemical_potential_cut_mev) + integral

    def state_many(self, epsilon_mev_fm3: Iterable[float]) -> dict[str, np.ndarray]:
        eps = np.asarray(tuple(epsilon_mev_fm3), dtype=float)
        if eps.ndim != 1 or eps.size == 0 or not np.all(np.isfinite(eps)):
            raise ValueError("epsilon query must be a nonempty finite one-dimensional sequence")
        if np.any(eps <= 0.0) or np.any(eps > self.epsilon_cut_mev_fm3 * (1.0 + 8*np.finfo(float).eps)):
            raise ValueError("epsilon query must lie in (0, epsilon_cut]")
        x = np.minimum(eps / self.epsilon_cut_mev_fm3, 1.0)
        pressure = self.pressure_cut_mev_fm3 * np.power(x, self.gamma_s) * self._shape_factor(x)
        slope_factor = self._slope_factor(x)
        cs2 = pressure / eps * slope_factor
        chemical_potential = np.exp(self._log_chemical_potential(x))
        number_density = (eps + pressure) / chemical_potential
        arrays = (pressure, number_density, chemical_potential, cs2)
        if not all(np.all(np.isfinite(a)) for a in arrays):
            raise ArithmeticError("surface continuation produced a nonfinite thermodynamic state")
        if np.any(pressure <= 0.0) or np.any(number_density <= 0.0) or np.any(chemical_potential <= 0.0):
            raise ArithmeticError("surface continuation produced a nonpositive state")
        if np.any(cs2 <= 0.0) or np.any(cs2 >= 1.0):
            raise ArithmeticError("surface continuation is non-hyperbolic or acausal")
        return {
            "epsilon_mev_fm3": eps,
            "pressure_mev_fm3": pressure,
            "number_density_fm3": number_density,
            "chemical_potential_mev": chemical_potential,
            "cs2": cs2,
        }

    def state_at_energy(self, epsilon_mev_fm3: float) -> dict[str, float]:
        state = self.state_many([float(epsilon_mev_fm3)])
        return {key: float(value[0]) for key, value in state.items()}

    def _validate_shape(self, *, min_energy_ratio: float) -> None:
        if not (math.isfinite(min_energy_ratio) and 0.0 < min_energy_ratio < 1.0):
            raise ValueError("min_energy_ratio must lie strictly between zero and one")
        x = np.geomspace(min_energy_ratio, 1.0, 512)
        try:
            state = self.state_many(self.epsilon_cut_mev_fm3 * x)
        except ArithmeticError as exc:
            raise ValueError(f"surface continuation shape is thermodynamically inadmissible: {exc}") from exc
        for key in ("pressure_mev_fm3", "number_density_fm3", "chemical_potential_mev"):
            values = state[key]
            if np.any(np.diff(values) <= 0.0):
                raise ValueError(f"surface continuation is not strictly monotone in {key}")

    def diagnostics(self, *, min_energy_ratio: float = 1.0e-12, points: int = 512) -> dict:
        if points < 16:
            raise ValueError("diagnostic grid requires at least 16 points")
        x = np.geomspace(min_energy_ratio, 1.0, points)
        state = self.state_many(self.epsilon_cut_mev_fm3 * x)
        loge = np.log(state["epsilon_mev_fm3"])
        logn = np.log(state["number_density_fm3"])
        first_law_lhs = np.gradient(logn, loge, edge_order=2)
        first_law_rhs = state["epsilon_mev_fm3"] / (
            state["epsilon_mev_fm3"] + state["pressure_mev_fm3"]
        )
        return {
            "model_id": f"c1_{self.family}_alpha_{self.alpha:+.3f}",
            "family": self.family,
            "alpha": self.alpha,
            "gamma_s": self.gamma_s,
            "gamma_eq_cut": (self.epsilon_cut_mev_fm3 + self.pressure_cut_mev_fm3)
            / self.pressure_cut_mev_fm3 * self.cs2_cut,
            "K_mev_fm3_power": self.pressure_cut_mev_fm3 / self.epsilon_cut_mev_fm3**self.gamma_s,
            "cutoff": {
                "epsilon_mev_fm3": self.epsilon_cut_mev_fm3,
                "pressure_mev_fm3": self.pressure_cut_mev_fm3,
                "number_density_fm3": self.number_density_cut_fm3,
                "chemical_potential_mev": self.chemical_potential_cut_mev,
                "cs2": self.cs2_cut,
            },
            "tail": {
                "min_energy_ratio": min_energy_ratio,
                "min_pressure_mev_fm3": float(state["pressure_mev_fm3"][0]),
                "min_epsilon_mev_fm3": float(state["epsilon_mev_fm3"][0]),
                "min_number_density_fm3": float(state["number_density_fm3"][0]),
                "min_chemical_potential_mev": float(state["chemical_potential_mev"][0]),
                "cs2_min": float(np.min(state["cs2"])),
                "cs2_max": float(np.max(state["cs2"])),
            },
            "closure": {
                "max_abs_mu_n_minus_eps_plus_p_rel": float(np.max(np.abs(
                    state["chemical_potential_mev"] * state["number_density_fm3"] /
                    (state["epsilon_mev_fm3"] + state["pressure_mev_fm3"]) - 1.0
                ))),
                "max_abs_first_law_dlnn_dlnE_rel": float(np.max(np.abs(
                    first_law_lhs / first_law_rhs - 1.0
                ))),
            },
        }


def build_extended_projected_rows(
    rows: Sequence[Mapping[str, object]],
    *,
    alpha: float = 0.0,
    family: str = "power_law",
    tail_points: int = 512,
    min_energy_ratio: float = 1.0e-12,
) -> list[dict[str, str]]:
    """Prepend a positive-state low-density tail; preserve all source rows verbatim."""
    if isinstance(tail_points, bool) or not isinstance(tail_points, int) or tail_points < 16:
        raise ValueError("tail_points must be an integer >= 16")
    model = ThermodynamicSurfaceContinuation.from_projected_rows(rows, alpha=alpha, family=family)
    model._validate_shape(min_energy_ratio=min_energy_ratio)
    x = np.geomspace(min_energy_ratio, 1.0, tail_points + 1)[:-1]
    state = model.state_many(model.epsilon_cut_mev_fm3 * x)
    aux = model.auxiliary_cutoff
    mn = model.neutron_mass_mev
    output: list[dict[str, str]] = []
    for i in range(tail_points):
        eps = float(state["epsilon_mev_fm3"][i])
        pressure = float(state["pressure_mev_fm3"][i])
        nb = float(state["number_density_fm3"][i])
        mu = float(state["chemical_potential_mev"][i])
        row = {
            "source_row": str(i - tail_points),
            "i_t": aux["i_t"],
            "i_nb": str(i - tail_points),
            "i_yq": aux["i_yq"],
            "n_b_fm3": format(nb, ".17g"),
            "n_b_m3": format(nb * FM3_TO_M3, ".17g"),
            "pressure_mev_fm3": format(pressure, ".17g"),
            "pressure_pa": format(pressure * MEV_FM3_TO_PA, ".17g"),
            "energy_mev_fm3": format(eps, ".17g"),
            "energy_j_m3": format(eps * MEV_FM3_TO_PA, ".17g"),
            "mu_b_mev": format(mu, ".17g"),
            "mu_b_j": format(mu * MEV_TO_J, ".17g"),
            "q1_p_over_nb_mev": format(pressure / nb, ".17g"),
            "q2_entropy_per_baryon": aux["q2_entropy_per_baryon"],
            "q3_mu_b_over_mn_minus_1": format(mu / mn - 1.0, ".17g"),
            "q4_mu_q_over_mn": aux["q4_mu_q_over_mn"],
            "q5_mu_l_over_mn": aux["q5_mu_l_over_mn"],
            "q6_f_over_nb_mn_minus_1": format(eps / nb / mn - 1.0, ".17g"),
            "q7_e_over_nb_mn_minus_1": format(eps / nb / mn - 1.0, ".17g"),
        }
        output.append(row)
    for source in rows:
        output.append({field: str(source[field]) for field in NORMALIZED_FIELDS})
    # The native reader preserves input order and requires all primary fields
    # to be positive; also ensure the added rows do not disturb source order.
    primary = ("n_b_m3", "pressure_pa", "energy_j_m3", "mu_b_j")
    for field in primary:
        values = np.asarray([float(row[field]) for row in output])
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0) or np.any(np.diff(values) <= 0.0):
            raise ValueError(f"extended normalized table is not positive and monotone in {field}")
    return output


def write_extended_projected_csv(
    source_path: str | Path,
    output_path: str | Path,
    *,
    alpha: float = 0.0,
    family: str = "power_law",
    tail_points: int = 512,
    min_energy_ratio: float = 1.0e-12,
) -> dict:
    source = Path(source_path)
    with source.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != NORMALIZED_FIELDS:
            raise ValueError("projected EOS header does not match the frozen 19-column normalized table")
        rows = list(reader)
    model = ThermodynamicSurfaceContinuation.from_projected_rows(rows, alpha=alpha, family=family)
    extended = build_extended_projected_rows(
        rows, alpha=alpha, family=family, tail_points=tail_points, min_energy_ratio=min_energy_ratio
    )
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=NORMALIZED_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(extended)
    import hashlib
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    output_digest = hashlib.sha256(target.read_bytes()).hexdigest()
    return {
        "source_path": str(source),
        "source_sha256": digest,
        "output_path": str(target),
        "output_sha256": output_digest,
        "source_row_count": len(rows),
        "output_row_count": len(extended),
        "tail_row_count": tail_points,
        "surface_pressure_pa": float(extended[0]["pressure_pa"]),
        "surface_energy_j_m3": float(extended[0]["energy_j_m3"]),
        "surface_number_density_m3": float(extended[0]["n_b_m3"]),
        "diagnostics": model.diagnostics(min_energy_ratio=min_energy_ratio, points=max(512, tail_points)),
    }
