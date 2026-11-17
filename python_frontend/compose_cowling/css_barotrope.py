"""Exact constant-sound-speed bridge followed by a separate core interpolant.

Physical units are MeV fm^-3, fm^-3 and MeV.  The CSS law is
epsilon(P)=epsilon_t+(P-P_t)/s and
n(P)=n_t*((epsilon+P)/(epsilon_t+P_t))**(1/(1+s)).
The core starts at the continuous endpoint; its derivative is one-sided.
"""
from __future__ import annotations

import numpy as np

from .direct_barotrope import DirectBarotrope, EOSState


class CSSCoreBarotrope:
    """Preserve the analytic bridge and the core's own pressure-energy knots."""

    def __init__(self, pressure, energy, number_density, chemical_potential,
                 *, join_index: int, sound_speed_sq: float):
        arrays = [np.asarray(v, dtype=float) for v in
                  (pressure, energy, number_density, chemical_potential)]
        p, e, n, mu = arrays
        if (any(a.ndim != 1 for a in arrays) or len(p) < 4
                or any(len(a) != len(p) for a in arrays)
                or any(not np.all(np.isfinite(a)) or np.any(a <= 0) for a in arrays)):
            raise ValueError("CSS/core arrays must be finite positive equal-length vectors")
        if any(np.any(np.diff(a) <= 0) for a in arrays):
            raise ValueError("CSS/core thermodynamic knots must increase strictly")
        if (isinstance(join_index, (bool, np.bool_))
                or not isinstance(join_index, (int, np.integer))
                or not 1 <= join_index < len(p) - 1):
            raise ValueError("invalid CSS/core join index")
        s = float(sound_speed_sq)
        if not np.isfinite(s) or not 0 < s <= 1:
            raise ValueError("CSS sound speed must satisfy 0 < cs2 <= 1")
        self.pressure, self.energy = p.copy(), e.copy()
        self.number_density_source = n.copy()
        self.join_index = int(join_index)
        self.join_pressure = float(p[join_index])
        self.sound_speed_sq = s
        self._p0, self._e0, self._n0 = map(float, (p[0], e[0], n[0]))
        closure = np.abs(e + p - mu * n) / np.abs(e + p)
        self.raw_first_law_max_relative_residual = float(np.max(closure))
        if self.raw_first_law_max_relative_residual > 2e-12:
            raise ValueError("CSS/core source Euler identity exceeds 2e-12 relative tolerance")
        exact = [self._css_state(float(value)) for value in p[:join_index + 1]]
        predicted = (np.array([v.energy for v in exact]),
                     np.array([v.number_density for v in exact]),
                     np.array([v.chemical_potential for v in exact]))
        for actual, expected in zip((e, n, mu), predicted):
            if np.any(np.abs(actual[:join_index + 1] - expected)
                      > 5e-9 * np.maximum(np.abs(actual[:join_index + 1]), np.abs(expected))):
                raise ValueError("CSS report disagrees with the serialized analytic bridge")
        # Source endpoint gates are those declared by the native CSS constructor.
        endpoint = exact[-1]
        for actual, expected in ((e[join_index], endpoint.energy),
                                 (n[join_index], endpoint.number_density),
                                 (mu[join_index], endpoint.chemical_potential)):
            if abs(actual - expected) > (2e-10 * max(1., abs(expected))
                                        + 64 * np.finfo(float).eps * max(1., abs(actual), abs(expected))):
                raise ValueError("CSS/core endpoint continuity contract failed")
        self.core = DirectBarotrope(p[join_index:], e[join_index:],
                                   n[join_index:], mu[join_index:])
        self.number_density_knots = np.r_[
            [v.number_density for v in exact[:-1]], self.core.number_density_knots]
        self.mu_knots = np.r_[
            [v.chemical_potential for v in exact[:-1]], self.core.mu_knots]

    def _css_state(self, p: float) -> EOSState:
        eps = self._e0 + (p - self._p0) / self.sound_speed_sq
        h0 = self._e0 + self._p0
        h = eps + p
        if not np.isfinite(h) or h <= 0:
            raise ValueError("CSS analytic bridge overflows its thermodynamic domain")
        try:
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                n = self._n0 * np.exp(np.log1p((h - h0) / h0) / (1 + self.sound_speed_sq))
                mu = h / n
        except FloatingPointError as exc:
            raise ValueError("CSS analytic bridge overflows its thermodynamic domain") from exc
        if not np.all(np.isfinite([eps, n, mu])) or min(eps, n, mu) <= 0:
            raise ValueError("CSS analytic bridge has non-finite or nonpositive state")
        return EOSState(p, eps, float(n), float(mu), self.sound_speed_sq)

    def evaluate(self, pressure: float) -> EOSState:
        p = float(pressure)
        if not np.isfinite(p) or not self.pressure[0] <= p <= self.pressure[-1]:
            raise ValueError("pressure query outside original CSS/core table domain")
        if p < self.join_pressure:
            return self._css_state(p)
        return self.core.evaluate(p)
