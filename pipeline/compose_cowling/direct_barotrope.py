"""Direct source-knot EOS provider: monotone PCHIP g(log P), first-law n reconstruction.

The original pressure and energy-density knots are never modified.  Baryon
number density and chemical potential are derived from the cold first law.
"""
from __future__ import annotations
import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.interpolate import PchipInterpolator

_GL_X, _GL_W = np.polynomial.legendre.leggauss(16)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


@dataclass(frozen=True)
class EOSState:
    pressure: float
    energy: float
    number_density: float
    chemical_potential: float
    sound_speed_sq: float


class DirectBarotrope:
    """A causal, thermodynamically closed continuous curve through raw P-epsilon knots."""

    def __init__(self, pressure: np.ndarray, energy: np.ndarray, number_density: np.ndarray, chemical_potential: np.ndarray):
        p = np.asarray(pressure, dtype=float)
        e = np.asarray(energy, dtype=float)
        n = np.asarray(number_density, dtype=float)
        mu_source = np.asarray(chemical_potential, dtype=float)
        if any(a.ndim != 1 for a in (p, e, n, mu_source)) or len(p) < 2 or not (len(p) == len(e) == len(n) == len(mu_source)):
            raise ValueError("P, epsilon and n_B must be equally sized one-dimensional arrays (N>=2)")
        if any(not np.all(np.isfinite(a)) for a in (p, e, n, mu_source)):
            raise ValueError("EOS arrays contain non-finite values")
        if np.any(p <= 0) or np.any(e <= 0) or np.any(n <= 0) or np.any(mu_source <= 0):
            raise ValueError("P, epsilon and n_B must be positive")
        if np.any(np.diff(p) <= 0) or np.any(np.diff(e) <= 0) or np.any(np.diff(n) <= 0):
            raise ValueError("P, epsilon and n_B must be strictly increasing")
        closure = np.abs(e + p - mu_source * n) / np.maximum(np.abs(e + p), 1e-300)
        self.raw_first_law_max_relative_residual = float(np.max(closure))
        if self.raw_first_law_max_relative_residual > 2e-12:
            raise ValueError("source knot first-law identity exceeds 2e-12 relative tolerance")
        if np.any(np.diff(mu_source) <= 0):
            raise ValueError("source chemical potential must be strictly increasing")
        g = e - p
        if np.any(g <= 0) or np.any(np.diff(g) <= 0):
            raise ValueError("g=epsilon-P must be positive and strictly increasing")
        self.pressure = p.copy()
        self.energy = e.copy()
        self.number_density_source = n.copy()
        self.u = np.log(p)
        self.g = g
        self._g_of_u = PchipInterpolator(self.u, g, extrapolate=False)
        self._dg_du = self._g_of_u.derivative()
        self.logn_knots = np.empty_like(p)
        self.logn_knots[0] = np.log(n[0])
        for i in range(len(p) - 1):
            self.logn_knots[i + 1] = self.logn_knots[i] + self._integrate_cell(i, self.u[i], self.u[i + 1])
        self.number_density_knots = np.exp(self.logn_knots)
        self.mu_knots = (e + p) / self.number_density_knots
        if not np.all(np.isfinite(self.mu_knots)) or np.any(self.mu_knots <= 0) or np.any(np.diff(self.mu_knots) <= 0):
            raise ValueError("first-law chemical potential is not finite, positive and strictly increasing")

    def _integrand(self, u: np.ndarray | float) -> np.ndarray:
        p = np.exp(u)
        g = self._g_of_u(u)
        gu = self._dg_du(u)
        return (p + gu) / (2.0 * p + g)

    def _integrate_cell(self, i: int, ua: float, ub: float) -> float:
        if ub == ua:
            return 0.0
        mid = 0.5 * (ua + ub)
        half = 0.5 * (ub - ua)
        points = mid + half * _GL_X
        return float(half * np.dot(_GL_W, self._integrand(points)))

    def evaluate(self, pressure: float) -> EOSState:
        p = float(pressure)
        if not np.isfinite(p) or p < self.pressure[0] or p > self.pressure[-1]:
            raise ValueError("pressure query outside original table domain")
        u = float(np.log(p))
        i = min(int(np.searchsorted(self.pressure, p, side="left")) - 1, len(self.pressure) - 2)
        i = max(i, 0)
        g = float(self._g_of_u(u))
        gu = float(self._dg_du(u))
        if gu < 0 and gu > -1e-12 * max(abs(g), 1.0):
            gu = 0.0
        eps = p + g
        logn = self.logn_knots[i] + self._integrate_cell(i, self.u[i], u)
        n = float(np.exp(logn))
        mu = (eps + p) / n
        chi = p / (p + gu)
        if not np.all(np.isfinite([eps, n, mu, chi])) or min(eps, n, mu, chi) <= 0 or chi > 1 + 2e-12:
            raise ArithmeticError("interpolated EOS violates positive causal thermodynamics")
        return EOSState(p, eps, n, mu, min(chi, 1.0))

    def evaluate_mu(self, chemical_potential: float) -> EOSState:
        mu = float(chemical_potential)
        if not np.isfinite(mu) or mu < self.mu_knots[0] or mu > self.mu_knots[-1]:
            raise ValueError("chemical-potential query outside original table domain")
        lo, hi = float(self.pressure[0]), float(self.pressure[-1])
        for _ in range(80):
            mid = lo + 0.5 * (hi - lo)
            if self.evaluate(mid).chemical_potential < mu:
                lo = mid
            else:
                hi = mid
            if hi - lo <= 8 * np.finfo(float).eps * max(abs(mid), np.finfo(float).tiny):
                break
        return self.evaluate(lo + 0.5 * (hi - lo))

    def audit(self) -> dict[str, Any]:
        chi_samples = []
        max_first_law = 0.0
        max_n_rel = float(np.max(np.abs(self.number_density_knots / self.number_density_source - 1.0)))
        for i in range(len(self.pressure) - 1):
            u = np.linspace(self.u[i], self.u[i + 1], 65)
            p = np.exp(u)
            gu = self._dg_du(u)
            chi = p / (p + gu)
            chi_samples.extend(chi.tolist())
        for p, e, n, mu in zip(self.pressure, self.energy, self.number_density_knots, self.mu_knots):
            max_first_law = max(max_first_law, abs(e + p - mu * n) / max(abs(e + p), 1e-300))
        return {
            "rows": int(len(self.pressure)),
            "intervals": int(len(self.pressure) - 1),
            "pressure_knots_preserved": True,
            "energy_knots_preserved": True,
            "source_number_density_max_relative_change": max_n_rel,
            "first_law_max_relative_residual": float(max_first_law),
            "raw_knot_first_law_max_relative_residual": self.raw_first_law_max_relative_residual,
            "min_chi_sampled": float(np.min(chi_samples)),
            "max_chi_sampled": float(np.max(chi_samples)),
            "g_strictly_increasing": bool(np.all(np.diff(self.g) > 0)),
            "quadrature": "16-point Gauss-Legendre per native log-pressure cell",
        }


def audit_csv(source: Path) -> dict[str, Any]:
    source = Path(source)
    raw = np.genfromtxt(source, delimiter=",", names=True, dtype=None, encoding="utf-8")
    model = DirectBarotrope(
        raw["pressure_mev_fm3"], raw["energy_mev_fm3"], raw["n_b_fm3"], raw["mu_b_mev"]
    )
    return {
        "schema": "direct-source-knot-p-epsilon-pchip-v1",
        "source_path": str(source.resolve()),
        "source_sha256": sha256(source),
        "representation": "PCHIP(g=epsilon-P versus ln(P)); n_B from cold first-law quadrature; mu_B=(epsilon+P)/n_B",
        "source_mutated": False,
        "metrics": model.audit(),
    }


def audit_directory(source_dir: Path, output_dir: Path) -> dict[str, Any]:
    source_dir, output_dir = Path(source_dir), Path(output_dir)
    files = [p for p in sorted(source_dir.glob("*.csv")) if p.name != "spectrum_summary_234.csv"]
    if not files:
        raise FileNotFoundError(f"no source EOS CSV files found in {source_dir}")
    receipts = []
    for source in files:
        receipt = audit_csv(source)
        target = output_dir / f"{source.stem}.direct_eos.json"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        receipts.append(receipt)
    manifest = {
        "schema": "direct-source-knot-eos-audit-manifest-v1",
        "n_tables": len(receipts),
        "source_tables_modified": 0,
        "provider": "monotone_PCHIP_epsilon_minus_pressure_of_logP_with_firstlaw_density_v1",
        "all_admissible": all(
            r["metrics"]["min_chi_sampled"] > 0
            and r["metrics"]["max_chi_sampled"] <= 1 + 2e-12
            and r["metrics"]["first_law_max_relative_residual"] < 5e-14
            for r in receipts
        ),
        "receipts": receipts,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "direct_eos_audit_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest
