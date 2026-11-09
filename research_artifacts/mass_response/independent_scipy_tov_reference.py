#!/usr/bin/env python3
"""Independent SciPy TOV reference using the native-knot cubic-Hermite P(mu) EOS.

Each raw table interval is a cubic Hermite pressure P(mu), with endpoint
derivatives equal to the raw baryon densities. The TOV ODE is integrated
independently in radius with SciPy DOP853, then compared with the frozen native
Fortran cell_solve TOV backend. No Cowling mode or radial-stability claim is
made by this experiment.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
MEV_FM3_TO_PA = 1.602176634e32
G_SI = 6.67430e-11
C_SI = 299792458.0
EPS = np.finfo(float).eps


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


@dataclass(frozen=True)
class NativeHermiteEOS:
    """Piecewise cubic P(mu) with exact raw-knot n=dP/dmu endpoints."""

    mu: np.ndarray
    pressure: np.ndarray
    density: np.ndarray
    source_sha256: str
    path: Path

    @classmethod
    def from_csv(cls, path: Path) -> "NativeHermiteEOS":
        rows = read_rows(path)
        mu = np.asarray([float(row["mu_b_j"]) for row in rows], dtype=float)
        pressure = np.asarray([float(row["pressure_pa"]) for row in rows], dtype=float)
        density = np.asarray([float(row["n_b_m3"]) for row in rows], dtype=float)
        if len(mu) < 3 or not (len(mu) == len(pressure) == len(density)):
            raise ValueError(f"invalid EOS row count: {path}")
        if not (np.all(np.isfinite(mu)) and np.all(np.isfinite(pressure))
                and np.all(np.isfinite(density))):
            raise ValueError(f"nonfinite EOS data: {path}")
        if (np.any(np.diff(mu) <= 0.0) or np.any(np.diff(pressure) <= 0.0)
                or np.any(density <= 0.0)):
            raise ValueError(f"EOS knots are not strictly admissible: {path}")
        return cls(mu, pressure, density, sha256(path), path)

    def _cell(self, index: int) -> tuple[float, float, float, float, float, float]:
        mu0, mu1 = self.mu[index:index + 2]
        p0, p1 = self.pressure[index:index + 2]
        n0, n1 = self.density[index:index + 2]
        h = mu1 - mu0
        a1 = h * n0
        a2 = 3.0 * (p1 - p0) - h * (2.0 * n0 + n1)
        a3 = 2.0 * (p0 - p1) + h * (n0 + n1)
        return float(mu0), float(h), float(p0), float(a1), float(a2), float(a3)

    def state_at_pressure(self, pressure: float) -> tuple[float, float, float, float]:
        """Return (epsilon, n, mu, cs2) for a pressure in the table domain."""
        if not math.isfinite(pressure) or pressure < self.pressure[0] or pressure > self.pressure[-1]:
            raise ValueError(f"pressure outside EOS domain: {pressure}")
        index = int(np.searchsorted(self.pressure, pressure, side="right") - 1)
        index = max(0, min(index, len(self.pressure) - 2))
        mu0, h, p0, a1, a2, a3 = self._cell(index)
        p1 = float(self.pressure[index + 1])
        if pressure <= p0:
            t = 0.0
        elif pressure >= p1:
            t = 1.0
        else:
            lo, hi = 0.0, 1.0
            dp = pressure - p0
            t = min(1.0, max(0.0, dp / (p1 - p0)))
            tolerance = 16.0 * EPS * max(abs(pressure), abs(p0), 1.0)
            for _ in range(80):
                f = t * (a1 + t * (a2 + t * a3)) - dp
                if abs(f) <= tolerance:
                    break
                if f > 0.0:
                    hi = t
                else:
                    lo = t
                derivative = a1 + t * (2.0 * a2 + 3.0 * t * a3)
                next_t = t - f / derivative
                if not (lo < next_t < hi) or not math.isfinite(next_t):
                    next_t = 0.5 * (lo + hi)
                if next_t == t:
                    break
                t = next_t
        mu = mu0 + h * t
        p = p0 + t * (a1 + t * (a2 + t * a3))
        n = (a1 + t * (2.0 * a2 + 3.0 * t * a3)) / h
        dn_dmu = (2.0 * a2 + 6.0 * a3 * t) / (h * h)
        epsilon = mu * n - p
        cs2 = n / (mu * dn_dmu)
        if not np.all(np.isfinite([epsilon, n, mu, cs2])) or min(epsilon, n, mu, cs2) <= 0.0:
            raise ArithmeticError("native Hermite EOS interpolation returned an inadmissible state")
        if cs2 > 1.0 + 64.0 * EPS:
            raise ArithmeticError(f"native Hermite EOS interpolation is acausal: cs2={cs2}")
        return float(epsilon), float(n), float(mu), float(cs2)


def scipy_tov(eos: NativeHermiteEOS, pc_pa: float, *, rtol: float = 2.0e-10,
              max_step_m: float = 2000.0, r0_m: float = 10.0,
              radius_limit_m: float = 1.0e8) -> dict:
    """Direct-r integration of the SI TOV equations with a finite-pressure event."""
    p_cut = float(eos.pressure[0])
    if not (p_cut < pc_pa < eos.pressure[-1]):
        raise ValueError("central pressure lies outside the finite EOS integration domain")
    eps_c = eos.state_at_pressure(pc_pa)[0]
    alpha = (2.0 * math.pi * G_SI / 3.0) * (
        (eps_c + pc_pa) / C_SI**2
    ) * ((eps_c + 3.0 * pc_pa) / C_SI**2)
    pressure_drop = alpha * r0_m**2
    if not (0.0 < pressure_drop < 0.01 * pc_pa):
        raise ValueError("centre expansion does not resolve a sub-percent pressure drop")
    m0 = (4.0 * math.pi / 3.0) * (eps_c / C_SI**2) * r0_m**3
    p0 = pc_pa - pressure_drop

    def rhs(radius: float, state: np.ndarray) -> np.ndarray:
        mass, pressure = map(float, state)
        pressure_for_eos = max(pressure, p_cut)
        epsilon = eos.state_at_pressure(pressure_for_eos)[0]
        compactness = 2.0 * G_SI * mass / (radius * C_SI**2)
        if not math.isfinite(compactness) or compactness >= 1.0:
            raise ArithmeticError("TOV compactness reached a trapped-surface denominator")
        phi = (G_SI / C_SI**2) * (
            mass / radius**2 + 4.0 * math.pi * radius * pressure_for_eos / C_SI**2
        ) / (1.0 - compactness)
        return np.asarray([
            4.0 * math.pi * radius**2 * epsilon / C_SI**2,
            -(epsilon + pressure_for_eos) * phi,
        ])

    def surface_event(radius: float, state: np.ndarray) -> float:
        return float(state[1] - p_cut)

    surface_event.terminal = True
    surface_event.direction = -1.0
    atol = np.asarray([1.0e10, max(1.0, p_cut * 1.0e-12)])
    solution = solve_ivp(
        rhs, (r0_m, radius_limit_m), np.asarray([m0, p0]),
        method="DOP853", rtol=rtol, atol=atol, max_step=max_step_m,
        events=surface_event,
    )
    if not solution.success:
        raise RuntimeError(f"SciPy DOP853 TOV integration failed: {solution.message}")
    if len(solution.t_events[0]) != 1:
        raise RuntimeError("SciPy TOV integration did not reach the declared finite-pressure surface")
    radius = float(solution.t_events[0][0])
    mass = float(solution.y_events[0][0][0])
    return {
        "mass_kg": mass,
        "radius_m": radius,
        "surface_pressure_pa": p_cut,
        "surface_residual_pa": float(solution.y_events[0][0][1] - p_cut),
        "accepted_steps": int(solution.t.size - 1),
        "rhs_evaluations": int(solution.nfev),
        "rtol": rtol,
        "max_step_m": max_step_m,
    }


def build_reference(output_dir: Path, *, relative_tolerance: float = 1.0e-6) -> dict:
    import sys
    sys.path.insert(0, str(ROOT / "python_frontend"))
    from compose_cowling.solver import solve_v5_tov_star

    eos_dir = ROOT / "oscillation_release_v1" / "data" / "eos_projected"
    case_specs = [
        ("css", 0.22679842504532394),
        ("css", 0.280428),
        ("maxwell", 1.1164034103761742),
        ("maxwell", 1.897632203650715),
        ("css", 10.0),
    ]
    results = []
    for method, pc_mev_fm3 in case_specs:
        path = eos_dir / f"GPPVA_FSU2__BBB_BHF2__{method}.csv"
        eos = NativeHermiteEOS.from_csv(path)
        pc_pa = pc_mev_fm3 * MEV_FM3_TO_PA
        independent = scipy_tov(eos, pc_pa)
        native = solve_v5_tov_star(
            path, pc_pa, p_cut_pa=float(eos.pressure[0]),
            r0_m=10.0, r_limit_m=1.0e8, rtol=1.0e-9, step_cap=0.05,
        )
        dm = abs(independent["mass_kg"] - float(native["mass_kg"])) / max(abs(float(native["mass_kg"])), 1.0)
        dr = abs(independent["radius_m"] - float(native["radius_m"])) / max(abs(float(native["radius_m"])), 1.0)
        results.append({
            "method": method,
            "pc_mev_fm3": pc_mev_fm3,
            "pc_pa": pc_pa,
            "eos_sha256": eos.source_sha256,
            "interpolation_id": "cubic_Hermite_P_of_mu_raw_n_v1",
            "scipy": independent,
            "native_fortran": {
                "mass_kg": float(native["mass_kg"]),
                "radius_m": float(native["radius_m"]),
                "surface_residual_pa": float(native.get("surface_residual_pa", math.nan)),
                "backend_build_key": native.get("backend_build_key"),
            },
            "relative_mass_difference": dm,
            "relative_radius_difference": dr,
        })
    max_mass = max(row["relative_mass_difference"] for row in results)
    max_radius = max(row["relative_radius_difference"] for row in results)
    if max(max_mass, max_radius) > relative_tolerance:
        raise AssertionError(
            f"SciPy/native TOV disagreement exceeds {relative_tolerance:g}: "
            f"mass={max_mass:g}, radius={max_radius:g}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    receipt = {
        "schema": "ns-cowling-independent-scipy-tov-reference-v1",
        "status": "PASS_SAME_EOS_INTERPOLATION_INDEPENDENT_RADIAL_ODE",
        "relative_tolerance": relative_tolerance,
        "max_relative_mass_difference": max_mass,
        "max_relative_radius_difference": max_radius,
        "cases": results,
        "scope": (
            "DOP853 direct-r SI TOV integration compared with frozen v5 Fortran "
            "cell_solve task=star. Both use the identical raw-knot cubic-Hermite "
            "P(mu) interpolation and finite positive surface-pressure cutoff. "
            "This is a structural TOV cross-check only, not a radial-stability "
            "or f/p1-spectrum validation."
        ),
    }
    (output_dir / "independent_scipy_tov_reference.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("dist/independent-scipy-tov"))
    parser.add_argument("--relative-tolerance", type=float, default=1.0e-6)
    args = parser.parse_args(argv)
    receipt = build_reference(args.output_dir, relative_tolerance=args.relative_tolerance)
    print(json.dumps({
        "status": receipt["status"],
        "case_count": len(receipt["cases"]),
        "max_relative_mass_difference": receipt["max_relative_mass_difference"],
        "max_relative_radius_difference": receipt["max_relative_radius_difference"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
