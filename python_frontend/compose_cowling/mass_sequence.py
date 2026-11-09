"""Stable-branch extraction and balanced paired-response diagnostics."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence
import math


@dataclass(frozen=True)
class StableBranch:
    points: tuple[dict, ...]
    turning_point_found: bool
    maximum_mass_kg: float | None
    maximum_mass_pc_pa: float | None
    branch_endpoint_mass_kg: float | None
    branch_endpoint_pc_pa: float | None
    status: str
    low_pressure_truncated: bool = False


def _valid_point(point: Mapping) -> bool:
    try:
        return (
            str(point.get("status", "ok")).lower() == "ok"
            and math.isfinite(float(point["pc_pa"]))
            and math.isfinite(float(point["mass_kg"]))
            and math.isfinite(float(point["radius_m"]))
            and float(point["pc_pa"]) > 0
            and float(point["mass_kg"]) > 0
            and float(point["radius_m"]) > 0
        )
    except (KeyError, TypeError, ValueError):
        return False


def extract_stable_branch(sequence: Iterable[Mapping], *,
    fractional_decrease_tolerance: float = 1.0e-4,
    sustained_decreases: int = 2) -> StableBranch:
    """Keep the initial connected rising-mass branch up to its first maximum.

    A maximum is accepted after two consecutive fractional decreases larger
    than the tolerance. A failed sample terminates the connected sequence;
    later samples are never bridged. An unresolved endpoint is not called Mmax.
    """
    if not math.isfinite(fractional_decrease_tolerance) or fractional_decrease_tolerance < 0:
        raise ValueError("fractional_decrease_tolerance must be finite and nonnegative")
    if sustained_decreases < 1:
        raise ValueError("sustained_decreases must be at least one")
    ordered = sorted((dict(p) for p in sequence), key=lambda p: int(p.get("scan_index", 0)))
    first = next((i for i, p in enumerate(ordered) if _valid_point(p)), None)
    if first is None:
        return StableBranch((), False, None, None, None, None, "no_valid_tov_samples", True)

    points, gap, previous_pc = [], False, None
    for point in ordered[first:]:
        if not _valid_point(point):
            gap = True
            break
        pc = float(point["pc_pa"])
        if previous_pc is not None and pc <= previous_pc:
            raise ValueError("central pressures must be strictly increasing")
        points.append(point)
        previous_pc = pc

    truncated_low = first > 0 or gap
    if len(points) < 2:
        end = points[-1] if points else None
        return StableBranch(tuple(points), False, None, None,
            float(end["mass_kg"]) if end else None, float(end["pc_pa"]) if end else None,
            "insufficient_connected_tov_samples", truncated_low)

    peak_index, run = None, 0
    for i in range(1, len(points)):
        m0, m1 = float(points[i-1]["mass_kg"]), float(points[i]["mass_kg"])
        if (m0 - m1) / m0 > fractional_decrease_tolerance:
            run += 1
        else:
            run = 0
        if run >= sustained_decreases:
            peak_index = i - sustained_decreases
            break

    if peak_index is not None:
        branch = tuple(points[:peak_index+1])
        peak = max(branch, key=lambda p: float(p["mass_kg"]))
        end = branch[-1]
        return StableBranch(branch, True, float(peak["mass_kg"]), float(peak["pc_pa"]),
            float(end["mass_kg"]), float(end["pc_pa"]), "first_turning_point_resolved", truncated_low)

    end = points[-1]
    return StableBranch(tuple(points), False, None, None, float(end["mass_kg"]),
        float(end["pc_pa"]), "sequence_gap_truncated" if gap else "rising_branch_truncated",
        truncated_low)


def find_mass_bracket(branch: StableBranch, target_mass_kg: float) -> dict | None:
    """Find adjacent stable-branch samples; final Pc is solved by Fortran."""
    if not math.isfinite(target_mass_kg) or target_mass_kg <= 0:
        raise ValueError("target_mass_kg must be finite and positive")
    for left, right in zip(branch.points[:-1], branch.points[1:]):
        m0, m1 = float(left["mass_kg"]), float(right["mass_kg"])
        if m1 <= m0 or not (m0 <= target_mass_kg <= m1):
            continue
        t = min(1.0, max(0.0, (target_mass_kg-m0)/(m1-m0)))
        lp = math.log(float(left["pc_pa"])) + t * (
            math.log(float(right["pc_pa"])) - math.log(float(left["pc_pa"])))
        lr = math.log(float(left["radius_m"])) + t * (
            math.log(float(right["radius_m"])) - math.log(float(left["radius_m"])))
        return {"pc_lo_pa": float(left["pc_pa"]), "pc_hi_pa": float(right["pc_pa"]),
            "pc_seed_pa": math.exp(lp), "radius_guess_m": math.exp(lr),
            "mass_lo_kg": m0, "mass_hi_kg": m1, "mass_fraction": t}
    return None


def _positive_finite(value) -> bool:
    try:
        x = float(value)
        return math.isfinite(x) and x > 0
    except (TypeError, ValueError):
        return False


def paired_sensitivity(rows: Sequence[Mapping], expected_methods: Sequence[str]
    ) -> tuple[list[dict], list[dict]]:
    """Compute six-method paired deviations; incomplete pairs are never imputed."""
    methods = tuple(expected_methods)
    if len(methods) != 6 or len(set(methods)) != 6:
        raise ValueError("the crossed matching design requires six unique methods")
    groups: dict[tuple, list[dict]] = {}
    for source in rows:
        row = dict(source)
        key = (str(row["crust"]), str(row["core"]), float(row["mass_msun"]))
        groups.setdefault(key, []).append(row)

    paired, deltas = [], []
    for (crust, core, mass), members in sorted(groups.items()):
        by_method, duplicate = {}, False
        for row in members:
            method = str(row.get("method", ""))
            if method in by_method:
                duplicate = True
            by_method[method] = row
        valid = {m: by_method[m] for m in methods if m in by_method
            and bool(by_method[m].get("valid", False))
            and _positive_finite(by_method[m].get("f_hz"))
            and _positive_finite(by_method[m].get("p1_hz"))}
        complete = not duplicate and len(valid) == 6
        record = {"crust": crust, "core": core, "mass_msun": mass,
            "pair_complete": complete, "n_valid_methods": len(valid),
            "missing_methods": ";".join(m for m in methods if m not in valid),
            "duplicate_method_rows": duplicate, "f_mean_hz": None, "p1_mean_hz": None,
            "f_range_percent": None, "p1_range_percent": None}
        if complete:
            fv = [float(valid[m]["f_hz"]) for m in methods]
            pv = [float(valid[m]["p1_hz"]) for m in methods]
            fm, pm = sum(fv)/6, sum(pv)/6
            record.update({"f_mean_hz": fm, "p1_mean_hz": pm,
                "f_range_percent": 100*(max(fv)-min(fv))/fm,
                "p1_range_percent": 100*(max(pv)-min(pv))/pm})
            for method in methods:
                row = valid[method]
                deltas.append({"crust": crust, "core": core, "mass_msun": mass, "method": method,
                    "f_hz": float(row["f_hz"]), "p1_hz": float(row["p1_hz"]),
                    "delta_f_percent": 100*(float(row["f_hz"])/fm-1),
                    "delta_p1_percent": 100*(float(row["p1_hz"])/pm-1)})
        paired.append(record)
    return paired, deltas
