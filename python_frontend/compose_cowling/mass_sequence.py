"""Initial rising-M(Pc) proxy, kinematic segments and paired responses."""
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
    """Keep the initial connected rising-M(Pc) segment as a branch proxy.

    This is NOT a radial-stability classification. A local mass maximum can
    change the number of unstable radial modes, and later stability recovery
    is EOS-dependent. A maximum is accepted after two consecutive fractional
    decreases larger than the tolerance. Failed samples are hard barriers.
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

    candidate_index, run = None, 0
    for i in range(1, len(points)):
        m0, m1 = float(points[i-1]["mass_kg"]), float(points[i]["mass_kg"])
        decreased = (m0 - m1) / m0 > fractional_decrease_tolerance
        if decreased:
            if candidate_index is None:
                candidate_index = i - 1
            run += 1
            if run >= sustained_decreases:
                branch = tuple(points[:candidate_index+1])
                peak = max(branch, key=lambda p: float(p["mass_kg"]))
                end = branch[-1]
                return StableBranch(branch, True, float(peak["mass_kg"]), float(peak["pc_pa"]),
                    float(end["mass_kg"]), float(end["pc_pa"]), "first_turning_point_resolved", truncated_low)
        elif candidate_index is not None:
            # A significant decline followed by recovery is not silently
            # bridged. Stop before the candidate and require grid refinement.
            branch = tuple(points[:candidate_index+1])
            end = branch[-1]
            return StableBranch(branch, False, None, None, float(end["mass_kg"]),
                float(end["pc_pa"]), "turning_point_candidate_unconfirmed", truncated_low)
        else:
            run = 0

    if candidate_index is not None:
        branch = tuple(points[:candidate_index+1])
        end = branch[-1]
        return StableBranch(branch, False, None, None, float(end["mass_kg"]),
            float(end["pc_pa"]), "turning_point_candidate_unconfirmed", truncated_low)

    end = points[-1]
    return StableBranch(tuple(points), False, None, None, float(end["mass_kg"]),
        float(end["pc_pa"]), "sequence_gap_truncated" if gap else "rising_branch_truncated",
        truncated_low)



def extract_candidate_mass_segments(sequence: Iterable[Mapping], *,
    fractional_slope_tolerance: float = 1.0e-4) -> list[dict]:
    """Partition each contiguous TOV scan block into monotonic M(Pc) segments.

    This is a kinematic diagnostic only. It does NOT classify radial stability:
    every returned segment is tagged unresolved until a radial pulsation
    eigenvalue/mode-count calculation is available. Turning points are shared
    endpoints of adjacent segments. Failed samples remain hard barriers.
    """
    if (not math.isfinite(fractional_slope_tolerance)
            or fractional_slope_tolerance < 0.0):
        raise ValueError("fractional_slope_tolerance must be finite and nonnegative")
    ordered = sorted((dict(p) for p in sequence),
                     key=lambda p: int(p.get("scan_index", 0)))
    blocks: list[list[dict]] = []
    current: list[dict] = []
    previous_scan = None
    for point in ordered:
        if not _valid_point(point):
            if current:
                blocks.append(current)
                current = []
            previous_scan = None
            continue
        scan_index = int(point["scan_index"])
        pc = float(point["pc_pa"])
        if current:
            if previous_scan is not None and scan_index <= previous_scan:
                raise ValueError("scan_index must be strictly increasing within a TOV block")
            if pc <= float(current[-1]["pc_pa"]):
                raise ValueError("central pressures must be strictly increasing within a TOV block")
            # A skipped index denotes a missing sample even if it was omitted
            # from the input iterable; do not bridge it silently.
            if previous_scan is not None and scan_index != previous_scan + 1:
                blocks.append(current)
                current = []
        current.append(point)
        previous_scan = scan_index
    if current:
        blocks.append(current)

    output: list[dict] = []
    for block_id, points in enumerate(blocks):
        if not points:
            continue
        if len(points) == 1:
            p = points[0]
            output.append({
                "sequence_block_id": block_id,
                "segment_in_block": 0,
                "direction": "unresolved_single_sample",
                "start_scan_index": int(p["scan_index"]),
                "end_scan_index": int(p["scan_index"]),
                "start_pc_pa": float(p["pc_pa"]),
                "end_pc_pa": float(p["pc_pa"]),
                "start_mass_kg": float(p["mass_kg"]),
                "end_mass_kg": float(p["mass_kg"]),
                "start_radius_m": float(p["radius_m"]),
                "end_radius_m": float(p["radius_m"]),
                "sample_count": 1,
                "start_turning_point_kind": "sequence_start",
                "end_turning_point_kind": "sequence_end_or_gap",
                "radial_stability_status": "unresolved_requires_radial_mode_spectrum",
            })
            continue

        raw_signs: list[int] = []
        for left, right in zip(points[:-1], points[1:]):
            m0, m1 = float(left["mass_kg"]), float(right["mass_kg"])
            fractional_slope = (m1 - m0) / max(abs(m0), abs(m1), 1.0e-300)
            if fractional_slope > fractional_slope_tolerance:
                raw_signs.append(1)
            elif fractional_slope < -fractional_slope_tolerance:
                raw_signs.append(-1)
            else:
                raw_signs.append(0)

        nonzero = [i for i, sign in enumerate(raw_signs) if sign != 0]
        if not nonzero:
            signs = [0] * len(raw_signs)
        else:
            # Treat a flat plateau as part of the preceding run. At a
            # maximum/minimum plateau this places the shared turning point at
            # the last plateau sample and makes the convention deterministic.
            signs = list(raw_signs)
            previous_sign = None
            for i, sign in enumerate(signs):
                if sign != 0:
                    previous_sign = sign
                elif previous_sign is not None:
                    signs[i] = previous_sign
            next_sign = None
            for i in range(len(signs)-1, -1, -1):
                if signs[i] != 0:
                    next_sign = signs[i]
                elif next_sign is not None:
                    signs[i] = next_sign

        start_node = 0
        current_sign = signs[0]
        start_kind = "sequence_start"
        segment_index = 0
        for edge_index in range(1, len(signs)):
            next_sign = signs[edge_index]
            if next_sign == current_sign:
                continue
            if current_sign == 0 or next_sign == 0:
                # A flat mass interval has no well-defined increasing/decreasing
                # side; keep it in a separate unresolved segment.
                boundary_kind = "flat_or_degenerate_turning_candidate"
            elif current_sign > 0 and next_sign < 0:
                boundary_kind = "local_maximum_candidate"
            elif current_sign < 0 and next_sign > 0:
                boundary_kind = "local_minimum_candidate"
            else:
                boundary_kind = "unresolved_direction_change"
            end_node = edge_index
            segment_points = points[start_node:end_node+1]
            first, last = segment_points[0], segment_points[-1]
            output.append({
                "sequence_block_id": block_id,
                "segment_in_block": segment_index,
                "direction": ("rising" if current_sign > 0 else
                              "falling" if current_sign < 0 else "flat"),
                "start_scan_index": int(first["scan_index"]),
                "end_scan_index": int(last["scan_index"]),
                "start_pc_pa": float(first["pc_pa"]),
                "end_pc_pa": float(last["pc_pa"]),
                "start_mass_kg": float(first["mass_kg"]),
                "end_mass_kg": float(last["mass_kg"]),
                "start_radius_m": float(first["radius_m"]),
                "end_radius_m": float(last["radius_m"]),
                "sample_count": len(segment_points),
                "start_turning_point_kind": start_kind,
                "end_turning_point_kind": boundary_kind,
                "radial_stability_status": "unresolved_requires_radial_mode_spectrum",
            })
            segment_index += 1
            start_node = end_node
            start_kind = boundary_kind
            current_sign = next_sign

        segment_points = points[start_node:]
        first, last = segment_points[0], segment_points[-1]
        output.append({
            "sequence_block_id": block_id,
            "segment_in_block": segment_index,
            "direction": ("rising" if current_sign > 0 else
                          "falling" if current_sign < 0 else "flat"),
            "start_scan_index": int(first["scan_index"]),
            "end_scan_index": int(last["scan_index"]),
            "start_pc_pa": float(first["pc_pa"]),
            "end_pc_pa": float(last["pc_pa"]),
            "start_mass_kg": float(first["mass_kg"]),
            "end_mass_kg": float(last["mass_kg"]),
            "start_radius_m": float(first["radius_m"]),
            "end_radius_m": float(last["radius_m"]),
            "sample_count": len(segment_points),
            "start_turning_point_kind": start_kind,
            "end_turning_point_kind": "sequence_end_or_gap",
            "radial_stability_status": "unresolved_requires_radial_mode_spectrum",
        })
    for global_index, row in enumerate(output):
        row["candidate_segment_id"] = f"B{row['sequence_block_id']:02d}-S{row['segment_in_block']:02d}"
    return output


def find_mass_bracket(branch: StableBranch, target_mass_kg: float) -> dict | None:
    """Find adjacent samples on the initial rising-M(Pc) proxy segment."""
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
